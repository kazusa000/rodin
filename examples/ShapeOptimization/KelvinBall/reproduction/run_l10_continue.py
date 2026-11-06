#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
"""Recover the rebooted L10 experiment on 28h and continue to global update 100."""
import argparse
import json
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import sys
import time
import traceback

from run_l10_n50 import UPSTREAM, PARAMETERS, THREADS, complete_history, git, pinned, sha, verify_blas, write
from restart_snapshot import export_restart
from recovery_validation import compare_fields, compare_native

BASE = 'b9d84fea5998f388ee903cf92f99c009c5a890fa'
OFFSET = 24
FINAL = 100


def audit(source):
    name = 'examples/ShapeOptimization/KelvinBall/KelvinBallOptimization.cpp'
    old = git(source, 'show', f'{BASE}:{name}')
    new = (source / name).read_text()
    marker = '    const Real h = meanElementSize(mesh);\n    const auto meshDiagnostics'
    if old.split(marker, 1)[1].strip() != new.split(marker, 1)[1].strip():
        raise ValueError('Numerical update-loop body changed beyond recovery/checkpoint I/O')
    if git(source, 'diff', BASE, 'HEAD', '--', 'src'):
        raise ValueError('Framework source changed')
    return dict(base=BASE, numerical_loop_after_checkpoint='unchanged', framework_src='unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'build', 'output', 'prefix', 'parent', 'blas-library-dir'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--commit', required=True)
    parser.add_argument('--compiled-campaign', type=Path)
    parser.add_argument('--framework-campaign', type=Path)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    if platform.node() != '28h':
        raise RuntimeError('All compilation, tests and numerical work require 28h')
    pinned(args.source, args.commit)
    for path in ((args.output,) if args.compiled_campaign else (args.build, args.output)):
        if path.exists() or path.resolve().is_relative_to(args.source.resolve()):
            raise ValueError('Require fresh external build/output paths')
    parent = json.loads((args.parent / 'campaign.json').read_text())
    interruption = json.loads((args.parent / 'reboot-interruption.json').read_text())
    if parent['commit'] != BASE or not interruption:
        raise ValueError('Unexpected parent experiment or missing reboot evidence')
    original = args.parent / 'optimization'
    rows = complete_history(original / 'kelvin-ball.csv', 101)
    if len(rows) != 24 or int(rows[-1]['iteration']) != 23:
        raise ValueError('Require the preserved parent states 0..23')
    target = float(rows[0]['target_volume'])
    if any(float(row['target_volume']) != target for row in rows):
        raise ValueError('Original target volume changed')
    snapshot = original / f'KelvinBall.Chamber.mesh.{OFFSET:06d}.h5'
    args.output.mkdir(parents=True)
    env = dict(os.environ, OMP_NUM_THREADS='12', OPENBLAS_NUM_THREADS='12', MKL_NUM_THREADS='12',
               KELVIN_SYMMETRY='cubic', PYTHONDONTWRITEBYTECODE='1')
    lib = args.prefix / 'lib/x86_64-linux-gnu'
    env['LD_LIBRARY_PATH'] = ':'.join(map(str, (args.blas_library_dir, lib, lib/'openmpi/lib', lib/'hdf5/serial')))
    record = dict(status='preparing', started_at=time.time(), pid=os.getpid(), host=platform.node(),
                  source=str(args.source), commit=args.commit, upstream=UPSTREAM, parent=str(args.parent),
                  parent_csv_sha256=sha(original/'kelvin-ball.csv'), saved_snapshot=str(snapshot),
                  saved_snapshot_sha256=sha(snapshot), original_target_volume=target,
                  parent_fully_evaluated_states=24, parent_completed_updates=23,
                  recovered_global_step=OFFSET, requested_final_global_step=FINAL,
                  requested_additional_updates=FINAL-OFFSET, continuation_states=FINAL-OFFSET+1,
                  recovery_policy='original fitted chamber and binary coordinates; independent SWIFT background reinitialized as a copy of the preserved chamber',
                  missing_original_swift_background=True, exact_old_swift_trajectory=False,
                  paired_native_checkpoints=True, recurring_monitor=False, independent_bem=False,
                  native_wall_time_limit_seconds=None, automatic_restart=False, stages=[],
                  source_audit=audit(args.source), fem_threads=THREADS, blas_threads=THREADS)

    def checkpoint():
        write(args.output/'campaign.json', record)

    def stop(signum, frame):
        raise KeyboardInterrupt(f'Received signal {signum}')
    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)

    def execute(name, command, budget=None):
        pinned(args.source, args.commit)
        directory = args.output/name
        directory.mkdir(parents=True)
        item = dict(name=name, command=list(map(str, command)), started_at=time.time(),
                    timeout_seconds=budget, status='running')
        record['stages'].append(item)
        record['status'] = name
        checkpoint()
        interrupted = None
        with (directory/'run.log').open('w') as log:
            process = subprocess.Popen(item['command'], cwd=directory, env=env,
                                       stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            item['pid'] = process.pid
            checkpoint()
            try:
                code = process.wait(timeout=budget)
            except BaseException as error:
                interrupted = error
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                code = process.returncode
                item['interrupted_by'] = type(error).__name__
        item.update(exit_code=code, finished_at=time.time(), status='completed' if code == 0 else ('expected_rejection' if name == 'validation/reject-invalid-tag' and code < 0 else 'failed'))
        (directory/'exit-code.txt').write_text(str(code)+'\n')
        checkpoint()
        if interrupted is not None:
            raise interrupted
        return code

    code = 1
    try:
        checkpoint()
        if execute('validation/recovery-guards', [sys.executable, '-B', '-m', 'unittest', 'discover',
                    '-s', args.source/'examples/ShapeOptimization/KelvinBall/reproduction',
                    '-p', 'test_recovery_validation.py', '-v'], 120):
            raise RuntimeError('Recovery guard tests failed')
        prep = args.output/'prepare'
        prep.mkdir()
        record['recovery_export'] = export_restart(snapshot, prep/'recovered.mesh', 10)
        record['source_audit']['diff_sha256'] = __import__('hashlib').sha256(
            git(args.source, 'diff', BASE, args.commit, '--', 'src', 'examples').encode()).hexdigest()
        checkpoint()
        mpi = lib/'openmpi'
        configure = ['cmake', '-S', args.source/'examples/ShapeOptimization/KelvinBall/reproduction',
                     '-B', args.build, '-DCMAKE_BUILD_TYPE=Release', f'-DCMAKE_PREFIX_PATH={args.prefix}',
                     f'-DMUMPS_ROOT={args.prefix}', f'-DMPI_C_COMPILER_INCLUDE_DIRS={mpi}/include;{mpi}/include/openmpi',
                     f'-DMPI_C_HEADER_DIR={mpi}/include', '-DMPI_C_LIB_NAMES=mpi', f'-DMPI_mpi_LIBRARY={mpi}/lib/libmpi.so']
        if args.framework_campaign:
            command = [sys.executable, '-B', args.source/'examples/ShapeOptimization/KelvinBall/reproduction/build_recovery_from_framework.py',
                       '--source', args.source, '--build', args.build, '--campaign', args.framework_campaign]
            if execute('compile', command, 4*3600):
                raise RuntimeError('Recovery compilation against preserved framework failed')
            record['reused_framework_provenance'] = json.loads((args.build/'reused-framework-provenance.json').read_text())
        elif args.compiled_campaign:
            compiled = json.loads(args.compiled_campaign.read_text())
            pinned(Path(compiled['source']), compiled['commit'])
            changed = git(args.source, 'diff', '--name-only', compiled['commit'], args.commit, '--', 'src', 'examples')
            if any(not path.startswith('examples/ShapeOptimization/KelvinBall/reproduction/')
                   for path in changed.splitlines()):
                raise RuntimeError('Reused binary numerical source differs')
            for name in ('configure', 'compile'):
                stage = next(item for item in compiled['stages'] if item['name'] == name)
                actual = int((args.compiled_campaign.parent/name/'exit-code.txt').read_text())
                if actual != 0 or stage['exit_code'] != 0:
                    raise RuntimeError('Reused build lacks successful exits')
            if sha(args.build/'bin/KelvinBall') != compiled['binary_sha256']:
                raise RuntimeError('Reused binary hash differs')
            record['compiled_build_provenance'] = dict(campaign=str(args.compiled_campaign),
                source=compiled['source'], commit=compiled['commit'], binary_sha256=compiled['binary_sha256'])
        else:
            for name, command in [('configure', configure), ('compile', ['cmake', '--build', args.build,
                                   '--target', 'KelvinBall', 'KelvinCheckpointTest', '-j', '4'])]:
                if execute(name, command, 4*3600):
                    raise RuntimeError(f'{name} failed')
        binary, test = args.build/'bin/KelvinBall', args.build/'bin/KelvinCheckpointTest'
        record['binary_sha256'] = sha(binary)
        record['math_libraries'] = verify_blas(binary, env, args.blas_library_dir)
        if execute('validation/checkpoint', [test], 120):
            raise RuntimeError('Synthetic checkpoint I/O failed')
        fixture = args.output/'validation/checkpoint'
        count = compare_native(fixture/'fixture-original.h5', fixture/'fixture-restored.h5')
        import h5py
        import numpy as np
        import shutil
        bad = fixture/'invalid-tags.h5'
        shutil.copy2(fixture/'fixture-original.h5', bad)
        with h5py.File(bad, 'r+') as data:
            del data['KelvinMMG/Tags0']
            data['KelvinMMG'].create_dataset('Tags0', data=np.array([2**63], dtype=np.uint64))
        if execute('validation/reject-invalid-tag', [test, '--load', bad], 120) == 0:
            raise RuntimeError('Out-of-range MMG tags were accepted')
        record['checkpoint_validation'] = dict(all_datasets_exact=True, datasets=count, invalid_tag_rejected=True)
        recovery = [f'--restart-mesh={prep}/recovered.mesh', f'--restart-coordinates={snapshot}',
                    f'--target-volume={target:.17g}', f'--global-offset={OFFSET}', '--save-checkpoint']
        def parameters(states):
            return [f'--iterations={states}' if p.startswith('--iterations=') else p for p in PARAMETERS] + recovery
        if execute('validation/recovery', [binary, *parameters(1)], 2*3600):
            raise RuntimeError('Saved-state reevaluation failed')
        validation = args.output/'validation/recovery'
        fields = compare_fields(original, validation, OFFSET)
        recovered = complete_history(validation/'kelvin-ball.csv', 1)
        if len(recovered) != 1 or float(recovered[0]['target_volume']) != target:
            raise RuntimeError('Recovered score or original target volume missing')
        record['recovery_validation'] = dict(status='passed', fields=fields,
                         fem={key:float(recovered[0][key]) for key in ('k','c','q','rho','volume','target_volume')},
                         global_step=OFFSET, original_full_fields_checked=True)
        checkpoint()
        if sha(snapshot) != record['saved_snapshot_sha256'] or sha(original/'kelvin-ball.csv') != record['parent_csv_sha256']:
            raise RuntimeError('Original parent files changed')
        native = execute('optimization', [binary, *parameters(FINAL-OFFSET+1)])
        history = complete_history(args.output/'optimization/kelvin-ball.csv', FINAL-OFFSET+1)
        record.update(native_exit_code=native, fully_evaluated_states=len(history),
                      completed_additional_updates=max(0,len(history)-1),
                      last_solved_global_step=OFFSET+max(0,len(history)-1),
                      requested_budget_completed=native==0 and len(history)==FINAL-OFFSET+1)
        if history:
            record['final_fem'] = {key:float(history[-1][key]) for key in ('k','c','q','rho','volume','target_volume')}
        record['status'] = 'completed' if record['requested_budget_completed'] else 'native_ended_early'
        code = 0 if record['requested_budget_completed'] else 1
    except BaseException as error:
        record.update(status='interrupted' if isinstance(error, KeyboardInterrupt) else 'failed',
                      error=str(error), error_type=type(error).__name__)
        traceback.print_exc()
    finally:
        record.update(exit_code=code, finished_at=time.time())
        checkpoint()
        (args.output/'exit-code.txt').write_text(str(code)+'\n')
    return code


if __name__ == '__main__':
    sys.exit(main())
