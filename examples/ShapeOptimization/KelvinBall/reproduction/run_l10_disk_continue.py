#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
"""Continue the completed global26 SWIFT reconstruction as states27..100 on28h."""
import argparse
import csv
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import sys
import time
import traceback
import json

from run_l10_n50 import PARAMETERS, complete_history, git, pinned, sha, verify_blas, write

BASE = 'a6b3240f33dd01043019d85ae76a000a808f86e1'
FIRST, FINAL = 27, 100


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'build', 'output', 'prefix', 'blas-library-dir',
                 'probe', 'parent', 'ooc-directory', 'framework-build', 'framework-source'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--commit', required=True)
    args = parser.parse_args()
    if platform.node() != '28h':
        raise RuntimeError('Compilation, validation and optimization require28h')
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    pinned(args.source, args.commit)
    for path in (args.build, args.output, args.ooc_directory):
        if path.exists() or path.resolve().is_relative_to(args.source.resolve()):
            raise ValueError('Require fresh external build/output/factor directories')
    if git(args.source, 'rev-parse', 'HEAD:src') != git(args.source, 'rev-parse', BASE+':src'):
        raise ValueError('Framework source differs from the tested disk-storage version')
    entry = 'examples/ShapeOptimization/KelvinBall/KelvinBallOptimization.cpp'
    marker = '    const Real h = meanElementSize(mesh);\n    const auto meshDiagnostics'
    old = git(args.source, 'show', BASE+':'+entry)
    current = (args.source/entry).read_text()
    if old.split(marker,1)[1] != current.split(marker,1)[1]:
        raise ValueError('Numerical optimization loop changed')
    probe = json.loads((args.probe/'campaign.json').read_text())
    if probe['status'] != 'completed_single_reconstruction':
        raise ValueError('The saved reconstruction did not complete')
    for path in (args.probe/'exit-code.txt', args.probe/'reconstruction/exit-code.txt'):
        if int(path.read_text()) != 0:
            raise ValueError('Require successful native/controller reconstruction exits')
    if not (args.probe/'reconstruction/RECONSTRUCTION_COMPLETE').is_file():
        raise ValueError('Reconstruction completion marker is missing')
    chamber = args.probe/'reconstruction/replayed-chamber.h5'
    background = args.probe/'reconstruction/replayed-background.h5'
    inputs = {str(chamber):sha(chamber), str(background):sha(background)}
    expected = ('bcabef1b2b45ade6f9211b8b0850d202cde724fc6930746e651fed41d151d174',
                '6fa21db12c76ca42581865533e5723a02ffae05cf2a438394a45c64b1e31942c')
    if tuple(inputs.values()) != expected:
        raise ValueError('The approved completed native mesh pair changed')
    parent_rows = complete_history(args.parent/'optimization/kelvin-ball.csv', 101)
    target = float(parent_rows[0]['target_volume'])
    if target != 0.1744360039783153:
        raise ValueError('Original target volume differs')
    args.output.mkdir(parents=True)
    args.ooc_directory.mkdir(parents=True)
    env = dict(os.environ, OMP_NUM_THREADS='12', OPENBLAS_NUM_THREADS='12',
               MKL_NUM_THREADS='12', KELVIN_SYMMETRY='cubic', PYTHONDONTWRITEBYTECODE='1')
    lib = args.prefix/'lib/x86_64-linux-gnu'
    env['LD_LIBRARY_PATH'] = ':'.join(map(str,(args.blas_library_dir,lib,lib/'openmpi/lib',lib/'hdf5/serial')))
    record = dict(status='preparing',started_at=time.time(),host=platform.node(),pid=os.getpid(),
                  source=str(args.source),commit=args.commit,origin='https://github.com/kazusa000/rodin.git',
                  branch='codex/fem-l10-n50-20261105',build=str(args.build),
                  restart_probe=str(args.probe),restart_input_sha256=inputs,
                  first_global_state=FIRST,final_global_state=FINAL,requested_states=FINAL-FIRST+1,
                  requested_new_updates_after_loaded_shape=FINAL-FIRST,
                  previous_complete_csv_global_state=25,previous_logged_fem_global_state=26,
                  starting_shape='completed global26 reconstruction, geometry for global27',
                  original_target_volume=target,initial_history_predictions='unavailable after process restart; diagnostic only',
                  memory_mumps_mb=16000,ooc_directory=str(args.ooc_directory),grid_parameters_changed=False,
                  optimization_loop_unchanged=True,initial_sphere_reset=False,
                  recurring_monitor_requested=True,independent_bem=False,automatic_restart=False,
                  native_wall_time_limit_seconds=None,stages=[])

    def checkpoint():
        write(args.output/'campaign.json',record)

    def interrupt(signum, frame):
        raise KeyboardInterrupt(f'Received signal{signum}')
    signal.signal(signal.SIGINT,interrupt)
    signal.signal(signal.SIGTERM,interrupt)

    def execute(name,command,budget=None,expected_rejection=False):
        pinned(args.source,args.commit)
        directory=args.output/name
        directory.mkdir(parents=True)
        item=dict(name=name,command=list(map(str,command)),status='running',started_at=time.time(),timeout_seconds=budget)
        record['stages'].append(item);record['status']=name;checkpoint()
        interrupted=None
        with (directory/'run.log').open('w') as log:
            process=subprocess.Popen(item['command'],cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            item['pid']=process.pid;checkpoint()
            try:
                code=process.wait(timeout=budget)
            except BaseException as error:
                interrupted=error
                if process.poll() is None:
                    os.killpg(process.pid,signal.SIGTERM)
                    try:process.wait(timeout=15)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid,signal.SIGKILL);process.wait()
                code=process.returncode;item['interrupted_by']=type(error).__name__
        accepted=(code!=0) if expected_rejection else (code==0)
        item.update(exit_code=code,finished_at=time.time(),status='expected_rejection' if expected_rejection and accepted else ('completed' if accepted else 'failed'))
        (directory/'exit-code.txt').write_text(str(code)+'\n');checkpoint()
        if interrupted is not None:raise interrupted
        if not accepted:raise RuntimeError(f'{name} failed with actual exit{code}')
        return code

    code=1
    try:
        execute('compile',[sys.executable,args.source/'examples/ShapeOptimization/KelvinBall/reproduction/build_ooc_from_framework.py',
                           '--source',args.source,'--build',args.build,'--framework-build',args.framework_build,
                           '--framework-source',args.framework_source],4*3600)
        binary=args.build/'bin/KelvinBall'
        record['binary_sha256']=sha(binary)
        record['math_libraries']=verify_blas(binary,env,args.blas_library_dir)
        test_filter=('Rodin_Solver_MUMPS.*:AllGeometries/*:Rodin_Adaptation_SWIFTExampleParameters.*:'
                     'Rodin_Adaptation_SWIFTReport.*:SWIFTLinearSolverTest.MumpsDiskBackendPreservesCenteredOperator')
        execute('validation/mumps',[args.build/'bin/KelvinMemoryTest','--gtest_filter='+test_filter],3600)
        parameters=[f'--iterations={FINAL-FIRST+1}' if p.startswith('--iterations=') else p for p in PARAMETERS]
        parameters += [f'--restart-chamber={chamber}',f'--restart-background={background}',
                       f'--target-volume={target:.17g}',f'--global-offset={FIRST}','--save-checkpoint',
                       '--swift-linear-out-of-core=1','--swift-linear-memory-mb=16000',
                       f'--swift-linear-ooc-directory={args.ooc_directory}','--swift-linear-ooc-prefix=l10-continue']
        record['native_parameters']=parameters;checkpoint()
        execute('validation/native-pair',[binary,*parameters,'--geometry-only'],600)
        if 'Native restart: preserved chamber and independent SWIFT background loaded directly' not in (args.output/'validation/native-pair/run.log').read_text():
            raise RuntimeError('Native paired restart was not used')
        execute('validation/reject-incomplete-pair',[binary,'--restart-chamber='+str(chamber)],60,True)
        execute('optimization',[binary,*parameters])
        rows=complete_history(args.output/'optimization/kelvin-ball.csv',FINAL-FIRST+1)
        record.update(fully_evaluated_states=len(rows),completed_updates_after_loaded_shape=max(0,len(rows)-1),
                      last_complete_global_state=FIRST+len(rows)-1 if rows else None)
        if len(rows)!=FINAL-FIRST+1:raise RuntimeError('The requested global100 target is incomplete')
        if any(sha(Path(path))!=digest for path,digest in inputs.items()):
            raise RuntimeError('Original restart pair changed')
        record['status']='completed_through_global100';code=0
    except BaseException as error:
        record.update(status='failed',error=f'{type(error).__name__}: {error}',traceback=traceback.format_exc())
        path=args.output/'optimization/kelvin-ball.csv'
        if path.exists():
            try:
                rows=complete_history(path,FINAL-FIRST+1)
                record.update(fully_evaluated_states=len(rows),last_complete_global_state=FIRST+len(rows)-1 if rows else None)
            except Exception as csv_error:record['csv_error']=str(csv_error)
    finally:
        record['finished_at']=time.time();checkpoint()
        (args.output/'exit-code.txt').write_text(str(code)+'\n')
    return code


if __name__=='__main__':
    raise SystemExit(main())
