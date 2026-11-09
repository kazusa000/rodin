#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
"""Build and test MUMPS disk storage, then replay one saved L10 reconstruction on 28h."""
import argparse
import json
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import time
import traceback

from run_l10_n50 import git, pinned, sha, verify_blas, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'build', 'output', 'prefix', 'blas-library-dir',
                 'checkpoint', 'advected', 'ooc-directory'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--commit', required=True)
    parser.add_argument('--framework-build', type=Path)
    parser.add_argument('--framework-source', type=Path)
    parser.add_argument('--memory-mb', type=int, default=16000)
    args = parser.parse_args()
    if platform.node() != '28h':
        raise RuntimeError('Compilation, tests and reconstruction require 28h')
    if args.memory_mb <= 0:
        raise ValueError('Require a positive MUMPS workspace limit')
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    pinned(args.source, args.commit)
    for path in (args.build, args.output, args.ooc_directory):
        if path.exists() or path.resolve().is_relative_to(args.source.resolve()):
            raise ValueError('Require fresh external build, output and factor directories')
    for path in (args.checkpoint/'READY', args.checkpoint/'background.h5', args.advected):
        if not path.is_file():
            raise ValueError('Missing completed checkpoint or saved transported field')
    args.output.mkdir(parents=True)
    args.ooc_directory.mkdir(parents=True)
    env = dict(os.environ, OMP_NUM_THREADS='12', OPENBLAS_NUM_THREADS='12',
               MKL_NUM_THREADS='12', KELVIN_SYMMETRY='cubic',
               PYTHONDONTWRITEBYTECODE='1')
    lib = args.prefix/'lib/x86_64-linux-gnu'
    env['LD_LIBRARY_PATH'] = ':'.join(map(str, (args.blas_library_dir, lib,
                                              lib/'openmpi/lib', lib/'hdf5/serial')))
    inputs = {str(path): sha(path) for path in
              (args.checkpoint/'background.h5', args.checkpoint/'state.h5', args.advected)}
    record = dict(status='preparing', host=platform.node(), pid=os.getpid(),
                  source=str(args.source), commit=args.commit, started_at=time.time(),
                  checkpoint=str(args.checkpoint), input_sha256=inputs,
                  memory_mb=args.memory_mb, ooc_directory=str(args.ooc_directory),
                  global_step=26, optimization_updates=0, grid_parameters_changed=False,
                  automatic_continuation=False, recurring_monitor=False, stages=[])

    def checkpoint():
        write(args.output/'campaign.json', record)

    def execute(name, command, budget=4*3600):
        pinned(args.source, args.commit)
        directory = args.output/name
        directory.mkdir()
        item = dict(name=name, command=list(map(str, command)), status='running',
                    started_at=time.time(), timeout_seconds=budget)
        record['stages'].append(item)
        record['status'] = name
        checkpoint()
        with (directory/'run.log').open('w') as log, (directory/'resources.jsonl').open('w') as samples:
            process = subprocess.Popen(item['command'], cwd=directory, env=env,
                                       stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
            item['pid'] = process.pid
            checkpoint()
            try:
                while process.poll() is None:
                    now = time.time()
                    status = {}
                    try:
                        for line in Path(f'/proc/{process.pid}/status').read_text().splitlines():
                            if line.startswith(('VmRSS:', 'VmHWM:', 'VmSwap:')):
                                key, value = line.split(':', 1)
                                status[key] = int(value.split()[0])
                    except FileNotFoundError:
                        pass
                    item['peak_rss_kib'] = max(item.get('peak_rss_kib', 0), status.get('VmHWM', 0))
                    factors = list(args.ooc_directory.iterdir())
                    size = 0
                    for path in factors:
                        try:
                            size += path.stat().st_size
                        except FileNotFoundError:
                            pass
                    item['peak_factor_file_bytes'] = max(item.get('peak_factor_file_bytes', 0), size)
                    samples.write(json.dumps(dict(time=now, **status,
                                                  factor_files=len(factors), factor_bytes=size))+'\n')
                    samples.flush()
                    if now-item['started_at'] > budget:
                        item['interrupted_by'] = 'timeout'
                        os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=15)
                    except subprocess.TimeoutExpired:
                        if item.get('interrupted_by'):
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait()
            except BaseException:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=15)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                raise
        item.update(exit_code=process.returncode, finished_at=time.time(),
                    status='completed' if process.returncode == 0 else 'failed')
        (directory/'exit-code.txt').write_text(str(process.returncode)+'\n')
        checkpoint()
        if process.returncode:
            raise RuntimeError(f'{name} failed with actual exit {process.returncode}')

    code = 1
    try:
        mpi = lib/'openmpi'
        configure_command = ['cmake', '-S', args.source/'examples/ShapeOptimization/KelvinBall/reproduction',
                             '-B', args.build, '-DCMAKE_BUILD_TYPE=Release',
                             '-DKELVIN_BUILD_MEMORY_TESTS=ON',
                             f'-DCMAKE_PREFIX_PATH={args.prefix}', f'-DMUMPS_ROOT={args.prefix}',
                             f'-DMPI_C_COMPILER_INCLUDE_DIRS={mpi}/include;{mpi}/include/openmpi',
                             f'-DMPI_C_HEADER_DIR={mpi}/include', '-DMPI_C_LIB_NAMES=mpi',
                             f'-DMPI_mpi_LIBRARY={mpi}/lib/libmpi.so']
        if args.framework_build:
            if not args.framework_source:
                raise ValueError('Immutable framework reuse requires its source')
            execute('compile', ['python3', args.source/'examples/ShapeOptimization/KelvinBall/reproduction/build_ooc_from_framework.py',
                                '--source', args.source, '--build', args.build,
                                '--framework-build', args.framework_build,
                                '--framework-source', args.framework_source])
        else:
            execute('configure', configure_command)
            execute('compile', ['cmake', '--build', args.build, '--target',
                                'KelvinBall', 'KelvinMemoryTest', '-j', '2'])
        binary = args.build/'bin/KelvinBall'
        record['binary_sha256'] = sha(binary)
        record['math_libraries'] = verify_blas(binary, env, args.blas_library_dir)
        execute('validation', [args.build/'bin/KelvinMemoryTest'], 3600)
        command = [binary, '--n=50', '--outer-radius=10', '--hmin-factor=0.05',
                   '--hmax-factor=10', '--regularization=20', '--step=0.2',
                   '--reconstruction=swift', '--mmg-adapt', '--thickness-min=0.05',
                   '--swift-linear-solver=mumps', '--swift-linear-threads=12',
                   '--swift-model-fit=10', '--swift-convergence-iterations-outer=30',
                   '--swift-convergence-iterations-inner=15',
                   '--swift-linear-out-of-core=1', f'--swift-linear-memory-mb={args.memory_mb}',
                   f'--swift-linear-ooc-directory={args.ooc_directory}',
                   '--swift-linear-ooc-prefix=l10-step26',
                   f'--replay-background={args.checkpoint/"background.h5"}',
                   f'--replay-advected={args.advected}']
        execute('reconstruction', command)
        log = (args.output/'reconstruction/run.log').read_text()
        if 'out_of_core=1  memory_mb=' + str(args.memory_mb) not in log:
            raise RuntimeError('No successful MUMPS disk solve was reported')
        if 'metric MUMPS: ok=0' in log:
            raise RuntimeError('A MUMPS factorization or solve failed')
        if record['stages'][-1].get('peak_factor_file_bytes', 0) == 0:
            raise RuntimeError('No numeric factor data was observed on disk')
        if not (args.output/'reconstruction/RECONSTRUCTION_COMPLETE').is_file():
            raise RuntimeError('No completed reconstruction marker')
        if any(sha(Path(path)) != digest for path, digest in inputs.items()):
            raise RuntimeError('Original replay inputs changed')
        if list(args.ooc_directory.iterdir()):
            raise RuntimeError('Factor files remained after a normal solver shutdown')
        record['status'] = 'completed_single_reconstruction'
        code = 0
    except BaseException as error:
        record.update(status='failed', error=f'{type(error).__name__}: {error}', traceback=traceback.format_exc())
    finally:
        record['finished_at'] = time.time()
        checkpoint()
        (args.output/'exit-code.txt').write_text(str(code)+'\n')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
