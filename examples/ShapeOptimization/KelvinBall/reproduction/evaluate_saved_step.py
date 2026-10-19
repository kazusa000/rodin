#!/usr/bin/env python3
"""Evaluate one saved fem-old Cubic state without starting optimization."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

from run_fem_old import (EXPORTER, VERIFIER, accepted_scan, complete_rows,
                         digest, pinned, save)

FEM_COMMIT = 'fed5689233fc035e99f9feb1fe4bad8d9e681eef'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['input', 'output', 'verifier', 'exporter', 'prior-bem']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--commit', required=True)
    parser.add_argument('--step', type=int, required=True)
    parser.add_argument('--input-sha256', required=True)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if platform.node() != '28h':
        raise ValueError('Computation is authorized only on 28h')
    controller = Path(__file__).resolve().parents[4]
    pinned(controller, args.commit)
    pinned(args.verifier, VERIFIER)
    pinned(args.exporter, EXPORTER)
    manifest = json.loads((args.input/'campaign.json').read_text())
    if manifest.get('commit') != FEM_COMMIT or manifest.get('symmetry') != 'cubic':
        raise ValueError('Wrong saved FEM experiment')
    rows = complete_rows(args.input/'optimization/kelvin-ball.csv')
    if args.step < 0 or args.step >= len(rows):
        raise ValueError('Requested step has no complete FEM state')
    snapshot = args.input/'optimization'/f'KelvinBall.Chamber.mesh.{args.step:06d}.h5'
    if digest(snapshot) != args.input_sha256:
        raise ValueError('Saved HDF5 differs from the selected geometry')
    if args.output.exists() or args.output.resolve().is_relative_to(controller):
        raise ValueError('Require a fresh output outside source')
    args.output.mkdir(parents=True)
    record = dict(state='exporting', controller_commit=args.commit,
        source_commit=FEM_COMMIT, verifier_commit=VERIFIER, exporter_commit=EXPORTER,
        input_snapshot=str(snapshot), input_sha256=args.input_sha256,
        selected_bem_step=args.step, selected_bem_geometry_fem_evaluated=True,
        prior_bem=str(args.prior_bem), solver_rtol=1e-6,
        optimization_started=False, started_at=time.time())
    checkpoint = lambda: save(args.output/'campaign.json', record)
    code = 1
    proc = None

    def interrupted(signum, _frame):
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                record['solver_exit_code'] = proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                record['solver_exit_code'] = proc.wait()
            (args.output/'solver/exit-code.txt').write_text(str(record['solver_exit_code'])+'\n')
        raise InterruptedError(f'Stopped by signal {signum}')

    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    try:
        checkpoint()
        sys.path.insert(0, str(args.exporter/'scripts'))
        os.environ['KELVIN_SYMMETRY'] = 'cubic'
        spec = importlib.util.spec_from_file_location('saved_export',
            args.exporter/'scripts/reproduce_smallbox_original.py')
        exporter = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(exporter)
        native = args.output/'native-export.mesh'
        mesh = args.output/'original-double.mesh'
        exporter.export_snapshot(snapshot, native)
        record['export_certificate'] = exporter.original_double_mesh(snapshot, native, mesh)
        if (record['export_certificate'].get('geometry_modified') is not False
                or record['export_certificate'].get('topology_unchanged') is not True
                or record['export_certificate'].get('gamma_panels') != int(rows[args.step]['interface_triangles'])):
            raise ValueError('Original shape/topology/panel count did not match the saved state')
        record['mesh_sha256'] = digest(mesh)
        if digest(snapshot) != args.input_sha256:
            raise ValueError('Original HDF5 changed during export')
        if args.prepare_only:
            record['state'] = 'prepared_only'
            code = 0
        else:
            record['state'] = 'waiting_for_prior_terra_bem'
            checkpoint()
            deadline = time.monotonic()+74*3600
            while not (args.prior_bem/'exit-code.txt').exists():
                if time.monotonic() >= deadline:
                    raise TimeoutError('Prior Terra queue exceeded 74 hours')
                time.sleep(30)
            record['prior_bem_exit_code'] = int((args.prior_bem/'exit-code.txt').read_text())
            pinned(args.verifier, VERIFIER)
            if digest(snapshot) != args.input_sha256 or digest(mesh) != record['mesh_sha256']:
                raise ValueError('BEM input changed during queue wait')
            command = [sys.executable, '-B', str(args.verifier/'scripts/bem_convergence_scan.py'),
                str(mesh), '--panels', 'full', '--epsilon-factors', '0.6', '--viscosity', '1',
                '--center', '0', '0', '0', '--backend', 'quotient-cg', '--subdivision-levels', '0',
                '--matrix-free-device', 'cpu', '--solver-rtol', '1e-6', '--solver-maxiter', '6000',
                '--matrix-free-row-chunk', '32', '--preconditioner-panel-block-size', '256',
                '--solver-progress', '--json', str(args.output/'scan.json'),
                '--csv', str(args.output/'scan.csv')]
            environment = dict(os.environ, KELVIN_SYMMETRY='cubic', OMP_NUM_THREADS='1',
                OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
            solver = args.output/'solver'
            solver.mkdir()
            with (solver/'run.log').open('w') as log:
                proc = subprocess.Popen(command, cwd=solver, env=environment,
                    stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                record.update(state='independent_bem', solver_pid=proc.pid, command=command)
                checkpoint()
                try:
                    actual = proc.wait(timeout=72*3600)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGTERM)
                    try:
                        actual = proc.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid, signal.SIGKILL)
                        actual = proc.wait()
                    record['timed_out'] = True
            record['solver_exit_code'] = actual
            (solver/'exit-code.txt').write_text(str(actual)+'\n')
            if actual:
                raise RuntimeError('Independent BEM failed; retain actual exit and log')
            record['independent_bem'] = accepted_scan(args.output, mesh)
            if record['independent_bem']['original_panels'] != record['export_certificate']['gamma_panels']:
                raise ValueError('BEM did not use all exported original panels')
            if digest(snapshot) != args.input_sha256:
                raise ValueError('Original HDF5 changed')
            pinned(controller, args.commit)
            record['state'] = 'completed'
            code = 0
    except Exception as error:
        record.update(state='failed', error=repr(error))
    finally:
        record.update(exit_code=code, finished_at=time.time())
        checkpoint()
        (args.output/'exit-code.txt').write_text(str(code)+'\n')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
