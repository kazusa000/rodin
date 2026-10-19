#!/usr/bin/env python3
"""Run the archived Terra numerical pipeline with an explicitly authorized Cubic seed."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

BASE = 'e004696706c169a775f039dc0782b07205c370da'
VERIFIER = '12f38a9a33dbbb0a1b7c0d0affcfff91620b927a'
EXPORTER = '2b3466878cd8c063a5c09ef5c5cd17a4bcf1ba75'
EXAMPLE = 'examples/ShapeOptimization/KelvinBall/'
NUMERICAL_ARGUMENTS = ['--n=20', '--outer-radius=2', '--regularization=100',
    '--step=0.5', '--advection-quadrature=3', '--reconstruction=mmg',
    '--mmg-adapt', '--mmg-snap=0.005', '--thickness-min=0',
    '--iterations=20', '--save-mesh']


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def git(source, *args):
    return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()


def pinned(source, commit):
    if git(source, 'rev-parse', 'HEAD') != commit or git(source, 'status', '--porcelain', '--untracked-files=all'):
        raise ValueError('Require the recorded clean commit: ' + str(source))


def save(path, value):
    path = Path(path)
    pending = path.with_name(path.name + '.pending')
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    pending.replace(path)


def complete_rows(path):
    if not path.exists():
        return []
    with path.open() as stream:
        rows = [r for r in csv.DictReader(stream)
                if None not in r and r and all(v not in (None, '') for v in r.values())]
    for i, row in enumerate(rows):
        if int(float(row['iteration'])) != i:
            raise ValueError('Nonconsecutive native CSV')
        for key, expected in [('outer_radius', 2), ('h', 2/19), ('dt', .5*2/19),
                              ('regularization_length', 100*2/19), ('advection_quadrature_order', 3)]:
            if not math.isclose(float(row[key]), expected, rel_tol=1e-12, abs_tol=1e-14):
                raise ValueError('Native parameter changed: ' + key)
        if not all(math.isfinite(float(row[k])) for k in ('k', 'c', 'q', 'rho')):
            raise ValueError('Nonfinite native diagnostics')
    return rows


def accepted_scan(directory, mesh):
    if (directory/'solver/exit-code.txt').read_text().strip() != '0':
        raise ValueError('BEM solver did not exit zero')
    rows = json.loads((directory/'scan.json').read_text())
    meta = json.loads((directory/'scan.metadata.json').read_text())
    if len(rows) != 1:
        raise ValueError('Expected exactly one full-panel BEM measurement')
    row = rows[0]
    for key, expected in [('status', 'ok'), ('symmetry', 'cubic'), ('subdivision_levels', 0),
                          ('requested_panels', 0), ('solver_converged', True), ('solver_rtol', 1e-6)]:
        if row.get(key) != expected:
            raise ValueError('BEM protocol mismatch: ' + key)
    if (row['panels'] != row['original_panels'] or row['full_panels'] != 24*row['original_panels']
            or row['epsilon_factor'] != .6 or row['backend'] != 'quotient-cg'):
        raise ValueError('BEM panel/operator protocol mismatch')
    rhs = row.get('rhs_relative_residuals')
    if not isinstance(rhs, list) or len(rhs) != 2:
        raise ValueError('Missing BEM RHS residuals')
    if not all(isinstance(x, (int, float)) and math.isfinite(x) and 0 <= x <= 1e-6
               for x in [row.get('linear_relative_residual'), *rhs]):
        raise ValueError('BEM residual exceeds actual 1e-6 gate')
    if meta.get('mesh_sha256') != digest(mesh) or meta.get('source_commit') != VERIFIER or meta.get('source_status'):
        raise ValueError('BEM input/source provenance mismatch')
    for key, expected in [('solver_rtol', 1e-6), ('viscosity', 1), ('center', [0, 0, 0]),
                          ('subdivision_levels', '0'), ('panels', 'full'), ('matrix_free_device', 'cpu')]:
        if meta['arguments'].get(key) != expected:
            raise ValueError('BEM argument provenance mismatch: ' + key)
    if (not all(math.isfinite(float(row[k])) for k in ('k', 'c', 'q', 'rho'))
            or row['k'] <= 0 or row['q'] <= 0
            or not math.isclose(row['rho'], abs(row['c'])/math.sqrt(row['k']*row['q']), rel_tol=1e-9, abs_tol=1e-14)):
        raise ValueError('BEM resistance values are inconsistent')
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'prefix', 'build', 'output', 'verifier', 'exporter', 'prior_bem']:
        parser.add_argument('--' + key.replace('_', '-'), type=Path, required=True)
    parser.add_argument('--commit', required=True)
    args = parser.parse_args()
    if platform.node() != '28h':
        raise ValueError('All computation is authorized only on 28h')
    pinned(args.source, args.commit)
    pinned(args.verifier, VERIFIER)
    pinned(args.exporter, EXPORTER)
    allowed = {EXAMPLE+p for p in ['Common.cpp', 'KelvinBallOptimization.cpp', 'SymmetryTest.cpp']}
    changed = git(args.source, 'diff', '--name-only', BASE, args.commit).splitlines()
    if any(p not in allowed and not p.startswith(EXAMPLE+'reproduction/') for p in changed):
        raise ValueError('Unapproved changes outside seed, its diagnostic/test, and reproduction wrapper')
    if git(args.source, 'diff', BASE, args.commit, '--', 'src'):
        raise ValueError('The numerical framework was changed')
    old_log = git(args.source, 'show', BASE+':'+EXAMPLE+'KelvinBallOptimization.cpp')
    actual_log = (args.source/EXAMPLE/'KelvinBallOptimization.cpp').read_text().strip()
    if old_log.replace(': "none")', ': "0.02*H3*H6")') != actual_log:
        raise ValueError('The old update algorithm changed beyond the seed diagnostic string')
    submodules = git(args.source, 'submodule', 'status', 'third-party/mmg', 'third-party/termcolor')
    if any(line.startswith(('-', '+', 'U')) for line in submodules.splitlines()):
        raise ValueError('Original MMG/termcolor dependencies must be initialized and pinned')
    for path in [args.output, args.build]:
        if path.exists() or path.resolve().is_relative_to(args.source.resolve()):
            raise ValueError('Use fresh build/output outside source')
    args.output.mkdir(parents=True)
    record = dict(source=str(args.source), commit=args.commit, original_terra_core=BASE,
        verifier_commit=VERIFIER, exporter_commit=EXPORTER, changed_files=changed,
        original_update_core_verified=True, submodules=submodules, symmetry='cubic',
        initial_shape='fixed 0.02*H3*H6 chiral perturbation; no seed selection',
        requested_states=20, requested_real_updates=19, parameters=NUMERICAL_ARGUMENTS,
        bem_feedback=False, started_at=time.time(), state='building', stages=[])
    environment = dict(os.environ, KELVIN_SYMMETRY='cubic', OMP_NUM_THREADS='6',
        OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    lib = args.prefix/'lib/x86_64-linux-gnu'
    environment['LD_LIBRARY_PATH'] = ':'.join(str(p) for p in [lib, lib/'openmpi/lib', lib/'hdf5/serial'])
    record['environment'] = {k:environment[k] for k in ['KELVIN_SYMMETRY', 'OMP_NUM_THREADS',
        'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'LD_LIBRARY_PATH']}
    def checkpoint():
        save(args.output/'campaign.json', record)
    def run(name, command, hours, env=None):
        directory = args.output/name
        directory.mkdir(parents=True)
        item = dict(name=name, command=[str(x) for x in command], started_at=time.time())
        record['stages'].append(item)
        with (directory/'run.log').open('w') as log:
            proc = subprocess.Popen(command, cwd=directory, env=env or environment,
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            item['pid'] = proc.pid
            checkpoint()
            try:
                code = proc.wait(timeout=hours*3600)
            except subprocess.TimeoutExpired:
                item['timed_out'] = True
                os.killpg(proc.pid, signal.SIGTERM)
                try: code = proc.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    code = proc.wait()
        item.update(exit_code=code, finished_at=time.time())
        (directory/'exit-code.txt').write_text(str(code)+'\n')
        checkpoint()
        return code
    code = 1
    try:
        configure = ['cmake', '-S', str(args.source/EXAMPLE/'reproduction'), '-B', str(args.build),
            '-DCMAKE_BUILD_TYPE=Release', f'-DCMAKE_PREFIX_PATH={args.prefix}', f'-DMUMPS_ROOT={args.prefix}',
            f'-DMPI_C_HEADER_DIR={lib/"openmpi/include"}',
            f'-DMPI_C_COMPILER_INCLUDE_DIRS={lib/"openmpi/include"};{lib/"openmpi/include/openmpi"}',
            '-DMPI_C_LIB_NAMES=mpi', f'-DMPI_mpi_LIBRARY={lib/"openmpi/lib/libmpi.so"}',
            f'-DCARLOS_MUMPS_SUPPORT_LIBRARIES={lib/"openmpi/lib/libmpi_cxx.so"}']
        if run('configure', configure, 1):
            raise RuntimeError('Configure failed; retain actual log/exit')
        if run('compile', ['cmake', '--build', str(args.build), '--target', 'KelvinBall',
                          'KelvinBallCubicGeometryTest', 'KelvinBallTetraGeometryTest', '-j', '2'], 4):
            raise RuntimeError('Compile failed; retain actual log/exit')
        for symmetry in ['Cubic', 'Tetra']:
            if run('validation/'+symmetry.lower(), [str(args.build/f'bin/KelvinBall{symmetry}GeometryTest')], .25):
                raise RuntimeError(symmetry+' geometry/seed control failed')
        record['symmetry_seed_tests_passed'] = True
        binary = args.build/'bin/KelvinBall'
        record['binary_sha256'] = digest(binary)
        record['state'] = 'optimizing'
        native = run('optimization', [str(binary), *NUMERICAL_ARGUMENTS], 12)
        record['native_exit_code'] = native
        rows = complete_rows(args.output/'optimization/kelvin-ball.csv')
        record.update(fully_evaluated_states=len(rows), completed_real_updates=max(0, len(rows)-1),
                      latest_fem_diagnostic=rows[-1] if rows else None)
        record['native_budget_completed'] = native == 0 and len(rows) == 20
        snapshots = sorted((args.output/'optimization').glob('KelvinBall.Chamber.mesh.*.h5'))
        if not snapshots:
            raise RuntimeError('No original saved geometry for independent evaluation')
        final = snapshots[-1]
        step = int(final.name.split('.')[-2])
        record.update(selected_bem_step=step, selected_bem_geometry_fem_evaluated=step < len(rows))
        record['state'] = 'exporting_saved_geometry'
        checkpoint()
        sys.path.insert(0, str(args.exporter/'scripts'))
        os.environ['KELVIN_SYMMETRY'] = 'cubic'
        os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
        spec = importlib.util.spec_from_file_location('original_export', args.exporter/'scripts/reproduce_smallbox_original.py')
        exporter = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(exporter)
        directory = args.output/'bem'
        directory.mkdir()
        native_ascii = directory/'native-export.mesh'
        exporter.export_snapshot(final, native_ascii)
        mesh = directory/'original-double.mesh'
        record['export_certificate'] = exporter.original_double_mesh(final, native_ascii, mesh)
        record['state'] = 'waiting_for_prior_terra_bem'
        checkpoint()
        deadline = time.monotonic()+74*3600
        while not (args.prior_bem/'exit-code.txt').exists():
            if time.monotonic() >= deadline:
                raise TimeoutError('Prior Terra BEM exit missing after bounded queue wait')
            time.sleep(30)
        record['prior_bem_exit_code'] = int((args.prior_bem/'exit-code.txt').read_text().strip())
        pinned(args.verifier, VERIFIER)
        bem_environment = dict(environment, OMP_NUM_THREADS='1')
        record['state'] = 'independent_bem'
        command = [sys.executable, '-B', str(args.verifier/'scripts/bem_convergence_scan.py'), str(mesh),
            '--panels', 'full', '--epsilon-factors', '0.6', '--viscosity', '1', '--center', '0', '0', '0',
            '--backend', 'quotient-cg', '--subdivision-levels', '0', '--matrix-free-device', 'cpu',
            '--solver-rtol', '1e-6', '--solver-maxiter', '6000', '--matrix-free-row-chunk', '32',
            '--preconditioner-panel-block-size', '256', '--solver-progress', '--json', str(directory/'scan.json'),
            '--csv', str(directory/'scan.csv')]
        if run('bem/solver', command, 72, bem_environment):
            raise RuntimeError('Independent final BEM failed; native outcome retained separately')
        record['independent_bem'] = accepted_scan(directory, mesh)
        pinned(args.source, args.commit)
        record['state'] = 'completed' if record['native_budget_completed'] else 'native_failed_final_bem_completed'
        code = 0 if record['native_budget_completed'] else 1
    except Exception as error:
        record.update(state='failed', error=repr(error))
    finally:
        try:
            rows = complete_rows(args.output/'optimization/kelvin-ball.csv')
            record.update(fully_evaluated_states=len(rows), completed_real_updates=max(0, len(rows)-1),
                          latest_fem_diagnostic=rows[-1] if rows else None)
        except Exception as error:
            record['csv_error'] = repr(error)
        record.update(exit_code=code, finished_at=time.time())
        checkpoint()
        (args.output/'exit-code.txt').write_text(str(code)+'\n')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
