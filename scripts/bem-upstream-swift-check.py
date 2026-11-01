#!/usr/bin/env python3
"""Bounded validation of native chamber BEM in upstream KelvinBall SWIFT.

Python only launches executables and compares artifacts; all numerical states,
derivatives, H1 projection, transport and reconstruction are Rodin C++.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def same_checkpoint(first, second):
    import h5py
    import numpy as np
    count = 0
    for name in ('chamber.h5', 'background.h5', 'state.h5'):
        with h5py.File(first / name) as a, h5py.File(second / name) as b:
            left, right = {}, {}
            def collect(root, result):
                def visit(key, obj):
                    if isinstance(obj, h5py.Dataset):
                        result[key] = obj[()]
                root.visititems(visit)
            collect(a, left)
            collect(b, right)
            if left.keys() != right.keys():
                raise ValueError(f'{name}: checkpoint dataset names changed')
            for key in left:
                if not np.array_equal(left[key], right[key]):
                    raise ValueError(f'{name}:{key}: checkpoint data changed')
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin', type=Path, required=True)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-rho', type=float, required=True)
    parser.add_argument('--mesh-sha256', required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(exist_ok=False)
    record = {'status': 'running', 'accepted_updates': 0,
              'design_domain': 'O24-chamber', 'gradient_domain': 'O24-chamber',
              'reconstruction_domain': 'O24-chamber', 'bem_domain': 'infinite-exterior',
              'step': 2, 'H1': 4, 'volume_projection_enabled': True,
              'volume_rejection_enabled': False, 'volume_tolerance': None,
              'cases': {}}
    save(root / 'validation.json', record)
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
           'MKL_NUM_THREADS': '1'}
    executable = args.bin.resolve() / 'KelvinBall'
    common = ['--bem-gradient', '--n=13', '--seed-amplitude=0.1', '--step=2',
              '--regularization=4', '--bem-order=3', '--bem-epsilon-factor=0.6',
              '--reconstruction=swift', '--mmg-adapt', '--save-checkpoints']

    def run(name, argv):
        directory = root / name
        directory.mkdir()
        command = [str(value) for value in argv]
        save(directory / 'command.json', command)
        start = time.monotonic()
        with (directory / 'run.log').open('w') as log:
            process = subprocess.Popen(command, cwd=directory, env=env,
                                       stdout=log, stderr=subprocess.STDOUT)
            (directory / 'native.pid').write_text(str(process.pid) + '\n')
            code = process.wait()
        (directory / 'run.exit').write_text(str(code) + '\n')
        record['cases'][name] = {'returncode': code, 'seconds': time.monotonic() - start}
        save(root / 'validation.json', record)
        if code:
            raise RuntimeError(f'{name}: native return code {code}; see {directory}/run.log')
        return directory

    try:
        actual_hash = hashlib.sha256(args.mesh.read_bytes()).hexdigest()
        if actual_hash != args.mesh_sha256:
            raise ValueError('Input checkpoint hash changed')
        record['input_sha256'] = actual_hash
        reference = run('reference', [args.bin.resolve() / 'KelvinBallBEMChamberValidation', '--reference'])
        record['reference'] = json.loads((reference / 'bem-chamber-validation.json').read_text())
        if not record['reference']['passed']:
            raise ValueError('Native reference / finite-difference test failed')
        initial = run('migration', [executable, *common, f'--mesh={args.mesh.resolve()}', '--geometry-only'])
        first = initial / 'checkpoint-000000'
        reloaded = run('roundtrip', [executable, *common, f'--restart-checkpoint={first}', '--geometry-only'])
        second = reloaded / 'checkpoint-000000'
        record['roundtrip_equal_datasets'] = same_checkpoint(first, second)
        fixed = run('fixed-state', [executable, *common, f'--restart-checkpoint={second}', '--state-only'])
        state = json.loads((fixed / 'bem-state-0.json').read_text())
        record['fixed_state'] = state
        if not math.isclose(state['rho'], args.expected_rho, rel_tol=1e-8, abs_tol=1e-18):
            raise ValueError('Last accepted rho changed after framework migration')
        update = run('update', [executable, *common, f'--restart-checkpoint={second}', '--iterations=2'])
        states = [json.loads((update / f'bem-state-{i}.json').read_text()) for i in (0, 1)]
        with (update / 'kelvin-ball.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        if len(rows) != 2:
            raise ValueError('Expected two native state evaluations')
        old, new = [s['rho'] for s in states]
        predicted = float(rows[1]['predicted_delta_rho'])
        slope = float(rows[0]['d_volume_theta'])
        fit = json.loads((update / 'swift-fit-latest.json').read_text())
        checks = {
            'restart_rho_matches': math.isclose(old, args.expected_rho, rel_tol=1e-8, abs_tol=1e-18),
            'strict_rho_increase': new > old,
            'weak_armijo': predicted > 0 and new >= old + 1e-4 * predicted,
            'closed_surface': all(s['open_edges'] == 0 for s in states),
            'state_residual': all(math.isfinite(s['residual']) and s['residual'] <= 1e-9 for s in states),
            'volume_projection': math.isfinite(slope) and abs(slope) <= 1e-10,
            'real_transport': float(rows[0]['advection_increment']) > 1e-12,
            'positive_volume': all(float(r['volume']) > 0 for r in rows),
            'swift_quality_budget': fit['quality_budget_satisfied'],
            'complete_candidate_checkpoint': (update / 'checkpoint-000001/READY').is_file(),
        }
        record.update(checks=checks, initial_rho=old, updated_rho=new,
                      predicted_delta_rho=predicted, volume_slope=slope, swift_fit=fit,
                      relative_volume_change=float(rows[1]['volume']) / float(rows[0]['volume']) - 1,
                      geometric_target_satisfied=fit['geometric_sup'] <= fit['geometric_target'])
        if not all(checks.values()):
            raise ValueError('Candidate failed: ' + ', '.join(k for k, v in checks.items() if not v))
        record.update(status='passed', accepted_updates=1,
                      last_accepted_checkpoint=str(update / 'checkpoint-000001'))
        save(root / 'validation.json', record)
        return 0
    except Exception as error:
        record.update(status='failed', error=str(error))
        save(root / 'validation.json', record)
        print(error, flush=True)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
