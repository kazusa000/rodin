#!/usr/bin/env python3
"""Compare native chamber reconstructions of one frozen transported P1 target.

This orchestrates native executables only. It never promotes a candidate or
starts an optimization continuation. The pre-transport control still includes
Eikonal redistancing and transfer; it is not an identity mesh transform.
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--expected-rho', type=float, required=True)
    parser.add_argument('--timeout', type=float, default=1800)
    args = parser.parse_args()
    if not all(p.is_absolute() for p in (args.executable, args.checkpoint, args.results)):
        parser.error('Use absolute paths')
    args.results.mkdir(exist_ok=False)
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}
    common = [str(args.executable), '--bem-gradient', '--n=13', '--seed-amplitude=0.1',
              '--step=2', '--regularization=4', '--mmg-retries=0']

    def run(name, flags):
        directory = args.results / name
        directory.mkdir()
        command = common + flags
        save(directory / 'command.json', command)
        start = time.monotonic()
        with (directory / 'run.log').open('w') as log:
            process = subprocess.Popen(command, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT)
            (directory / 'pid').write_text(str(process.pid) + '\n')
            timed_out = False
            try:
                code = process.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                process.terminate()
                try:
                    code = process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    code = process.wait()
        record = {'returncode': code, 'timed_out': timed_out, 'seconds': time.monotonic() - start}
        save(directory / 'result.json', record)
        (directory / 'run.exit').write_text(str(code) + '\n')
        return directory, record

    save(args.results / 'inputs.json', {
        'checkpoint': str(args.checkpoint),
        'checkpoint_sha256': hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),
        'binary_sha256': hashlib.sha256(args.executable.read_bytes()).hexdigest(),
        'step': 2, 'H1': 4, 'volume_projection_enabled': True,
        'volume_rejection_enabled': False, 'volume_tolerance': None,
        'threads': 1, 'accepted_updates': 0,
    })
    export, result = run('transport', ['--iterations=2', '--bem-export-transport',
        '--reconstruction=mmg', '--mmg-snap=0.001', f'--mesh={args.checkpoint}'])
    if result['returncode']:
        return 1
    old = json.loads((export / 'bem-state-0.json').read_text())
    with (export / 'kelvin-ball.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError('Expected one export state')
    row = rows[0]
    slope = float(row['d_volume_theta'])
    predicted = float(row['dt']) * float(row['d_rho_theta'])
    if not (math.isclose(old['rho'], args.expected_rho, rel_tol=1e-9) and
            abs(slope) <= 1e-10 and predicted > 0 and old['residual'] <= 1e-9):
        raise ValueError('Export violated the restart/state/volume-projection contract')
    results = {'initial': old, 'predicted_delta_rho': predicted,
               'volume_slope': slope, 'accepted_updates': 0, 'cases': []}
    for target in ('target.sol', 'transport-initial.sol'):
        for backend in ('mmg', 'swift'):
            name = ('advected' if target == 'target.sol' else 'pretransport') + '-' + backend
            flags = ['--reconstruct-only', f'--reconstruction={backend}',
                     f'--mesh={export / "transport.mesh"}', f'--target={export / target}']
            if backend == 'mmg':
                flags += ['--mmg-snap=0.001']
            directory, reconstruction = run(name, flags)
            case = {'name': name, 'reconstruction': reconstruction,
                    'target_sha256': hashlib.sha256((export / target).read_bytes()).hexdigest(),
                    'background_sha256': hashlib.sha256((export / 'transport.mesh').read_bytes()).hexdigest()}
            if (directory / 'swift-fit.json').exists():
                case['swift_fit'] = json.loads((directory / 'swift-fit.json').read_text())
            if reconstruction['returncode'] == 0:
                scored, scoring = run(name + '-state', ['--iterations=1', '--state-only',
                    '--reconstruction=mmg', f'--mesh={directory / "candidate.mesh"}'])
                case['scoring'] = scoring
                if scoring['returncode'] == 0:
                    state = json.loads((scored / 'bem-state-0.json').read_text())
                    with (scored / 'kelvin-ball.csv').open() as stream:
                        newrow = next(csv.DictReader(stream))
                    case['state'] = state
                    case['relative_rho_change'] = state['rho'] / old['rho'] - 1
                    case['relative_volume_change'] = float(newrow['volume']) / float(row['volume']) - 1
                    case['checks'] = {
                        'closed_interface': state['open_edges'] == 0,
                        'state_residual': state['residual'] <= 1e-9,
                        'strict_rho_increase': state['rho'] > old['rho'],
                        'weak_armijo': state['rho'] >= old['rho'] + 1e-4 * predicted,
                        'positive_volume': float(newrow['volume']) > 0,
                        'swift_fit_certified': backend != 'swift' or all(case['swift_fit'][k]
                            for k in ('quality_budget_satisfied', 'geometric_target_reached')),
                    }
                    case['eligible_advected_candidate'] = target == 'target.sol' and all(case['checks'].values())
            results['cases'].append(case)
            save(args.results / 'comparison.json', results)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
