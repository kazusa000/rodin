#!/usr/bin/env python3
# Copyright Carlos BRITO PACHECO 2021 - 2026.
# Distributed under the Boost Software License, Version 1.0.
"""Bounded chamber continuation; numerical states/gradients/updates remain native."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess


def save(path, record):
    """Atomically publish a JSON record at path."""
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def validate_update(root, expected=None):
    """Check two native state evaluations and their single actual MMG update."""
    u = root / 'update'
    states = [json.loads((u / f'bem-state-{i}.json').read_text()) for i in [0, 1]]
    with (u / 'kelvin-ball.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 2:
        raise ValueError('Expected exactly two chamber evaluations')
    old, new = [s['rho'] for s in states]
    predicted = float(rows[1]['predicted_delta_rho'])
    slope = float(rows[0]['d_volume_theta'])
    increment = float(rows[0]['advection_increment'])
    volumes = [float(r['volume']) for r in rows]
    if not all(math.isfinite(x) for x in [old, new, predicted, slope, increment, *volumes]):
        raise ValueError('Nonfinite native update diagnostics')
    checks = {
        'strict_rho_increase': new > old,
        'weak_armijo': predicted > 0 and new >= old + 1e-4 * predicted,
        'closed_rotated_interface': all(s['open_edges'] == 0 for s in states),
        'native_state_residual': all(math.isfinite(s['residual']) and s['residual'] <= 1e-9 for s in states),
        'chamber_state_space': all(s['design_domain'] == 'O24-chamber' and 0 < s['independent_dofs'] <= s['state_dofs'] for s in states),
        'volume_projection': abs(slope) <= 1e-10,
        'real_transport': increment > 1e-12,
        'chamber_checkpoints': all((u / f'bem-chamber-{i}.mesh').is_file() for i in [0, 1]),
        'restart_rho_matches': expected is None or math.isclose(old, expected, rel_tol=1e-8, abs_tol=1e-18),
        'positive_volume': all(v > 0 for v in volumes),
    }
    record = {
        'passed': all(checks.values()), 'design_domain': 'O24-chamber',
        'gradient_domain': 'O24-chamber', 'reconstruction_domain': 'O24-chamber',
        'bem_domain': 'infinite-exterior', 'checks': checks,
        'initial_rho': old, 'updated_rho': new, 'predicted_delta_rho': predicted,
        'state_evaluations': 2, 'real_updates': 1,
        'volume_projection_enabled': True, 'volume_rejection_enabled': False,
        'volume_tolerance': None, 'relative_chamber_volume_change': volumes[1] / volumes[0] - 1,
    }
    save(root / 'chamber-update-validation.json', record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-update', type=Path, help='Check an existing RESULTS/update directory only')
    parser.add_argument('--executable', type=Path)
    parser.add_argument('--results', type=Path)
    parser.add_argument('--checkpoint', type=Path)
    parser.add_argument('--expected-rho', type=float)
    parser.add_argument('--updates', type=int, default=20, help='Maximum additional accepted updates')
    parser.add_argument('--step', type=float, default=0.10, help='Native step factor')
    parser.add_argument('--regularization', type=float, default=4, help='Native H1 length factor')
    parser.add_argument('--max-backtracks', type=int, default=0, help='Maximum step halvings per update')
    parser.add_argument('--first-backtrack', type=int, default=0, help='Skip previously rejected larger steps on the first update only')
    args = parser.parse_args()
    if args.check_update:
        record = validate_update(args.check_update)
        print(json.dumps(record, indent=2))
        return 0 if record['passed'] else 2
    if (not all([args.executable, args.results, args.checkpoint]) or
            args.expected_rho is None or not math.isfinite(args.expected_rho) or
            args.expected_rho <= 0 or args.updates <= 0):
        parser.error('Need executable, new results, checkpoint, positive expected rho and update count')
    if not all(math.isfinite(x) and x > 0 for x in [args.step, args.regularization]):
        parser.error('Step and H1 regularization factors must be finite and positive')
    if not 0 <= args.first_backtrack <= args.max_backtracks:
        parser.error('Need 0 <= first-backtrack <= max-backtracks')
    for p in [args.executable, args.results, args.checkpoint]:
        if not p.is_absolute():
            parser.error('Use absolute paths')
    if not args.executable.is_file() or not args.checkpoint.is_file():
        parser.error('Missing executable or checkpoint')
    args.results.mkdir(exist_ok=False)
    accepted = args.results / 'accepted-0.mesh'
    shutil.copy2(args.checkpoint, accepted)
    history = {
        'complete': False, 'status': 'running', 'requested_updates': args.updates,
        'accepted_updates': 0, 'starting_rho': args.expected_rho, 'last_rho': args.expected_rho,
        'last_accepted_checkpoint': str(accepted), 'design_domain': 'O24-chamber',
        'gradient_domain': 'O24-chamber', 'reconstruction_domain': 'O24-chamber',
        'bem_domain': 'infinite-exterior', 'initial_step_factor': args.step,
        'regularization_factor': args.regularization,
        'max_backtracks': args.max_backtracks, 'first_backtrack': args.first_backtrack,
        'step_reduction': 0.5, 'rejected_trials': [],
        'volume_projection_enabled': True, 'volume_rejection_enabled': False, 'volume_tolerance': None,
    }
    history_path = args.results / 'history.json'
    save(history_path, history)
    for i in range(args.updates):
        iteration = args.results / f'iteration-{i}'
        first_trial = args.first_backtrack if i == 0 else 0
        native_completed = False
        last_native_returncode = None
        for trial in range(first_trial, args.max_backtracks + 1):
            factor = args.step * 0.5 ** trial
            step = iteration / f'trial-{trial}' if args.max_backtracks else iteration
            u = step / 'update'
            u.mkdir(parents=True)
            history.update(current_iteration=i, current_trial=trial, current_step_factor=factor)
            save(history_path, history)
            command = [str(args.executable), '--bem-gradient', '--n=13', '--seed-amplitude=0.1',
                       '--iterations=2', f'--step={factor}', f'--regularization={args.regularization}',
                       '--reconstruction=mmg', '--mmg-snap=0.001', '--mmg-retries=0',
                       '--save-mesh', f'--mesh={accepted}']
            save(u / 'command.json', command)
            try:
                with (u / 'run.log').open('w') as log:
                    process = subprocess.run(command, cwd=u, stdout=log, stderr=subprocess.STDOUT,
                                             env={**os.environ, 'OMP_NUM_THREADS': '1'})
                save(u / 'native-result.json', {'returncode': process.returncode})
                if process.returncode != 0:
                    last_native_returncode = process.returncode
                    history['rejected_trials'].append({
                        'iteration': i, 'trial': trial, 'step_factor': factor,
                        'native_returncode': process.returncode,
                        'native_result': str(u / 'native-result.json')})
                    save(history_path, history)
                    print(f"rejected_update={i + 1} trial={trial} step={factor} native_returncode={process.returncode}", flush=True)
                    if args.max_backtracks == 0:
                        history.update(status='native_failure', stop_reason=f'native_returncode_{process.returncode}')
                        save(history_path, history)
                        return process.returncode if process.returncode > 0 else 128 - process.returncode
                    continue
                native_completed = True
                check = validate_update(step, history['last_rho'])
                if not check['passed']:
                    failed = [k for k, v in check['checks'].items() if not v]
                    history['rejected_trials'].append({
                        'iteration': i, 'trial': trial, 'step_factor': factor,
                        'reason': failed, 'validation': str(step / 'chamber-update-validation.json')})
                    save(history_path, history)
                    print(f"rejected_update={i + 1} trial={trial} step={factor} reason={','.join(failed)}", flush=True)
                    # Projection, restart and artifact contract errors cannot be
                    # repaired by changing the transport step, so remain fatal.
                    if not set(failed) <= {'strict_rho_increase', 'weak_armijo'}:
                        history.update(status='candidate_rejected', stop_reason=','.join(failed))
                        save(history_path, history)
                        return 2
                    continue
                accepted = args.results / f'accepted-{i + 1}.mesh'
                shutil.copy2(u / 'bem-chamber-1.mesh', accepted)
                history.update(accepted_updates=i + 1, last_rho=check['updated_rho'],
                               last_accepted_checkpoint=str(accepted), last_step_factor=factor,
                               last_accepted_trial=trial,
                               last_accepted_validation=str(step / 'chamber-update-validation.json'))
                save(history_path, history)
                print(f"accepted_update={i + 1} rho={history['last_rho']:.17g} step={factor} trial={trial}", flush=True)
                break
            except Exception as error:
                history.update(status='observation_error', stop_reason=str(error))
                save(history_path, history)
                raise
        else:
            if not native_completed and last_native_returncode is not None:
                history.update(status='line_search_native_failure',
                               stop_reason='native_failure_at_all_tested_steps')
                save(history_path, history)
                return last_native_returncode if last_native_returncode > 0 else 128 - last_native_returncode
            history.update(status='line_search_exhausted' if args.max_backtracks else 'candidate_rejected',
                           stop_reason='objective_rejected_at_all_tested_steps' if args.max_backtracks else
                           ','.join(history['rejected_trials'][-1]['reason']))
            save(history_path, history)
            return 2
    history.update(complete=True, status='requested_updates_complete', stop_reason='update_limit')
    save(history_path, history)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
