#!/usr/bin/env python3
# Copyright Carlos BRITO PACHECO 2021 - 2026.
# Distributed under the Boost Software License, Version 1.0.
"""Subprocess regression checks for bounded acceptance; no physical solver fixtures."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

DRIVER = Path(__file__).resolve().parents[1] / 'bem-gradient-chamber-run.py'
FAKE = '''#!/usr/bin/env python3
import csv,json,pathlib,sys
p=pathlib.Path('.')
config=json.loads(pathlib.Path(__file__).with_name('policy.json').read_text())
if config.get('native_exit'):raise SystemExit(config['native_exit'])
args=dict(x[2:].split('=',1) for x in sys.argv[1:] if '=' in x)
old=float(pathlib.Path(args['mesh']).read_text());factor=float(args['step'])
new=old*1.1 if factor<=config['accept_at'] else old*.9
s=dict(rho=old,open_edges=0,residual=1e-15,design_domain='O24-chamber',state_dofs=100,independent_dofs=80)
(p/'bem-state-0.json').write_text(json.dumps(s));s['rho']=new;(p/'bem-state-1.json').write_text(json.dumps(s))
(p/'bem-chamber-0.mesh').write_text(str(old));(p/'bem-chamber-1.mesh').write_text(str(new))
with (p/'kelvin-ball.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['predicted_delta_rho','d_volume_theta','advection_increment','volume']);w.writerow(['nan',config.get('volume_slope',0),.01,1]);w.writerow([old*.1,0,0,2])
'''


class ChamberRun(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fake = self.root / 'native'
        self.fake.write_text(FAKE)
        self.fake.chmod(0o700)
        self.seed = self.root / 'seed.mesh'
        self.seed.write_text('1e-6')
        self.results = self.root / 'results'

    def run_driver(self, policy, *options):
        (self.root / 'policy.json').write_text(json.dumps(policy))
        cmd = [sys.executable, str(DRIVER), '--executable', str(self.fake),
               '--results', str(self.results), '--checkpoint', str(self.seed),
               '--expected-rho', '1e-6', '--updates', '1', '--step', '1.0',
               '--regularization', '4', '--max-backtracks', '2', *options]
        result = subprocess.run(cmd, capture_output=True, text=True)
        history = self.results / 'history.json'
        return result.returncode, json.loads(history.read_text()) if history.exists() else None

    def test_second_halving_accepts_from_same_checkpoint(self):
        code, h = self.run_driver({'accept_at': .25})
        self.assertEqual(code, 0)
        self.assertTrue(h['complete'])
        self.assertEqual(h['last_step_factor'], .25)
        commands = [json.loads((self.results / f'iteration-0/trial-{i}/update/command.json').read_text()) for i in range(3)]
        self.assertEqual([next(x for x in c if x.startswith('--step=')) for c in commands], ['--step=1.0', '--step=0.5', '--step=0.25'])
        self.assertEqual(len({next(x for x in c if x.startswith('--mesh=')) for c in commands}), 1)
        self.assertEqual((self.results / 'accepted-0.mesh').read_text(), '1e-6')
        self.assertAlmostEqual(float((self.results / 'accepted-1.mesh').read_text()), 1.1e-6)
        self.assertEqual(len(h['rejected_trials']), 2)
        # Large volume drift is observed, never used as a rejection cap.
        record = json.loads((self.results / 'iteration-0/trial-2/chamber-update-validation.json').read_text())
        self.assertEqual(record['relative_chamber_volume_change'], 1)
        self.assertTrue(record['passed'])

    def test_exhaustion_preserves_last_accepted(self):
        code, h = self.run_driver({'accept_at': 0})
        self.assertEqual(code, 2)
        self.assertEqual(h['status'], 'line_search_exhausted')
        self.assertEqual(h['accepted_updates'], 0)
        self.assertFalse((self.results / 'accepted-1.mesh').exists())
        self.assertEqual((self.results / 'accepted-0.mesh').read_text(), '1e-6')
        self.assertEqual(len(list(self.results.glob('iteration-0/trial-*'))), 3)

    def test_skip_known_failure_only_on_first_iteration(self):
        code, h = self.run_driver({'accept_at': .5}, '--first-backtrack', '1', '--updates', '2')
        self.assertEqual(code, 0)
        self.assertEqual(h['accepted_updates'], 2)
        self.assertFalse((self.results / 'iteration-0/trial-0').exists())
        self.assertTrue((self.results / 'iteration-1/trial-0/update/run.log').exists())
        self.assertTrue((self.results / 'iteration-1/trial-1/update/run.log').exists())
        cmd = json.loads((self.results / 'iteration-1/trial-0/update/command.json').read_text())
        self.assertIn('--mesh=' + str(self.results / 'accepted-1.mesh'), cmd)

    def test_native_failures_are_bounded_and_preserved(self):
        code, h = self.run_driver({'accept_at': .5, 'native_exit': 7})
        self.assertEqual(code, 7)
        self.assertEqual(h['status'], 'line_search_native_failure')
        self.assertEqual(len(h['rejected_trials']), 3)
        self.assertTrue(all(t['native_returncode'] == 7 for t in h['rejected_trials']))
        self.assertFalse((self.results / 'accepted-1.mesh').exists())

    def test_default_native_failure_still_stops_immediately(self):
        code, h = self.run_driver({'accept_at': .5, 'native_exit': 7}, '--max-backtracks', '0')
        self.assertEqual(code, 7)
        self.assertEqual(h['status'], 'native_failure')
        self.assertEqual(len(h['rejected_trials']), 1)

    def test_lost_projection_is_not_retried(self):
        code, h = self.run_driver({'accept_at': 1, 'volume_slope': .01})
        self.assertEqual(code, 2)
        self.assertEqual(h['status'], 'candidate_rejected')
        self.assertFalse((self.results / 'iteration-0/trial-1').exists())

    def test_invalid_backtrack_count_creates_no_run(self):
        code, h = self.run_driver({'accept_at': .5}, '--first-backtrack', '3')
        self.assertEqual(code, 2)
        self.assertIsNone(h)
        self.assertFalse(self.results.exists())


if __name__ == '__main__':
    unittest.main()
