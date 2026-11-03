# Distributed under the Boost Software License, Version 1.0.
"""Check final-state selection under partial output and failed native runs."""
import csv
import tempfile
from pathlib import Path
import unittest

from run_carlos_l10 import complete_history, final_snapshot


class FinalStateTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.path = self.directory / "kelvin-ball.csv"
        self.keys = ["iteration", "k", "c", "q", "rho", "volume", "target_volume",
                     "outer_radius", "initial_grid_spacing", "dt", "regularization_length", "thickness_min"]

    def save(self, iterations=(0, 1), **override):
        with self.path.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=self.keys)
            writer.writeheader()
            for step in iterations:
                row = dict(iteration=step, k=20, c=.001, q=30, rho=.00004, volume=.174,
                           target_volume=.174, outer_radius=10, initial_grid_spacing=10/59,
                           dt=.2*10/59, regularization_length=20*10/59, thickness_min=.05)
                row.update(override)
                writer.writerow(row)

    def test_selects_final_solved_not_later_unsolved_geometry(self):
        self.save()
        for step in (1, 2):
            (self.directory / f"KelvinBall.Chamber.mesh.{step:06d}.h5").touch()
        self.assertEqual(final_snapshot(self.directory, complete_history(self.path))[0], 1)

    def test_ignores_partial_final_csv_line(self):
        self.save()
        with self.path.open("a") as stream:
            stream.write("2,20,.001")
        self.assertEqual(len(complete_history(self.path)), 2)

    def test_rejects_malformed_complete_row(self):
        self.save()
        with self.path.open("a") as stream:
            stream.write("2,20,.001\n")
        with self.assertRaises(ValueError):
            complete_history(self.path)

    def test_rejects_nonfinite_metric(self):
        self.save(rho="nan")
        with self.assertRaises(ValueError):
            complete_history(self.path)

    def test_rejects_wrong_actual_spacing(self):
        self.save(initial_grid_spacing=1/6)
        with self.assertRaises(ValueError):
            complete_history(self.path)

    def test_rejects_duplicate_state(self):
        self.save(iterations=(0, 0))
        with self.assertRaises(ValueError):
            complete_history(self.path)

    def test_rejects_missing_final_snapshot(self):
        self.save()
        with self.assertRaises(ValueError):
            final_snapshot(self.directory, complete_history(self.path))

    def test_rejects_no_completed_state(self):
        self.assertEqual(complete_history(self.path), [])
        with self.assertRaises(ValueError):
            final_snapshot(self.directory, [])


if __name__ == "__main__":
    unittest.main()
