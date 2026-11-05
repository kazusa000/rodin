# Distributed under the Boost Software License, Version 1.0.
import csv
import io
from pathlib import Path
import tempfile
import unittest

from run_l10_n50 import complete_history


class HistoryTest(unittest.TestCase):
    def row(self, iteration=0):
        h = 10/49
        return dict(iteration=iteration, outer_radius=10, initial_grid_spacing=h,
                    hmin_factor=.05, hmax_factor=10, dt=.2*h,
                    regularization_length=20*h, thickness_min=.05,
                    k=20, c=.01, q=50, rho=.0003, volume=.17, target_volume=.17)

    def read(self, rows, tail=""):
        output = io.StringIO()
        writer = csv.DictWriter(output, fieldnames=self.row().keys())
        writer.writeheader()
        writer.writerows(rows)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"history.csv"
            path.write_text(output.getvalue()+tail)
            return complete_history(path)

    def test_partial_unsolved_line_excluded(self):
        self.assertEqual(len(self.read([self.row()], "1,10,0.204")), 1)

    def test_iteration_gap_rejected(self):
        with self.assertRaises(ValueError):
            self.read([self.row(), self.row(2)])

    def test_previous_run_parameters_rejected(self):
        row = self.row()
        row["initial_grid_spacing"] = 10/59
        with self.assertRaises(ValueError):
            self.read([row])

    def test_nonfinite_score_rejected(self):
        row = self.row()
        row["rho"] = "nan"
        with self.assertRaises(ValueError):
            self.read([row])

    def test_incomplete_flushed_row_rejected(self):
        row = self.row()
        row["q"] = ""
        with self.assertRaises(ValueError):
            self.read([row])


if __name__ == "__main__":
    unittest.main()
