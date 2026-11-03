# Distributed under the Boost Software License, Version 1.0.
import json
import unittest
from unittest.mock import patch

from test_run_carlos_l10 import FinalStateTest
from evaluate_carlos_step30 import validate_input, FEM, BINARY_SHA, SNAPSHOT_SHA


class InputTest(FinalStateTest):
    def prepare(self):
        simulation = self.directory / "optimization"
        simulation.mkdir()
        self.path = simulation / "kelvin-ball.csv"
        self.save(iterations=range(31))
        (simulation / "KelvinBall.Chamber.mesh.000030.h5").write_bytes(b"checksum certified separately")
        campaign = self.directory / "campaign.json"
        campaign.write_text(json.dumps(dict(commit=FEM, binary_sha256=BINARY_SHA, status="optimization")))
        return campaign

    def test_accepts_completed_step_while_fem_continues(self):
        campaign = self.prepare()
        with patch("evaluate_carlos_step30.sha", return_value=SNAPSHOT_SHA):
            row, snapshot = validate_input(campaign)
        self.assertEqual(int(row["iteration"]), 30)
        self.assertTrue(snapshot.name.endswith("000030.h5"))

    def test_rejects_different_preview_geometry(self):
        campaign = self.prepare()
        with patch("evaluate_carlos_step30.sha", return_value="different"):
            with self.assertRaisesRegex(ValueError, "displayed shape"):
                validate_input(campaign)

    def test_rejects_another_fem_source(self):
        campaign = self.prepare()
        campaign.write_text(json.dumps(dict(commit="other", binary_sha256=BINARY_SHA)))
        with self.assertRaisesRegex(ValueError, "frozen"):
            validate_input(campaign)


if __name__ == "__main__":
    unittest.main()
