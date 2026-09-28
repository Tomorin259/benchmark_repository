"""Ensure integrated kinetics summaries follow the current full-gene evaluation."""

import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class KineticsSummaryTests(unittest.TestCase):
    def test_committed_summary_matches_current_evaluation(self):
        subprocess.run(
            [sys.executable, str(ROOT / "evaluation/summary/refresh_kinetics_summary.py")],
            check=True,
        )

    def test_current_organoid_reference_and_labelled_input_contract(self):
        config = json.loads((ROOT / "evaluation/kinetics_recovery/config.json").read_text())
        branches = config["datasets"]["organoid"]["branches"]
        self.assertTrue(branches["1"]["reference_reverse_time"])
        self.assertFalse(branches["2"]["reference_reverse_time"])
        dynamo = config["methods"]["Dynamo"]["organoid"]
        self.assertEqual(dynamo["n_obs"], 1373)
        self.assertTrue(dynamo["label_data"])
        self.assertEqual(dynamo["layers"], {"alpha": "cell_wise_alpha"})


if __name__ == "__main__":
    unittest.main()
