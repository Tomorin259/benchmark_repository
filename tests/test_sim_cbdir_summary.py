"""Protect Figure 7 against historical simulation scores and pooled averages."""

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from evaluation.summary.refresh_real_cbdir_summary import check_summary


class SimulationSummaryTests(unittest.TestCase):
    def test_committed_summary_matches_reference30(self):
        subprocess.run([sys.executable, "evaluation/summary/refresh_sim_cbdir_summary.py"],
                       cwd=Path(__file__).resolve().parents[1], check=True)

    def test_historical_and_pooled_scores_are_rejected(self):
        expected = dict(Metric="Sim-CBDir", Method="TIvelo", Mean="0.01407165424752784")
        for value in ("0.01116133333333332", "-0.2045037569606416"):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "summary.csv"
                with path.open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(expected))
                    writer.writeheader()
                    writer.writerow(dict(expected, Mean=value))
                before = path.read_bytes()
                with self.assertRaises(SystemExit):
                    check_summary(path, list(expected), [expected], metric="Sim-CBDir")
                self.assertEqual(path.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
