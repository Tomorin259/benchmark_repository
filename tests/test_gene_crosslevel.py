"""Check the complete S9 data replay against the registered assembly inputs."""

import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class CrosslevelTests(unittest.TestCase):
    def test_s9_replay(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "replay"
            subprocess.run([sys.executable, "evaluation/gene/reproduce_crosslevel.py",
                            "--output-root", str(out)], cwd=ROOT, check=True)
            for name in ("dataset_scores", "frontiers", "leave_one_out", "plot_data", "upper_means"):
                filename = name + ".csv"
                formal = ROOT / "evaluation/gene/results/crosslevel" / filename
                self.assertEqual((out / filename).read_bytes(), formal.read_bytes(), filename)
            for name, value, expected_count in (("plot_data.csv", "value", 423),
                                                ("upper_means.csv", "mean", 54)):
                scores = []
                for path in (out / name, ROOT / "figures/assembly/suppfig09" / name):
                    with path.open(newline="") as handle:
                        rows = list(csv.DictReader(handle))
                    self.assertEqual(len(rows), expected_count)
                    table = {tuple(v for k, v in row.items() if k != value): float(row[value])
                             for row in rows}
                    self.assertEqual(len(table), expected_count)
                    scores.append(table)
                self.assertEqual(scores[0].keys(), scores[1].keys())
                for key, actual in scores[0].items():
                    self.assertAlmostEqual(actual, scores[1][key], delta=1e-12, msg=str(key))

    def test_output_must_be_new(self):
        with tempfile.TemporaryDirectory() as temp:
            command = [sys.executable, "evaluation/gene/reproduce_crosslevel.py", "--output-root"]
            existing = subprocess.run([*command, temp], cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(existing.returncode, 0)
            self.assertIn("must be a new directory", existing.stderr)
            self.assertEqual(list(Path(temp).iterdir()), [])

            out = Path(temp) / "shared"
            runs = [subprocess.Popen([*command, str(out)], cwd=ROOT, stderr=subprocess.PIPE, text=True)
                    for _ in range(2)]
            results = [(run.communicate()[1], run.returncode) for run in runs]
            self.assertEqual(sorted(code for _, code in results), [0, 2])
            self.assertIn("must be a new directory", next(msg for msg, code in results if code))
            self.assertEqual(list(Path(temp).iterdir()), [out])
            self.assertEqual(len(list(out.glob("*.csv"))), 5)
