"""Protect the complete equal-dataset Cell summary used by author charts."""
import csv
import tempfile
import unittest
from pathlib import Path

from evaluation.summary.refresh_cell_summary import build_summary, quantile, refresh


class CellSummaryTests(unittest.TestCase):
    def test_linear_quartiles(self):
        self.assertEqual(quantile([1, 2, 3, 4], .25), 1.75)
        self.assertEqual(quantile([1, 2, 3, 4], .75), 3.25)

    def test_committed_statistics_match_canonical_cell(self):
        rows = refresh()
        self.assertEqual(len(rows), 100)
        self.assertEqual({r["n_datasets"] for r in rows}, {14})

    def test_missing_dataset_is_rejected(self):
        root = Path(__file__).resolve().parents[1]
        with (root / "evaluation/cell/results/cell_metrics.csv").open(newline="") as handle:
            reader = csv.DictReader(handle)
            fields, rows = reader.fieldnames, list(reader)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "incomplete.csv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows[:-1])
            with self.assertRaises(ValueError):
                build_summary(path)


if __name__ == "__main__":
    unittest.main()
