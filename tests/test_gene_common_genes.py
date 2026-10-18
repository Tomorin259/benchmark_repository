"""Protect finite nine-method intersection and equal-dataset S9 aggregation."""

import csv
import statistics
import subprocess
import sys
import unittest
from pathlib import Path

from evaluation.gene.refresh_common_gene_summary import common_gene_tables, EXPECTED_COMMON_COUNTS

ROOT = Path(__file__).resolve().parents[1]


class CommonGeneTests(unittest.TestCase):
    def test_committed_tables_match_finite_intersection(self):
        subprocess.run([sys.executable, "evaluation/gene/refresh_common_gene_summary.py"], cwd=ROOT, check=True)
        source = ROOT / "evaluation/gene/results/gene_cbdir_reference30_all.csv"
        tables = common_gene_tables(source)
        self.assertEqual(list(EXPECTED_COMMON_COUNTS.values()), [6, 6, 6, 7, 2, 4, 4])
        self.assertEqual(len(tables["membership.csv"]), 35)
        self.assertEqual(len(tables["dataset_summary.csv"]), 63)
        macro = {row["method"]: row for row in tables["macro_summary.csv"]}
        self.assertEqual(macro["veloVI"]["rank"], 5)
        dataset_means = [row["mean"] for row in tables["dataset_summary.csv"] if row["method"] == "veloVI"]
        self.assertEqual(macro["veloVI"]["mean"], statistics.mean(dataset_means))
        with source.open(newline="") as handle:
            source_rows = list(csv.DictReader(handle))
        members = {(row["dataset"], row["gene"]) for row in tables["membership.csv"]}
        pooled = [float(row["new_cbdir"]) for row in source_rows
                  if row["method"] == "veloVI" and (row["dataset"], row["gene"]) in members]
        self.assertGreater(abs(macro["veloVI"]["mean"] - statistics.mean(pooled)), 1e-4)


if __name__ == "__main__":
    unittest.main()
