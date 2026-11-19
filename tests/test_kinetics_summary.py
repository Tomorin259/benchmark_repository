"""Ensure integrated kinetics summaries follow the current full-gene evaluation."""

import csv
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

    def test_figure7_uses_both_systems_on_expanded_reference_genes(self):
        folder = ROOT / "evaluation/summary/results/figure7_summary"
        with (folder / "fig7_gene_metric_summary.csv").open() as stream:
            alpha = {row["Method"]: row for row in csv.DictReader(stream)
                     if row["Metric"] == "Kinetic-Spearman-alpha"}
        with (folder / "fig7_gene_metric_wps_values.csv").open() as stream:
            display = {row["Figure7_Method"]: row for row in csv.DictReader(stream)}
        for method, mean, count in (("TSvelo", 0.3263157762832308, 105),
                                    ("Dynamo(label data)", 0.6446651075033082, 513)):
            with self.subTest(method=method):
                self.assertEqual(int(alpha[method]["N_dataset_or_system"]), 2)
                self.assertEqual(int(alpha[method]["N_raw_genes_or_rows"]), count)
                self.assertAlmostEqual(float(alpha[method]["Mean"]), mean, places=12)
                self.assertEqual(int(display[method]["Kinetic_Spearman_alpha_n_gene_rows"]), count)
                self.assertAlmostEqual(float(display[method]["Kinetic_Spearman_alpha_mean_plus_1"]),
                                       mean + 1, places=12)

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
