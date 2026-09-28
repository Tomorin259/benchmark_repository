"""Cross-version roundoff must not hide genuinely stale Real-CBDir summaries."""

import csv
import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from evaluation.summary.refresh_real_cbdir_summary import check_summary


class RealCBDirSummaryTests(unittest.TestCase):
    def check_row(self, expected, actual, succeeds):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "summary.csv"
            fields = list(expected)
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow(actual)
            before = path.read_bytes()
            if succeeds:
                check_summary(path, fields, [expected])
            else:
                with self.assertRaises(SystemExit):
                    check_summary(path, fields, [expected])
            self.assertEqual(path.read_bytes(), before)

    def test_one_ulp_difference_is_accepted_without_rewriting(self):
        expected = dict(Metric="Real-CBDir", Method="DeepVelo_scv", Mean="0.0293057947653121",
                        Std="0.08813477152040747", N_dataset_or_system="7")
        actual = dict(expected, Std=repr(math.nextafter(float(expected["Std"]), math.inf)))
        self.check_row(expected, actual, succeeds=True)

    def test_stale_scores_and_nonfinite_values_are_rejected(self):
        expected = dict(Metric="Real-CBDir", Method="DeepVelo", Mean="0.12601342756275158")
        for value in ("0.1967575", "0.12611342756275158", "nan", "inf", "invalid"):
            with self.subTest(value=value):
                self.check_row(expected, dict(expected, Mean=value), succeeds=False)

    def test_counts_labels_and_other_metrics_remain_exact(self):
        expected = dict(Metric="Real-CBDir", Method="DeepVelo", Mean="0.1", N="7")
        self.check_row(expected, dict(expected, N="8"), succeeds=False)
        self.check_row(expected, dict(expected, Method="cellDancer"), succeeds=False)
        expected = dict(Metric="Kinetic-Spearman-alpha", Method="DeepVelo", Mean="0.3150070141443387")
        self.check_row(expected, dict(expected, Mean="0.3150070141443386"), succeeds=False)

    def test_formatted_wps_scores_remain_exact(self):
        expected = dict(Figure2_Method="DeepVelo", Real_CBDir_mean_plus_1="1.1197")
        self.check_row(expected, dict(expected, Real_CBDir_mean_plus_1="1.1198"), succeeds=False)

    def test_missing_rows_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "summary.csv"
            path.write_text("Metric,Method,Mean\n")
            with self.assertRaises(SystemExit):
                check_summary(path, ["Metric", "Method", "Mean"],
                              [dict(Metric="Real-CBDir", Method="DeepVelo", Mean="0.1")])


if __name__ == "__main__":
    unittest.main()
