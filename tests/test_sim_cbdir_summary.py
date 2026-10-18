"""Protect Figure 7 CBDir/CosSim against historical scores and pooled averages."""

import csv
import hashlib
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from evaluation.summary.refresh_real_cbdir_summary import check_summary
from evaluation.summary import refresh_sim_cbdir_summary as simulation


class SimulationSummaryTests(unittest.TestCase):
    def test_committed_summary_matches_reference30(self):
        subprocess.run([sys.executable, "evaluation/summary/refresh_sim_cbdir_summary.py"],
                       cwd=Path(__file__).resolve().parents[1], check=True)

    def test_refresh_preserves_current_precision(self):
        root = Path(__file__).resolve().parents[1]
        summaries = ["fig7_gene_metric_dataset_means.csv", "fig7_gene_metric_summary.csv", "fig7_gene_metric_wps_values.csv"]
        with tempfile.TemporaryDirectory() as directory:
            temporary = Path(directory)
            paths = ["evaluation/gene/results/gene_cbdir_reference30_all.csv",
                     "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv"]
            paths += ["evaluation/summary/results/figure7_summary/" + name for name in summaries]
            for relative in paths:
                destination = temporary / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(root / relative, destination)
            with (temporary / "artifacts.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["path", "sha256"], delimiter="\t")
                writer.writeheader()
                for relative in paths[2:]:
                    writer.writerow(dict(path=relative, sha256=hashlib.sha256((temporary / relative).read_bytes()).hexdigest()))
            result_path = temporary / "evaluation/summary/results/figure7_summary"
            before = {name: (result_path / name).read_bytes() for name in summaries}
            with patch.object(simulation, "ROOT", temporary), patch.object(simulation, "RESULTS", result_path):
                simulation.main(False)
                simulation.main(True)
            self.assertEqual(before, {name: (result_path / name).read_bytes() for name in summaries})

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

    def test_historical_velovi_cosine_and_drift_are_rejected(self):
        expected = dict(Metric="Sim-CosSim", Method="veloVI", Mean="0.07653622874170384")
        for value in ("0.039307", "0.076536228743"):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "summary.csv"
                with path.open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(expected))
                    writer.writeheader()
                    writer.writerow(dict(expected, Mean=value))
                with self.assertRaises(SystemExit):
                    check_summary(path, list(expected), [expected], metric=("Sim-CBDir", "Sim-CosSim"))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "summary.csv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(expected))
                writer.writeheader()
                writer.writerow(dict(expected, Mean="0.0765362287417"))
            check_summary(path, list(expected), [expected], metric="Sim-CosSim")


if __name__ == "__main__":
    unittest.main()
