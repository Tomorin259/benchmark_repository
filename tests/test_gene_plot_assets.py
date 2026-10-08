"""Contracts for CSV-only summaries and registered native phase inputs."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

from evaluation.gene.plot_gene_assets import (
    CBDIR_CSV, COSSIM_CSV, ORDER, REAL_CONFIGS, ROOT, SIM_CONFIGS,
    load_scores, source_inputs, read_phase_inputs, sha256,
)


class GeneAssetTests(unittest.TestCase):
    def test_scores_keep_reference30_cbdir_and_exact_cossim_join(self):
        actual = load_scores()
        original = pd.read_csv(CBDIR_CSV)
        key = ["dataset", "method", "gene"]
        pd.testing.assert_series_equal(actual.set_index(key).new_cbdir.sort_index(),
                                       original.set_index(key).new_cbdir.sort_index())
        compact = pd.read_csv(COSSIM_CSV).rename(columns={
            "Dataset": "dataset", "Method": "method", "Gene": "gene", "CosSim": "cossim",
        })
        pd.testing.assert_series_equal(actual[actual.dataset.isin(SIM_CONFIGS)].set_index(key).cossim.sort_index(),
                                       compact.set_index(key).cossim.sort_index())
        self.assertTrue(actual.loc[actual.dataset.isin(REAL_CONFIGS), "cossim"].isna().all())
        # The compact CBDir column intentionally differs; it must never replace new_cbdir.
        paired = actual.merge(compact[key + ["CBDir"]], on=key)
        self.assertGreater(float(np.abs(paired.new_cbdir - paired.CBDir).max()), 0.01)

    def test_input_paths_and_layers_match_formal_score_provenance(self):
        scores = load_scores()
        for case in scores.dataset.unique():
            inputs = source_inputs(case, scores)
            self.assertEqual(len(inputs), len(ORDER))
            self.assertTrue(all(not Path(item["path"]).is_absolute() for item in inputs))
        corrupted = scores.copy()
        corrupted.loc[0, "velocity_u_key"] = "wrong_role"
        with self.assertRaisesRegex(ValueError, "disagrees"):
            source_inputs(corrupted.loc[0, "dataset"], corrupted)

    def test_duplicate_or_missing_cossim_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "cossim.csv"
            compact = pd.read_csv(COSSIM_CSV)
            pd.concat([compact, compact.iloc[:1]]).to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                load_scores(cossim_path=path)
            compact.iloc[1:].to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "different simulation gene keys"):
                load_scores(cossim_path=path)
            path.unlink()
            with self.assertRaises(FileNotFoundError):
                load_scores(cossim_path=path)

    def test_phase_inputs_preserve_native_coordinates_and_h5ad_bytes(self):
        from anndata import AnnData

        scores = load_scores()
        case, gene = "bonemarrow", "LMNA"
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            originals = {}
            for item in source_inputs(case, scores):
                path = root / item["path"]
                path.parent.mkdir(parents=True, exist_ok=True)
                data = AnnData(X=np.zeros((2, 1)),
                               obs=pd.DataFrame({"clusters": pd.Categorical(["a", "b"])}, index=["c0", "c1"]),
                               var=pd.DataFrame(index=[gene]))
                data.layers["Ms"] = np.array([[1.0], [3.0]])
                data.layers["Mu"] = np.array([[2.0], [4.0]])
                data.layers[item["velocity_s"]] = np.array([[np.nan], [1e12]])
                data.layers[item["velocity_u"]] = np.array([[np.inf], [-1e12]])
                data.uns["clusters_colors"] = ["#ff0000", "#0000ff"]
                data.write_h5ad(path)
                originals[path] = sha256(path)
            inputs, records = read_phase_inputs(case, scores, root, [gene])
            self.assertEqual(len(records), 9)
            for data in inputs:
                np.testing.assert_array_equal(data.layers["Ms"], [[1.0], [3.0]])
                np.testing.assert_array_equal(data.layers["Mu"], [[2.0], [4.0]])
                np.testing.assert_array_equal(data.layers["velocity"], [[0.0], [1e10]])
                np.testing.assert_array_equal(data.layers["velocity_u"], [[0.0], [-1e10]])
                self.assertNotIn("neighbors", data.uns)
            self.assertEqual(originals, {path: sha256(path) for path in originals})

    def test_summaries_run_without_h5ad_root_or_external_plot_scripts(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "panels"
            result = subprocess.run([
                sys.executable, str(ROOT / "evaluation/gene/plot_gene_assets.py"),
                "--summaries-only", "--cases", "bonemarrow", "murk", "--dpi", "30",
                "--output-root", str(output),
            ], cwd=folder, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output / "real/bonemarrow/summary_cbdir.png").is_file())
            self.assertTrue((output / "simulation/murk/summary_cbdir_cossim.png").is_file())
            provenance = json.loads((output / "provenance.json").read_text())
            self.assertEqual(provenance["inputs"], {})
            self.assertEqual(set(provenance["score_sources"]), {
                "evaluation/gene/results/gene_cbdir_reference30_all.csv",
                "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv",
            })
            self.assertFalse(list(output.rglob("*.txt")))
            # A second run cannot overwrite a successful export.
            again = subprocess.run(result.args, cwd=folder, capture_output=True, text=True)
            self.assertNotEqual(again.returncode, 0)
            self.assertIn("must be empty", again.stderr)


if __name__ == "__main__":
    unittest.main()
