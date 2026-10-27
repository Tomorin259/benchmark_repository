"""Native MultiVelo velocity, mask and projection must describe one output."""
import unittest
import tempfile
import csv
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from evaluation.cell.prepare_cell_metrics import prepare_velocity, install_native_2d_velocity
from evaluation.cell.run_cell_metrics import write_forebrain_cbdir_audit


class NativeVelocityContract(unittest.TestCase):
    def test_all_scored_input_maps_bind_the_worker_source_digest(self):
        root = Path(__file__).resolve().parents[1] / "evaluation/cell"
        for scope in ["", "gene_matched", "simulation"]:
            with (root / scope / "inputs.tsv").open() as handle:
                inputs = list(csv.DictReader(handle, delimiter="\t"))
            sources = {(r["evaluation_dataset"], r["operational_id"]): r for r in inputs}
            with (root / scope / "results/cell_metrics_audit.csv").open() as handle:
                for audit in csv.DictReader(handle):
                    row = sources[audit["dataset"], audit["method"]]
                    self.assertEqual(row["sha256"], audit["result_sha256"])
                    self.assertEqual(row["result_path"], audit["result_path"])

    def test_simulation_only_scope_does_not_require_a_forebrain_method_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_forebrain_cbdir_audit([], root, {"datasets": {"murk": {}}})
            self.assertEqual(list(root.iterdir()), [])

    def test_native_selection_updates_values_mask_and_mode_together(self):
        obj = ad.AnnData(np.zeros((3, 2)))
        native = np.asarray([[1., 2.], [3., 4.], [5., 6.]])
        obj.layers["velocity"] = np.full((3, 2), -99.)
        obj.layers["velo_s_norm"] = native.copy()
        obj.var["velocity_genes"] = [False, True]
        obj.var["velo_s_norm_genes"] = [True, False]
        obj.uns["velocity_params"] = {"mode": "stochastic"}
        obj.uns["velo_s_norm_params"] = {"mode": "dynamical"}
        prepare_velocity(obj, "velo_s_norm")
        np.testing.assert_array_equal(obj.layers["velocity"], native)
        np.testing.assert_array_equal(obj.var["velocity_genes"], [True, False])
        self.assertEqual(obj.uns["velocity_params"]["mode"], "dynamical")
        np.testing.assert_array_equal(obj.layers["velo_s_norm"], native)

    def test_native_selection_requires_matching_metadata(self):
        for missing in ["mask", "params"]:
            obj = ad.AnnData(np.zeros((2, 2)))
            obj.layers["velo_s_norm"] = np.ones((2, 2))
            if missing != "mask": obj.var["velo_s_norm_genes"] = [True, True]
            if missing != "params": obj.uns["velo_s_norm_params"] = {"mode": "dynamical"}
            with self.assertRaises(KeyError): prepare_velocity(obj, "velo_s_norm")

    def test_native_velocity_rejects_legacy_projection(self):
        obj = ad.AnnData(np.zeros((2, 2)))
        obj.obsm["X_umap"] = np.zeros((2, 2))
        obj.obsm["velocity_umap"] = np.ones((2, 2))
        row = pd.Series({"velocity_layer_key": "velo_s_norm", "velocity_embedding_key": "velocity_umap", "alignment_mode": "direct"})
        with self.assertRaisesRegex(ValueError, "matching projection"):
            install_native_2d_velocity(obj, row, "pancreas", obj.obsm["X_umap"])


if __name__ == "__main__": unittest.main()
