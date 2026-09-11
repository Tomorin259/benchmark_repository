#!/usr/bin/env python3
"""Small behavior tests for the cell-level evaluation boundary."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from evaluation.cell.metrics import consistency_scores, resolve_artifact
from evaluation.cell.reference_neighbors import (
    DEFAULT_REFERENCE_NEIGHBORS_ROOT,
    REFERENCE_PARAMS,
    apply_reference_neighbors,
    load_reference_neighbors,
)
from evaluation.cell.run_cell_metrics import DEFAULT_INPUTS
from evaluation.cell.validate_results import validate_results


class CellEvaluationTests(unittest.TestCase):
    def test_cli_defaults_are_repository_anchored(self):
        repository_root = Path(__file__).resolve().parents[1]
        self.assertEqual(
            DEFAULT_INPUTS, repository_root / "evaluation" / "cell" / "inputs.tsv"
        )

    def test_cell_inputs_embed_deepvelo_role_contract(self):
        inputs_path = Path(__file__).resolve().parents[1] / "evaluation" / "cell" / "inputs.tsv"
        with inputs_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, ["config_id", "operational_id", "dataset_id", "evaluation_dataset", "result_path", "velocity_layer_key", "source_embedding_key", "velocity_embedding_key", "evaluation_embedding", "alignment_mode", "n_reference", "n_method", "n_common", "source_embedding_max_error", "affine_relative_rmse", "source_nonfinite_velocity_values", "aligned_nonfinite_velocity_values", "file_size_bytes", "sha256", "evaluation_role", "pred_unspliced"])
            rows = list(reader)
        self.assertEqual(len(rows), 280)
        self.assertEqual({row["evaluation_role"] for row in rows}, {"cell"})
        deepvelo = [row for row in rows if row["operational_id"] == "deepvelo"]
        self.assertEqual(len(deepvelo), 14)
        self.assertEqual({row["pred_unspliced"] for row in deepvelo}, {"false"})
        self.assertEqual(
            {row["pred_unspliced"] for row in rows if row["operational_id"] != "deepvelo"},
            {""},
        )

    def test_canonical_bundle_validates_tables_and_reference_assets(self):
        validate_results()

    def test_reference_contract_comes_from_cell_config(self):
        self.assertEqual(
            DEFAULT_REFERENCE_NEIGHBORS_ROOT,
            Path(__file__).resolve().parents[1] / "evaluation" / "cell" / "reference_neighbors",
        )
        self.assertEqual(REFERENCE_PARAMS["n_neighbors"], 30)
        self.assertEqual(REFERENCE_PARAMS["n_pcs"], 30)
        self.assertEqual(REFERENCE_PARAMS["mode_neighbors"], "distances")
        self.assertEqual(REFERENCE_PARAMS["n_recurse_neighbors"], 1)

    def test_reference_neighbor_default_asset_is_repository_local(self):
        reference = load_reference_neighbors("pancreas")
        self.assertEqual(reference["indices"].shape, (3696, 30))
        self.assertEqual(Path(reference["path"]).resolve().parent, DEFAULT_REFERENCE_NEIGHBORS_ROOT)
        self.assertFalse(str(reference["metadata"].get("source", "")).startswith("/"))

    def test_consistency_uses_neighbor_cosines_and_skips_self(self):
        velocity = np.asarray(
            [[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float
        )
        neighbors = [[0, 1, 2], [1, 0, 2], [2, 0, 1]]
        scores = consistency_scores(velocity, neighbors, n_neighbors=2)
        self.assertTrue(np.allclose(scores, [0.5, 0.5, 0.0]))

    def test_resolve_artifact_rejects_absolute_and_parent_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "result.h5ad").write_bytes(b"test")
            self.assertEqual(
                resolve_artifact(root, "result.h5ad"), root / "result.h5ad"
            )
            with self.assertRaises(ValueError):
                resolve_artifact(root, "/tmp/result.h5ad")
            with self.assertRaises(ValueError):
                resolve_artifact(root, "../result.h5ad")

    def test_reference_neighbors_align_by_cell_name(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            graph = sparse.csr_matrix(
                np.asarray(
                    [[1.0, 0.5, 0.0], [0.5, 1.0, 0.25], [0.0, 0.25, 1.0]]
                )
            )
            np.savez_compressed(
                root / "toy.npz",
                obs_names=np.asarray(["c0", "c1", "c2"]),
                indices=np.asarray([[0, 1, 2], [1, 0, 2], [2, 1, 0]]),
                distances_data=graph.data,
                distances_indices=graph.indices,
                distances_indptr=graph.indptr,
                distances_shape=np.asarray(graph.shape),
                connectivities_data=graph.data,
                connectivities_indices=graph.indices,
                connectivities_indptr=graph.indptr,
                connectivities_shape=np.asarray(graph.shape),
                metadata=np.asarray('{"dataset": "toy"}'),
            )
            adata = ad.AnnData(np.zeros((3, 1)))
            adata.obs_names = ["c2", "c0", "c1"]
            applied, detail = apply_reference_neighbors(adata, "toy", root)
            self.assertTrue(applied, detail)
            self.assertEqual(adata.uns["neighbors"]["params"]["n_neighbors"], 30)
            self.assertTrue(
                np.array_equal(
                    adata.uns["neighbors"]["indices"],
                    np.asarray([[0, 2, 1], [1, 2, 0], [2, 1, 0]]),
                )
            )


if __name__ == "__main__":
    unittest.main()
