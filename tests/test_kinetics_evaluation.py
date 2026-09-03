#!/usr/bin/env python3
"""Small synthetic tests for the kinetic evaluation boundary."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from evaluation.kinetics.metrics import (
    calculate_metrics,
    dynamo_rpe1_raw,
    organoid_raw_from_layer,
    process_organoid,
    process_rpe1,
    select_organoid_branch,
    tsvelo_rpe1_raw,
)


def write_h5ad(path, cells, positions, values):
    var_names = ["g1", "g2"]
    adata = AnnData(
        X=np.zeros((len(cells), len(var_names))),
        obs=pd.DataFrame(
            {"Cell_cycle_possition": positions},
            index=pd.Index(cells, dtype=str),
        ),
        var=pd.DataFrame(index=pd.Index(var_names, dtype=str)),
    )
    adata.layers["alpha"] = np.asarray(values, dtype=float)
    adata.write_h5ad(path)


class KineticsEvaluationTests(unittest.TestCase):
    def test_rpe1_processing_and_metrics(self):
        raw = pd.DataFrame(
            [[1.0, 2.0, 4.0, 8.0], [8.0, 4.0, 2.0, 1.0]],
            index=["g1", "g2"],
        )
        settings = {"smooth_window": 2, "epsilon": 1e-6}
        processed = process_rpe1(raw, settings, clip_negative=False)
        scores = calculate_metrics(processed, processed)

        self.assertEqual(processed.shape, raw.shape)
        self.assertEqual(list(scores["Gene"]), ["g1", "g2"])
        self.assertTrue(np.allclose(scores["Spearman_Corr"], 1.0))

    def test_dynamo_sparse_layer_is_binned_and_zero_genes_are_removed(self):
        matrix = sparse.csr_matrix(
            np.array([[1.0, 0.0, 2.0], [0.0, 0.0, 0.0], [3.0, 1.0, 0.0]])
        )
        with tempfile.TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            path = tmp_path / "dynamo.h5ad"
            with h5py.File(path, "w") as handle:
                layers = handle.create_group("layers")
                layer = layers.create_group("cell_wise_alpha")
                layer.attrs["encoding-type"] = "csr_matrix"
                layer.attrs["shape"] = matrix.shape
                layer.create_dataset("data", data=matrix.data)
                layer.create_dataset("indices", data=matrix.indices)
                layer.create_dataset("indptr", data=matrix.indptr)
                var = handle.create_group("var")
                var.create_dataset(
                    "_index", data=np.asarray(["g1", "g2", "g3"], dtype="S")
                )
                obs = handle.create_group("obs")
                obs.create_dataset(
                    "Cell_cycle_possition", data=np.asarray([0.0, 1.0, 2.0])
                )

            source = {
                "result": "dynamo.h5ad",
                "layer": "cell_wise_alpha",
                "position_key": "Cell_cycle_possition",
                "drop_all_zero_genes": True,
            }
            reference = pd.DataFrame(index=["g1", "g2"])
            settings = {"n_obs": 3, "positions": 3}
            raw = dynamo_rpe1_raw(source, reference, settings, tmp_path)

        self.assertEqual(list(raw.index), ["g1", "g2"])
        self.assertEqual(raw.shape, (2, 2))
        self.assertEqual(raw.loc["g1", 0], 0.5)
        self.assertEqual(raw.loc["g2", 1], 1.0)

    def test_tsvelo_lineages_average_duplicate_cells(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            write_h5ad(
                tmp_path / "lineage0.h5ad",
                ["cell1", "cell2"],
                [0, 1],
                [[1, 2], [3, 4]],
            )
            write_h5ad(
                tmp_path / "lineage1.h5ad",
                ["cell2", "cell3"],
                [1, 2],
                [[4, 5], [6, 7]],
            )
            source = {
                "results": ["lineage0.h5ad", "lineage1.h5ad"],
                "layer": "alpha",
                "position_key": "Cell_cycle_possition",
            }
            reference = pd.DataFrame(index=["g1", "g2"])
            settings = {"positions": 3}
            raw = tsvelo_rpe1_raw(source, reference, settings, tmp_path)

        self.assertTrue(np.allclose(raw.loc["g1"], [1.0, 3.5, 6.0]))
        self.assertTrue(np.allclose(raw.loc["g2"], [2.0, 4.5, 7.0]))

    def test_organoid_branch_selection_and_processing(self):
        adata = AnnData(
            X=np.zeros((4, 2)),
            obs=pd.DataFrame(
                {
                    "monocle_branch_id": ["1", "1", "2", "2"],
                    "som_cluster_id": ["4", "9", "1", "2"],
                    "monocle_pseudotime": [3.5, 3.0, 6.5, 7.5],
                },
                index=["c1", "c2", "c3", "c4"],
            ),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        adata.layers["alpha"] = np.asarray(
            [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0], [4.0, 8.0]]
        )
        settings = {
            "branch_key": "monocle_branch_id",
            "cluster_key": "som_cluster_id",
            "pseudotime_key": "monocle_pseudotime",
            "bins": 2,
            "smooth_window": 1,
            "baseline_bins": 1,
            "epsilon": 1e-6,
            "branches": {
                "1": {"clusters": ["4", "9"], "start": 3.5, "reverse_time": True}
            },
        }
        selected, aligned_time = select_organoid_branch(adata, settings, "1")
        raw = organoid_raw_from_layer(selected, aligned_time, "alpha", settings)
        processed = process_organoid(raw, settings)

        self.assertEqual(list(selected.obs_names), ["c1", "c2"])
        self.assertTrue(np.allclose(aligned_time.to_numpy(), [0.0, 0.5]))
        self.assertEqual(processed.shape, (2, 2))


if __name__ == "__main__":
    unittest.main()
