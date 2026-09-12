"""Unit tests for the minimal signed-T ensemble contract."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from ensemble.core import (
    H5AD_SCHEMA_VERSION,
    METHOD_NAMES,
    VARIANT_MEAN4,
    VARIANT_NAMES,
    VARIANT_ROW_GEOMETRIC_MEDIAN4,
    aggregate_offdiagonal,
    align_method_values,
    complete_signed_operator,
    install_fixed_neighbors,
    load_config,
    load_fixed_neighbor_asset,
    project_with_scvelo,
    protocol_metadata,
    validate_signed_operator,
)
from ensemble.geometric_median import validate_row_geometric_median4
from ensemble.metrics import build_rankings, write_rankings


class EnsembleTests(unittest.TestCase):
    def test_config_has_four_methods_and_fourteen_datasets(self):
        config = load_config()
        self.assertEqual(tuple(config["methods"]), METHOD_NAMES)
        self.assertEqual(len(config["datasets"]), 14)
        self.assertEqual(set(config["datasets"]), {
            "bonemarrow",
            "capofibro",
            "dentategyrus",
            "dentate_lamanno",
            "forebrain",
            "gastrulation",
            "hematopoiesis",
            "mousebrain",
            "neuron",
            "organoid",
            "pancreas",
            "retina",
            "rpe1",
            "rpe1_gen",
        })
        self.assertEqual(set(VARIANT_NAMES), {
            VARIANT_MEAN4,
            VARIANT_ROW_GEOMETRIC_MEDIAN4,
        })

    def test_mean_and_geometric_median_are_permutation_invariant(self):
        values = np.asarray([
            [0.10, -0.40, 0.80],
            [0.50, 0.20, -0.10],
            [-0.20, 0.60, 0.30],
            [0.90, -0.80, 0.00],
        ])
        template = sparse.csr_matrix(
            ([1.0, 1.0, 1.0], ([0, 0, 1], [1, 2, 2])), shape=(3, 3)
        )
        for variant in VARIANT_NAMES:
            expected, detail = aggregate_offdiagonal(template, values, variant)
            permutation = np.asarray([2, 0, 3, 1])
            actual, permuted_detail = aggregate_offdiagonal(
                template, values[permutation], variant
            )
            np.testing.assert_allclose(actual, expected)
            self.assertEqual(detail["n_methods"], 4)
            self.assertEqual(permuted_detail["n_methods"], 4)

    def test_row_geometric_median_invariants(self):
        checks = validate_row_geometric_median4()
        self.assertTrue(all(checks.values()))

    def test_union_alignment_and_signed_l1_diagonal_completion(self):
        matrices = [
            sparse.csr_matrix(([0.2, -0.3], ([0, 0], [1, 2])), shape=(3, 3)),
            sparse.csr_matrix(([0.4], ([0], [1])), shape=(3, 3)),
            sparse.csr_matrix(([0.1], ([1], [2])), shape=(3, 3)),
            sparse.csr_matrix(([0.2], ([2], [0])), shape=(3, 3)),
        ]
        template, values = align_method_values(matrices)
        self.assertEqual(template.shape, (3, 3))
        self.assertEqual(values.shape[0], 4)
        operator, detail = complete_signed_operator(
            template, values.mean(axis=0)
        )
        validate_signed_operator(operator, "test operator")
        dense = operator.toarray()
        np.testing.assert_allclose(np.diag(dense), [0.775, 0.975, 0.95])
        self.assertEqual(detail["n_rows_rescaled"], 0)
        np.testing.assert_allclose(
            np.asarray(np.abs(operator).sum(axis=1)).ravel(), 1.0
        )

    def test_fixed_neighbor_asset_reorders_by_cell_name(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            n_obs = 30
            graph = sparse.eye(n_obs, format="csr")
            indices = np.asarray(
                [np.r_[row, np.delete(np.arange(n_obs), row)] for row in range(n_obs)]
            )
            np.savez_compressed(
                root / "toy.npz",
                obs_names=np.asarray([f"c{i}" for i in range(n_obs)]),
                indices=indices,
                distances_data=graph.data,
                distances_indices=graph.indices,
                distances_indptr=graph.indptr,
                distances_shape=np.asarray(graph.shape),
                connectivities_data=graph.data,
                connectivities_indices=graph.indices,
                connectivities_indptr=graph.indptr,
                connectivities_shape=np.asarray(graph.shape),
                metadata=np.asarray(json.dumps({
                    "dataset": "toy",
                    "parameters": {
                        "n_neighbors": 30,
                        "n_pcs": 30,
                        "use_rep": "X_pca",
                        "metric": "euclidean",
                        "random_state": 0,
                    },
                })),
            )
            asset = load_fixed_neighbor_asset(root / "toy.npz")
            data = ad.AnnData(np.zeros((n_obs, 1)))
            data.obs_names = [f"c{i}" for i in reversed(range(n_obs))]
            detail = install_fixed_neighbors(data, asset)
            mapped = data.uns["neighbors"]["indices"]
            np.testing.assert_array_equal(mapped[:, 0], np.arange(n_obs))
            self.assertEqual(mapped.shape, (n_obs, 30))
            self.assertEqual(detail["n_obs"], n_obs)
            self.assertEqual(
                data.uns["neighbors"]["params"]["n_neighbors"], 30
            )
            self.assertFalse(
                data.uns["neighbors"]["fixed_reference_neighbors"]
                ["reference_pca"]["installed"]
            )

    def test_protocol_metadata_is_frozen_and_portable(self):
        metadata = protocol_metadata()
        self.assertEqual(metadata["scvelo_version"], "0.3.4")
        self.assertEqual(
            metadata["transition_matrix"]["use_negative_cosines"], True
        )
        self.assertEqual(
            metadata["velocity_graph"], {"vkey": "velocity", "xkey": "Ms"}
        )
        self.assertEqual(metadata["fixed_neighbors"]["n_neighbors"], 30)

    def test_explicit_transition_projection(self):
        data = ad.AnnData(np.asarray([
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ]))
        data.layers["Ms"] = np.asarray([
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ])
        data.obsm["X_umap"] = np.asarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        data.obs["label"] = ["a", "b", "c"]
        operator = sparse.csr_matrix(np.asarray([
            [0.0, 0.75, 0.25],
            [0.25, 0.0, 0.75],
            [0.65, 0.35, 0.0],
        ]))
        result, invalid = project_with_scvelo(
            data,
            operator,
            variant=VARIANT_MEAN4,
            logical_inputs={"reference": {"locator": "velocity:toy"}},
        )
        self.assertLessEqual(invalid, 3)
        self.assertIn("ensemble_transition_signed", result.obsp)
        self.assertIn("ensemble_velocity_umap", result.obsm)
        self.assertEqual(
            result.uns["ensemble_signed_params"]["schema_version"],
            H5AD_SCHEMA_VERSION,
        )
        self.assertEqual(result.obsm["ensemble_velocity_umap"].shape, (3, 2))
        validate_signed_operator(
            result.obsp["ensemble_transition_signed"], "projected operator"
        )

    def test_no_deprecated_centerweighted_variant(self):
        repository_root = Path(__file__).resolve().parents[1]
        for relative in [
            "ensemble/core.py",
            "ensemble/metrics.py",
            "ensemble/config.json",
            "ensemble/README.md",
        ]:
            text = (repository_root / relative).read_text(encoding="utf-8")
            self.assertNotIn("Center" + "Weighted", text)
            self.assertNotIn("center" + "weighted", text)


    def test_ranking_strict_ties_and_variant_isolation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = {"datasets": {"bonemarrow": {}}}
            (root / "config.json").write_text(json.dumps(config), encoding="utf-8")
            with (root / "scores.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["dataset", "variant", "cbdir", "batc"])
                writer.writeheader()
                writer.writerows([
                    {"dataset": "bonemarrow", "variant": "signed_t_mean4", "cbdir": "0.5", "batc": "0.0"},
                    {"dataset": "bonemarrow", "variant": "signed_t_row_geometric_median4", "cbdir": "0.0", "batc": "0.5"},
                ])
            methods = [f"m{i:02d}" for i in range(20)]
            with (root / "cell.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["dataset", "method", "cbdir", "batc"])
                writer.writeheader()
                for i, method in enumerate(methods):
                    writer.writerow({"dataset": "bonemarrow", "method": method, "cbdir": 0.5 if i == 0 else (0.6 if i == 1 else 0.1), "batc": i / 20})
            rows = build_rankings([{"dataset": "bonemarrow", "variant": "signed_t_mean4", "cbdir": "0.5", "batc": "0.0"}, {"dataset": "bonemarrow", "variant": "signed_t_row_geometric_median4", "cbdir": "0.0", "batc": "0.5"}], list(csv.DictReader((root / "cell.csv").open())), {"bonemarrow": "bonemarrow"})
            mean_cbdir = next(r for r in rows if r["variant"] == "signed_t_mean4" and r["metric"] == "cbdir")
            center_cbdir = next(r for r in rows if r["variant"] == "signed_t_row_geometric_median4" and r["metric"] == "cbdir")
            self.assertEqual((mean_cbdir["rank"], mean_cbdir["tie_count"], mean_cbdir["methods_above"]), (2, 1, 1))
            self.assertEqual(center_cbdir["rank"], 21)
            self.assertEqual(sum(r["rank"] == 1 for r in rows if r["variant"] == "signed_t_mean4"), 0)

    def test_real_rankings_are_complete_and_byte_stable(self):
        repository_root = Path(__file__).resolve().parents[1]
        scores = repository_root / "ensemble/results/scores.csv"
        cell = repository_root / "evaluation/cell/results/cell_metrics.csv"
        config = repository_root / "ensemble/config.json"
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.csv"
            first_summary = Path(directory) / "first.md"
            second = Path(directory) / "second.csv"
            second_summary = Path(directory) / "second.md"
            write_rankings(scores, cell, config, first, first_summary)
            write_rankings(scores, cell, config, second, second_summary)
            first_bytes, second_bytes = first.read_bytes(), second.read_bytes()
            self.assertEqual(first_bytes, second_bytes)
            self.assertEqual(first_summary.read_bytes(), second_summary.read_bytes())
            with first.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 56)
            self.assertEqual({r["out_of"] for r in rows}, {"21"})
            self.assertEqual({(r["variant"], r["dataset"], r["metric"]) for r in rows}, {(v, d, m) for v in ("signed_t_mean4", "signed_t_row_geometric_median4") for d in load_config()["datasets"] for m in ("cbdir", "batc")})


if __name__ == "__main__":
    unittest.main()
