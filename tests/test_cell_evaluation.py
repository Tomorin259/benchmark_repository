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

from evaluation.cell.batc import compute_batc
from evaluation.cell.metrics import consistency_scores, resolve_artifact
from evaluation.cell.reference_neighbors import (
    DEFAULT_REFERENCE_NEIGHBORS_ROOT,
    REFERENCE_PARAMS,
    apply_reference_neighbors,
    load_reference_neighbors,
)
from evaluation.cell.run_cell_metrics import DEFAULT_INPUTS, extended_edge_records
from evaluation.cell.validate_results import validate_results


class CellEvaluationTests(unittest.TestCase):
    def test_optional_batc_edge_audit_preserves_branch_max_score(self):
        rng = np.random.default_rng(42)
        coordinates = np.concatenate([
            np.column_stack((np.linspace(0, 0.9, 30), 0.01 * rng.normal(size=30))),
            np.column_stack((np.linspace(1.1, 2, 30), 0.03 + 0.01 * rng.normal(size=30))),
            np.column_stack((
                np.linspace(1.1, 2, 30),
                0.2 + np.linspace(0, 1, 30) + 0.01 * rng.normal(size=30),
            )),
        ])
        adata = ad.AnnData(np.zeros((len(coordinates), 1)))
        adata.obs["clusters"] = ["A"] * 30 + ["B"] * 30 + ["C"] * 30
        adata.obsm["X_umap"] = coordinates
        adata.obsm["velocity_umap"] = np.tile([1.0, 0.0], (len(coordinates), 1))
        options = dict(n_bins=5, min_per_bin=2, n_samples=80, store_in_adata=False)
        edges = [("A", "B"), ("A", "C")]
        score = compute_batc(adata, edges, **options)
        detailed_score, detail = compute_batc(
            adata, edges, return_details=True, **options
        )
        self.assertEqual(score, detailed_score)
        self.assertEqual(score, detail["overall"])
        self.assertEqual(
            sum(item["branch_winner_source_cells"] for item in detail["edges"]), 30
        )
        self.assertEqual(
            {(item["source"], item["target"]) for item in detail["edges"]},
            set(edges),
        )

    def test_extended_edge_audit_counts_nonfinite_native_values(self):
        indicators = {"metrics": {
            "cbdir_edge_raw": {"A->B": [0.5, None], "B->C": []},
            "batc_edge_details": [{
                "source": "A", "target": "B",
                "candidate_mean_source_cells": 0.25,
                "source_cells": 2,
                "candidate_nonfinite_source_cells": 1,
                "branch_winner_source_cells": 2,
            }],
        }}
        records = extended_edge_records("toy", "method", indicators)
        self.assertEqual(len(records), 3)
        self.assertIsNone(records[0]["score"])
        self.assertEqual(records[0]["raw_nonfinite_source_cells"], 1)
        self.assertEqual(records[0]["evaluable_source_cells"], 2)
        self.assertIsNone(records[1]["score"])
        self.assertEqual(records[2]["metric"], "batc_candidate_source_mean")
        self.assertEqual(records[2]["branch_winner_source_cells"], 2)

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

    def test_gene_matched_cell_inputs_follow_gene_result_sources(self):
        root = Path(__file__).resolve().parents[1]
        with (root / "methods/configurations.tsv").open(newline="", encoding="utf-8") as handle:
            names = {
                row["config_id"]: row["manuscript_name"]
                for row in csv.DictReader(handle, delimiter="\t")
            }
        with (root / "evaluation/gene/results/gene_cbdir_reference30_all.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            gene_rows = list(csv.DictReader(handle))
        with (root / "evaluation/cell/gene_matched/inputs.tsv").open(
            newline="", encoding="utf-8"
        ) as handle:
            inputs = list(csv.DictReader(handle, delimiter="\t"))
        with (root / "evaluation/cell/inputs.tsv").open(
            newline="", encoding="utf-8"
        ) as handle:
            canonical = list(csv.DictReader(handle, delimiter="\t"))

        reference = [row for row in inputs if row["evaluation_role"] == "reference_only"]
        scored = [row for row in inputs if row["evaluation_role"] == "gene_matched_cell"]
        self.assertEqual(len(reference), 7)
        self.assertEqual(len(scored), 21)
        self.assertEqual(len({(r["evaluation_dataset"], r["operational_id"]) for r in inputs}), 28)
        for row in reference:
            matches = [
                item for item in canonical
                if item["evaluation_dataset"] == row["evaluation_dataset"]
                and item["operational_id"] == "scvelo_deterministic"
            ]
            self.assertEqual(len(matches), 1)
            self.assertEqual(row["result_path"], matches[0]["result_path"])
        for row in scored:
            matches = {
                (item["source_h5ad"], item["velocity_s_key"])
                for item in gene_rows
                if item["status"] == "ok"
                and item["dataset"].casefold() == row["dataset_id"].casefold()
                and item["method"] == names[row["config_id"]]
            }
            self.assertEqual(len(matches), 1)
            source, velocity_key = matches.pop()
            self.assertEqual(
                row["result_path"],
                "velocity:" + source.replace("artifact://velocity/", "", 1),
            )
            self.assertEqual(row["velocity_layer_key"], velocity_key)
            if row["operational_id"] == "deepvelo":
                self.assertEqual(row["pred_unspliced"], "true")

    def test_gene_matched_cell_results_are_complete_and_audited(self):
        root = Path(__file__).resolve().parents[1]
        result_root = root / "evaluation/cell/gene_matched/results"

        def rows(path):
            with path.open(newline="", encoding="utf-8") as handle:
                return list(csv.DictReader(handle))

        with (root / "evaluation/cell/gene_matched/inputs.tsv").open(
            newline="", encoding="utf-8"
        ) as handle:
            manifest = list(csv.DictReader(handle, delimiter="\t"))
        sources = {
            (row["evaluation_dataset"], row["operational_id"]): row["result_path"]
            for row in manifest if row["evaluation_role"] == "gene_matched_cell"
        }
        primary = rows(result_root / "cell_metrics.csv")
        supplementary = rows(result_root / "cell_metrics_supplementary.csv")
        audit = rows(result_root / "cell_metrics_audit.csv")
        parity = rows(result_root / "tivelo_parity.csv")
        for table in (primary, supplementary, audit):
            self.assertEqual(len(table), 21)
            self.assertEqual({(r["dataset"], r["method"]) for r in table}, set(sources))
        for row in audit:
            key = row["dataset"], row["method"]
            self.assertEqual(row["result_path"], sources[key])
            self.assertEqual(len(row["result_sha256"]), 64)
            self.assertEqual(row["status"], "ok")
        self.assertEqual(len(parity), 42)
        self.assertTrue(all(row["raw_values_match"] == "True" for row in parity))
        self.assertTrue(all(float(row["aggregate_abs_difference"]) <= 1e-7 for row in parity))

        canonical = {
            (row["dataset"], row["method"]): row
            for row in rows(root / "evaluation/cell/results/cell_metrics.csv")
        }
        for row in primary:
            if row["method"] == "unitvelo":
                for metric in ("consistency", "iccoh", "cbdir", "transcosine", "batc"):
                    self.assertEqual(float(row[metric]), float(canonical[row["dataset"], row["method"]][metric]))

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
