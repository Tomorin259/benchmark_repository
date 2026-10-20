#!/usr/bin/env python3
"""Small behavior tests for the cell-level evaluation boundary."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
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
from evaluation.cell.run_cell_metrics import (
    DEFAULT_INPUTS, evaluation_protocol, extended_edge_records, pair_is_current,
    write_tables,
)
from evaluation.cell.reference_embeddings import (
    array_sha256, load_forebrain_frames, PCA_COORDINATE_SHA256,
    UMAP_COORDINATE_SHA256,
)
from evaluation.cell.validate_results import _validate_forebrain_cbdir, validate_results


def forebrain_pair_fixture():
    """Reconstruct the 20 frozen pair records without external H5AD inputs."""
    root = DEFAULT_INPUTS.parent
    config = json.loads((root / "config.json").read_text())

    def table(name, delimiter=","):
        with (root / name).open(encoding="utf-8", newline="") as handle:
            return [row for row in csv.DictReader(handle, delimiter=delimiter)
                    if row.get("dataset", row.get("evaluation_dataset")) == "forebrain"]

    inputs = {row["operational_id"]: row for row in table("inputs.tsv", "\t")}
    primary = {row["method"]: row for row in table("results/cell_metrics.csv")}
    supplementary = {row["method"]: row for row in table("results/cell_metrics_supplementary.csv")}
    audit = {row["method"]: row for row in table("results/cell_metrics_audit.csv")}
    edges = table("results/forebrain_cbdir_pca_edge_details.csv")
    frozen_manifest = json.loads((root / "results/forebrain_cbdir_pca_manifest.json").read_text())
    file_hashes = {
        name: hashlib.sha256((root / relative).read_bytes()).hexdigest()
        for name, relative in {
            "config": "config.json", "inputs": "inputs.tsv",
            "prepare_script": "prepare_cell_metrics.py",
            "official_script": "score_tivelo_metrics.py",
            "frame_loader": "reference_embeddings.py", "metrics": "metrics.py",
            "batc": "batc.py", "reference_neighbors_code": "reference_neighbors.py",
            "reference_neighbors": "reference_neighbors/forebrain.npz",
            "frame_manifest": "reference_embeddings/manifest.json",
            "umap_frame": "reference_embeddings/forebrain_umap.npz",
        }.items()
    }
    file_hashes["reference_moments"] = frozen_manifest["reference_moments_sha256"]
    protocols = {}
    pairs = {}
    for method in config["method_order"]:
        protocol = evaluation_protocol(inputs[method], config, file_hashes)
        protocols["forebrain", method] = protocol
        pairs[method] = {
            "dataset": "forebrain", "method": method, "status": "ok",
            "evaluation_protocol": protocol,
            "metrics": {
                **{key: float(primary[method][key]) for key in config["result_tables"]["primary"][2:]},
                "highdim_cbdir": float(supplementary[method]["highdim_cbdir"]),
            },
            "counts": {key: int(value) for key, value in supplementary[method].items()
                       if key.endswith(("_edges", "_cells"))},
            "metadata": {
                **audit[method], "invalid_velocity_values_replaced_with_zero": 0,
                "scvelo_version": "0.2.5",
                "metric_frames": {
                    "cbdir": {"space": "pca", "dimensions": 30,
                              "coordinate_sha256": PCA_COORDINATE_SHA256},
                    "batc": {"space": "umap", "dimensions": 2,
                             "coordinate_sha256": UMAP_COORDINATE_SHA256},
                },
            },
            "edge_details": [], "tivelo_parity": [],
            "extended_edge_details": [{
                "dataset": "forebrain", "method": method, "metric": "cbdir",
                "source": row["source"], "target": row["target"],
                "score": float(row["score"]) if row["score"] else None,
                "evaluable_source_cells": int(row["evaluable_source_cells"]),
                "raw_nonfinite_source_cells": 0, "branch_winner_source_cells": "",
            } for row in edges if row["method"] == method],
        }
    return config, pairs, protocols


class CellEvaluationTests(unittest.TestCase):
    def test_forebrain_frames_preserve_the_frozen_coordinate_bytes(self):
        frames = load_forebrain_frames()
        self.assertEqual(frames["X_pca"].shape, (1720, 30))
        self.assertEqual(frames["X_umap"].shape, (1720, 2))
        self.assertEqual(frames["X_pca"].dtype, np.dtype("float32"))
        self.assertEqual(array_sha256(frames["X_pca"]), PCA_COORDINATE_SHA256)
        self.assertEqual(array_sha256(frames["X_umap"]), UMAP_COORDINATE_SHA256)

    def test_forebrain_frames_reject_a_different_cell_order(self):
        root = Path(__file__).resolve().parents[1] / "evaluation/cell"
        frames = load_forebrain_frames()
        with tempfile.TemporaryDirectory() as directory:
            changed = Path(directory)
            np.savez_compressed(changed / "forebrain_umap.npz",
                obs_names=frames["obs_names"][::-1], X_umap=frames["X_umap"])
            manifest = json.loads((root / "reference_embeddings/manifest.json").read_text())
            manifest["datasets"]["forebrain"]["umap_file_sha256"] = hashlib.sha256(
                (changed / "forebrain_umap.npz").read_bytes()).hexdigest()
            (changed / "manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "cell orders differ"):
                load_forebrain_frames(root / "reference_neighbors", changed)

    def test_resume_rejects_old_umap_pairs_and_changed_protocols(self):
        config = {"datasets": {"forebrain": {"cbdir": {
            "embedding": "pca", "dimensions": 30,
            "coordinate_sha256": PCA_COORDINATE_SHA256,
        }}}}
        row = {"config_id": "C01", "evaluation_dataset": "forebrain", "sha256": "source",
               "velocity_layer_key": "velocity"}
        protocol = evaluation_protocol(row, config, {"config": "new"})
        old = {"status": "ok", "extended_edge_details": []}
        self.assertFalse(pair_is_current(old, protocol))
        current = {**old, "evaluation_protocol": protocol}
        self.assertTrue(pair_is_current(current, protocol))
        changed = evaluation_protocol(row, config, {"config": "different"})
        self.assertFalse(pair_is_current(current, changed))
        self.assertFalse(pair_is_current(current, None))

    def test_dataset_cli_does_not_require_unselected_reference_materials(self):
        root = DEFAULT_INPUTS.parent
        with DEFAULT_INPUTS.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames
            rows = [row for row in reader if row["config_id"] == "C01"
                    and row["evaluation_dataset"] in {"forebrain", "pancreas"}]
        selected = next(row for row in rows if row["evaluation_dataset"] == "forebrain")
        official = {"highdim_cbdir": 0.5, "transcosine": 0.6,
                    "highdim_cbdir_raw": {"0->1": [0.5]},
                    "transcosine_raw": {"0->1": [0.6]}}
        indicators = {
            "metadata": {"result_path": selected["result_path"],
                         "result_sha256": selected["sha256"],
                         "velocity_layer_key": selected["velocity_layer_key"],
                         "invalid_velocity_values_replaced_with_zero": 0},
            "metrics": {**official, "consistency": 0.1, "iccoh": 0.2,
                        "cbdir": 0.3, "batc": 0.4,
                        "cbdir_edge_raw": {"0->1": [0.3]}, "batc_edge_details": []},
        }
        with tempfile.TemporaryDirectory() as directory:
            temporary = Path(directory)
            inputs = temporary / "inputs.tsv"
            with inputs.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            moments = temporary / "moments"
            neighbors = temporary / "neighbors"
            moments.mkdir()
            neighbors.mkdir()
            (moments / "forebrain.npz").write_bytes(b"stub reference moments")
            shutil.copyfile(root / "reference_neighbors/forebrain.npz", neighbors / "forebrain.npz")
            prepare = temporary / "prepare.py"
            score = temporary / "score.py"
            prelude = (
                "import json, sys\nfrom pathlib import Path\n"
                "def option(name): return Path(sys.argv[sys.argv.index(name) + 1])\n"
            )
            prepare.write_text(
                prelude + f"payload = json.loads({json.dumps(json.dumps(indicators))})\n"
                "option('--prepared-output').write_bytes(b'stub prepared input')\n"
                "option('--indicators-output').write_text(json.dumps(payload))\n"
            )
            score.write_text(
                prelude + f"payload = json.loads({json.dumps(json.dumps(official))})\n"
                "option('--output').write_text(json.dumps(payload))\n"
            )
            command = [
                sys.executable, str(root / "run_cell_metrics.py"),
                "--artifact-root", str(temporary), "--inputs", str(inputs),
                "--reference-neighbors-root", str(neighbors),
                "--reference-moments-root", str(moments),
                "--output-root", str(temporary / "output"),
                "--indicators-python", sys.executable, "--tivelo-python", sys.executable,
                "--prepare-script", str(prepare), "--official-script", str(score),
                "--dataset", "forebrain",
            ]
            completed = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            with (temporary / "output/results/cell_metrics.csv").open(newline="") as handle:
                result = list(csv.DictReader(handle))
            self.assertEqual([(row["dataset"], row["method"]) for row in result],
                             [("forebrain", "scvelo_deterministic")])
            self.assertFalse((moments / "pancreas.npz").exists())
            self.assertFalse((neighbors / "pancreas.npz").exists())
            canonical = temporary / "canonical_cell"
            canonical.mkdir()
            shutil.copyfile(root / "run_cell_metrics.py", canonical / "run_cell_metrics.py")
            command[1] = str(canonical / "run_cell_metrics.py")
            command[command.index("--output-root") + 1] = str(canonical)
            command.extend(["--config", str(root / "config.json")])
            protected = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(protected.returncode, 0)
            self.assertIn("A partial run must use a separate --output-root", protected.stderr)
            self.assertFalse((canonical / "results").exists())
            self.assertFalse((canonical / ".work").exists())

    def test_forebrain_rerun_writes_valid_complete_audit_and_excludes_zero_empty_edges(self):
        config, pairs, protocols = forebrain_pair_fixture()
        root = DEFAULT_INPUTS.parent
        with tempfile.TemporaryDirectory() as directory:
            pair_root = Path(directory) / "pairs"
            result_root = Path(directory) / "results"
            pair_root.mkdir()
            for case in ("frozen", "zero_and_empty"):
                with self.subTest(case=case):
                    if case == "zero_and_empty":
                        payload = pairs[config["method_order"][0]]
                        edges = payload["extended_edge_details"]
                        edges[0]["score"] = 0.0
                        edges[1].update(score=None, evaluable_source_cells=0)
                        payload["metrics"]["cbdir"] = np.mean([edge["score"] for edge in edges[2:]])
                    for method, payload in pairs.items():
                        (pair_root / f"forebrain__{method}.json").write_text(json.dumps(payload))
                    write_tables(pair_root, result_root, config, protocols)
                    with (result_root / "cell_metrics.csv").open(newline="") as handle:
                        primary = list(csv.DictReader(handle))
                    _validate_forebrain_cbdir(config, result_root,
                        root / "reference_neighbors", root / "reference_embeddings", primary)
                    manifest = json.loads((result_root / "forebrain_cbdir_pca_manifest.json").read_text())
                    self.assertEqual(manifest["edge_rows"], 120)
                    self.assertEqual(manifest["pca_coordinate_sha256"], PCA_COORDINATE_SHA256)
                    self.assertEqual(manifest["umap_coordinate_sha256"], UMAP_COORDINATE_SHA256)
                    self.assertEqual(manifest["evaluation_files"],
                                     protocols["forebrain", config["method_order"][0]]["files"])
                    self.assertEqual({row["method"]: row["sha256"] for row in manifest["inputs"]},
                        {method: payload["metadata"]["result_sha256"] for method, payload in pairs.items()})
                    with (result_root / "forebrain_cbdir_pca_edge_details.csv").open(newline="") as handle:
                        audit = list(csv.DictReader(handle))
                    self.assertEqual(len(audit), 120)
                    if case == "zero_and_empty":
                        method_edges = [row for row in audit if row["method"] == config["method_order"][0]]
                        self.assertEqual([row["included_in_mean"] for row in method_edges[:2]],
                                         ["False", "False"])

    def test_partial_forebrain_rerun_does_not_leave_complete_audit(self):
        config, pairs, protocols = forebrain_pair_fixture()
        with tempfile.TemporaryDirectory() as directory:
            pair_root = Path(directory) / "pairs"
            result_root = Path(directory) / "results"
            pair_root.mkdir()
            for method, payload in pairs.items():
                (pair_root / f"forebrain__{method}.json").write_text(json.dumps(payload))
            write_tables(pair_root, result_root, config, protocols)
            manifest = result_root / "forebrain_cbdir_pca_manifest.json"
            self.assertTrue(manifest.exists())
            (pair_root / f"forebrain__{config['method_order'][-1]}.json").unlink()
            write_tables(pair_root, result_root, config, protocols)
            self.assertFalse(manifest.exists())
            self.assertFalse((result_root / "forebrain_cbdir_pca_edge_details.csv").exists())

    def test_unitvelo_pca_cbdir_uses_dimensions_beyond_the_display(self):
        try:
            import unitvelo
        except ImportError:
            self.skipTest("UnitVelo is available in the indicators environment")
        adata = ad.AnnData(np.zeros((2, 1)))
        adata.obs["clusters"] = ["A", "B"]
        adata.uns["neighbors"] = {"indices": np.asarray([[0, 1], [1, 0]])}
        adata.obsm["X_pca"] = np.asarray([[0., 0., 0.], [0., 0., 1.]])
        adata.obsm["velocity_umap"] = np.asarray([[1., 0.], [1., 0.]])
        adata.obsm["cbdir_pca_velocity"] = np.asarray([[0., 0., 1.], [0., 0., 1.]])
        raw = unitvelo.eval_utils.cross_boundary_correctness(
            adata, "clusters", "cbdir_pca_velocity", [("A", "B")],
            return_raw=True, x_emb="X_pca",
        )
        np.testing.assert_allclose(raw[("A", "B")], [1.0])

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
        self.assertEqual(score, compute_batc(
            adata, edges, embedding_key="X_umap", velocity_key="velocity_umap", **options
        ))
        adata.obsm["velocity"] = np.zeros_like(coordinates)
        self.assertEqual(score, compute_batc(adata, edges, **options))
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

    def test_simulation_cell_results_match_gene_sources_and_reference(self):
        root = Path(__file__).resolve().parents[1]
        sim = root / "evaluation/cell/simulation"

        def csv_rows(path, delimiter=","):
            with path.open(newline="", encoding="utf-8") as handle:
                return list(csv.DictReader(handle, delimiter=delimiter))

        inputs = csv_rows(sim / "inputs.tsv", "\t")
        primary = csv_rows(sim / "results/cell_metrics.csv")
        supplementary = csv_rows(sim / "results/cell_metrics_supplementary.csv")
        audit = csv_rows(sim / "results/cell_metrics_audit.csv")
        parity = csv_rows(sim / "results/tivelo_parity.csv")
        paired = csv_rows(sim / "results/cross_level_summary.csv")
        gene = csv_rows(root / "evaluation/gene/results/gene_cbdir_reference30_all.csv")
        cossim = csv_rows(root / "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv")
        cossim_by_key = {
            (row["Dataset"], row["Gene"], row["Method"]): float(row["CosSim"])
            for row in cossim
        }
        self.assertEqual(len(cossim_by_key), len(cossim))
        registry = csv_rows(root / "methods/configurations.tsv", "\t")
        names = {row["config_id"]: row["manuscript_name"] for row in registry}
        gene_names = {**names, "C04": "MultiVelo", "C19": "Graphvelo_scv"}

        keys = {(row["evaluation_dataset"], row["operational_id"]) for row in inputs}
        self.assertEqual(len(inputs), len(keys))
        self.assertEqual({dataset for dataset, _ in keys}, {"murk", "branching", "all"})
        self.assertEqual(len(keys), 27)
        for row in inputs:
            matches = {
                (item["source_h5ad"], item["velocity_s_key"])
                for item in gene
                if item["status"] == "ok"
                and item["dataset"] == row["evaluation_dataset"]
                and item["method"] == gene_names[row["config_id"]]
            }
            self.assertEqual(len(matches), 1)
            source, velocity_key = matches.pop()
            self.assertEqual(row["result_path"], source.replace("artifact://velocity/", "velocity:", 1))
            self.assertEqual(row["velocity_layer_key"], velocity_key)

        for table in (primary, supplementary, audit):
            self.assertEqual(len(table), 27)
            self.assertEqual({(r["dataset"], r["method"]) for r in table}, keys)
        source_by_key = {(r["evaluation_dataset"], r["operational_id"]): r["result_path"] for r in inputs}
        for row in audit:
            self.assertEqual(row["result_path"], source_by_key[row["dataset"], row["method"]])
            self.assertEqual(len(row["result_sha256"]), 64)
            self.assertEqual(row["status"], "ok")
            self.assertEqual(int(row["invalid_velocity_values_replaced_with_zero"]), 0)
        for row in primary:
            self.assertTrue(all(np.isfinite(float(row[metric])) for metric in (
                "consistency", "iccoh", "cbdir", "transcosine", "batc"
            )))
        self.assertEqual(len(parity), 54)
        self.assertTrue(all(row["raw_values_match"] == "True" for row in parity))
        self.assertTrue(all(float(row["aggregate_abs_difference"]) <= 1e-7 for row in parity))

        cell = {(r["dataset"], r["method"]): r for r in primary}
        self.assertEqual(len(paired), 36)
        self.assertEqual({r["category"] for r in paired}, {
            "normal + MURK", "normal + branching", "hybrid: MURK", "hybrid: branching"
        })
        for row in paired:
            method = row["operational_id"]
            self.assertAlmostEqual(float(row["cell_cbdir"]), float(cell[row["dataset"], method]["cbdir"]), places=14)
            self.assertAlmostEqual(float(row["cell_batc"]), float(cell[row["dataset"], method]["batc"]), places=14)
            self.assertEqual(int(row["n_genes"]), {
                "normal + MURK": 10, "normal + branching": 20,
                "hybrid: MURK": 20, "hybrid: branching": 20,
            }[row["category"]])
            source = source_by_key[row["dataset"], method].replace(
                "velocity:", "artifact://velocity/", 1
            )
            config_id = next(item["config_id"] for item in inputs
                             if (item["evaluation_dataset"], item["operational_id"])
                             == (row["dataset"], method))
            matching = [item for item in gene if
                        item["dataset"] == row["dataset"]
                        and item["method"] == gene_names[config_id]
                        and item["source_h5ad"] == source]
            if row["category"] == "hybrid: MURK":
                matching = [item for item in matching if item["gene"].startswith((
                    "OneSidedBoost_", "OneSidedBoost2_"
                ))]
            elif row["category"] == "hybrid: branching":
                matching = [item for item in matching if item["gene"].startswith((
                    "MultiForward_", "BoostALL_"
                ))]
            values = np.asarray([float(item["new_cbdir"]) for item in matching])
            self.assertEqual(len(values), int(row["n_genes"]))
            self.assertAlmostEqual(float(row["gene_cbdir_mean"]), float(values.mean()), places=14)
            self.assertAlmostEqual(float(row["gene_cbdir_sd"]), float(values.std(ddof=1)), places=14)
            cos_values = np.asarray([
                cossim_by_key[item["dataset"], item["gene"], item["method"]]
                for item in matching
            ])
            self.assertTrue(np.isfinite(cos_values).all())
            self.assertAlmostEqual(float(row["gene_cossim_mean"]), float(cos_values.mean()), places=14)
            self.assertAlmostEqual(float(row["gene_cossim_sd"]), float(cos_values.std(ddof=1)), places=14)

        for dataset in ("murk", "branching", "all"):
            with np.load(sim / "reference_neighbors" / f"{dataset}.npz", allow_pickle=False) as cell_ref:
                with np.load(root / "evaluation/gene/results/reference_neighbors" / f"{dataset}.npz", allow_pickle=False) as gene_ref:
                    np.testing.assert_array_equal(cell_ref["obs_names"], gene_ref["obs_names"])
                    np.testing.assert_array_equal(cell_ref["indices"], gene_ref["indices"])

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
                resolve_artifact(root, "result.h5ad"), (root / "result.h5ad").resolve()
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
