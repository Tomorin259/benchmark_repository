#!/usr/bin/env python3
"""Reproject frozen ensemble operators with GraphVelo's displacement formula."""

from __future__ import annotations

import argparse
import csv
from importlib.metadata import version
import json
from pathlib import Path
import sys

import anndata as ad
from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ensemble.benchmark import resolve_artifact, write_scores
from ensemble.core import (
    DEFAULT_CONFIG, VARIANT_NAMES, VARIANT_SPECS, install_fixed_neighbors,
    load_config, load_fixed_neighbor_asset, project_with_scvelo, sha256,
    validate_signed_operator,
)
from ensemble.metrics import score_global_metrics
from evaluation.cell.reference_embeddings import load_forebrain_frames


ROOT = Path(__file__).resolve().parents[1]
BASELINE_ATOL = 5e-8


def project_graphvelo(operator, coordinates) -> np.ndarray:
    """Apply GraphVelo 0.1.11 project_velocity with an external signed basis.

    Reference: xing-lab-pitt/GraphVelo, graphvelo/graph_velocity.py.
    Each row is sum_j T_ij * (Y_j - Y_i); diagonal weights contribute zero.
    """
    matrix = sparse.csr_matrix(operator, dtype=np.float64, copy=True)
    coordinates = np.asarray(coordinates, dtype=np.float64)
    if coordinates.ndim != 2 or matrix.shape != (len(coordinates), len(coordinates)):
        raise ValueError("operator axes must match the coordinate rows")
    if not np.isfinite(matrix.data).all() or not np.isfinite(coordinates).all():
        raise ValueError("projection inputs must be finite")
    matrix.setdiag(0.0)
    matrix.eliminate_zeros()
    projected = np.zeros(coordinates.shape, dtype=np.float64)
    for cell in range(matrix.shape[0]):
        start, end = matrix.indptr[cell:cell + 2]
        neighbors = matrix.indices[start:end]
        projected[cell] = matrix.data[start:end] @ (
            coordinates[neighbors] - coordinates[cell]
        )
    return projected


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ensemble-root", type=Path, required=True)
    parser.add_argument(
        "--source-manifest", type=Path, required=True,
        help="fixed_neighbors_all14 release manifest with output H5AD hashes",
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--dataset", action="append")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    datasets = args.dataset or list(config["datasets"])
    if len(set(datasets)) != len(datasets) or set(datasets).difference(config["datasets"]):
        raise ValueError("--dataset must contain unique configured dataset keys")
    release = json.loads(args.source_manifest.read_text(encoding="utf-8"))
    baseline_path = ROOT / "ensemble/results/scores.csv"
    with baseline_path.open(newline="", encoding="utf-8") as handle:
        official = {(r["dataset"], r["variant"]): r for r in csv.DictReader(handle)}
    neighbor_root = ROOT / config["reference_neighbors"]
    forebrain = load_forebrain_frames(neighbor_root)
    args.output_root.mkdir(parents=True, exist_ok=False)
    projection_root = args.output_root / "projections"
    projection_root.mkdir()
    metadata = {
        "schema_version": 1,
        "status": "exploratory",
        "datasets": datasets,
        "variants": list(VARIANT_NAMES),
        "source_artifact_id": release["artifact_id"],
        "source_manifest_sha256": sha256(args.source_manifest),
        "config_sha256": sha256(args.config),
        "baseline_scores_sha256": sha256(baseline_path),
        "code_sha256": {
            name: sha256(ROOT / name)
            for name in ("ensemble/reproject.py", "ensemble/core.py",
                         "ensemble/metrics.py", "evaluation/cell/reference_embeddings.py")
        },
        "versions": {"python": sys.version.split()[0], **{
            name: version(name) for name in
            ("scvelo", "anndata", "numpy", "scipy", "h5py", "pandas", "scikit-learn")
        }},
        "projection": {
            "reference": "GraphVelo 0.1.11 GraphVelo.project_velocity",
            "source_url": "https://github.com/xing-lab-pitt/GraphVelo/blob/main/graphvelo/graph_velocity.py",
            "external_basis": "off-diagonal of the frozen ensemble_transition_signed",
            "formula": "sum_j T_ij * (Y_j - Y_i)",
            "normalization": "retain signed weights and native directional mass",
            "training": "none",
            "randomness": "none; stored operators, coordinates and neighbors",
        },
        "baseline_atol": BASELINE_ATOL,
        "forebrain_frames": forebrain["metadata"],
        "inputs": [],
        "metric_details": [],
    }
    scores, comparisons = [], []
    for dataset in datasets:
        spec = config["datasets"][dataset]
        asset_path = neighbor_root / spec["reference_neighbor_file"]
        asset = load_fixed_neighbor_asset(asset_path)
        bases = ("pca", "umap") if dataset == "forebrain" else ("umap",)
        arrays = {}
        for variant in VARIANT_NAMES:
            entry = release["outputs"][f"output:{variant}:{dataset}:h5ad"]
            path = resolve_artifact(args.ensemble_root, entry["relative_path"])
            before = path.stat()
            if before.st_size != entry["size_bytes"] or sha256(path) != entry["sha256"]:
                raise ValueError(f"source differs from release manifest: {entry['relative_path']}")
            with h5py.File(path, "r") as handle:
                params = handle["uns/ensemble_signed_params"]
                schema = int(params["schema_version"][()])
                if schema not in (6, 7) or params["aggregation_variant"].asstr()[()] != variant:
                    raise ValueError(f"unexpected ensemble metadata: {entry['relative_path']}")
                obs = read_elem(handle["obs"])[[spec["label"]]]
                names = obs.index.to_numpy(dtype=str)
                operator = read_elem(handle["obsp/ensemble_transition_signed"])
                stored_neighbors = handle["uns/neighbors/indices"][:]
                coordinates = {"umap": handle["obsm/X_umap"][:]}
                baseline = {"umap": handle["obsm/ensemble_velocity_umap"][:]}
                if dataset == "forebrain":
                    positions = pd.Index(forebrain["obs_names"]).get_indexer(names)
                    if np.any(positions < 0) or len(names) != len(forebrain["obs_names"]):
                        raise ValueError("Forebrain cells differ from the frozen frames")
                    coordinates = {basis: forebrain[f"X_{basis}"][positions] for basis in bases}
                    if not np.array_equal(handle["obsm/X_pca"][:], coordinates["pca"]):
                        raise ValueError("Forebrain PCA differs from the frozen frame")
                    if schema == 6:
                        if not np.array_equal(handle["obsm/X_umap"][:], coordinates["pca"]):
                            raise ValueError("schema-6 Forebrain field is not the accepted PCA30 alias")
                        baseline["pca"] = baseline.pop("umap")
                        moments = read_elem(handle["layers/Ms"])
                        reference = ad.AnnData(obs=obs.copy(), shape=moments.shape)
                        reference.layers["Ms"] = moments
                        for basis in bases:
                            reference.obsm[f"X_{basis}"] = coordinates[basis]
                        replay, invalid = project_with_scvelo(
                            reference, operator, variant=variant,
                            logical_inputs={}, bases=("umap",),
                        )
                        if invalid:
                            raise ValueError("Forebrain baseline replay produced invalid rows")
                        baseline["umap"] = replay.obsm["ensemble_velocity_umap"].copy()
                    else:
                        if not np.array_equal(handle["obsm/X_umap"][:], coordinates["umap"]):
                            raise ValueError("Forebrain UMAP differs from the frozen frame")
                        baseline["pca"] = handle["obsm/ensemble_velocity_pca"][:]

            validate_signed_operator(operator, entry["relative_path"])
            prepared = ad.AnnData(obs=obs.copy())
            install_fixed_neighbors(prepared, asset)
            if not np.array_equal(prepared.uns["neighbors"]["indices"], stored_neighbors):
                raise ValueError(f"stored reference neighbors differ: {dataset}/{variant}")
            for basis in bases:
                prepared.obsm[f"X_{basis}"] = coordinates[basis]
                prepared.obsm[f"ensemble_velocity_{basis}"] = baseline[basis]
            old = score_global_metrics(prepared, dataset, spec)
            frozen = official[dataset, variant]
            for metric in ("cbdir", "batc"):
                if not np.isclose(old[metric], float(frozen[metric]), rtol=0, atol=BASELINE_ATOL):
                    raise ValueError(f"baseline mismatch: {dataset}/{variant}/{metric}: "
                                     f"{old[metric]} != {frozen[metric]}")
            if not arrays:
                arrays = {
                    "cell_ids": names, "labels": obs[spec["label"]].to_numpy(dtype=str),
                    "neighbors": stored_neighbors,
                    **{f"X_{basis}": coordinates[basis] for basis in bases},
                }
            elif (not np.array_equal(arrays["cell_ids"], names)
                  or not np.array_equal(arrays["labels"], obs[spec["label"]].to_numpy(dtype=str))
                  or any(not np.array_equal(arrays[f"X_{b}"], coordinates[b]) for b in bases)):
                raise ValueError(f"variant cells, labels or coordinates differ: {dataset}")
            for basis in bases:
                projected = project_graphvelo(operator, coordinates[basis])
                prepared.obsm[f"ensemble_velocity_{basis}"] = projected
                arrays[f"{variant}_{basis}"] = projected
            new = score_global_metrics(prepared, dataset, spec)
            for key in ("cbdir_edge_evaluable_source_cells", "batc_evaluable_source_cells"):
                if old[key] != new[key]:
                    raise ValueError(f"metric coverage changed: {dataset}/{variant}/{key}")
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError(f"source changed during evaluation: {entry['relative_path']}")
            scores.append({
                "dataset": dataset, "variant": variant,
                "label": VARIANT_SPECS[variant]["label"],
                "cbdir": new["cbdir"], "batc": new["batc"],
            })
            comparison = {"dataset": dataset, "variant": variant, "n_obs": len(obs)}
            for metric in ("cbdir", "batc"):
                comparison.update({
                    f"scvelo_{metric}": old[metric],
                    f"graphvelo_projection_{metric}": new[metric],
                    f"delta_{metric}": new[metric] - old[metric],
                    f"baseline_replay_error_{metric}": old[metric] - float(frozen[metric]),
                })
            comparisons.append(comparison)
            metadata["inputs"].append({
                "dataset": dataset, "variant": variant, **entry,
                "schema_version": schema, "neighbor_sha256": asset["sha256"],
            })
            metadata["metric_details"].append({
                "dataset": dataset, "variant": variant,
                "cbdir_basis": spec.get("cbdir_basis", "umap"),
                "scvelo": old, "graphvelo_projection": new,
            })
            print(f"{dataset}/{variant}: "
                  f"CBDir {old['cbdir']:.6f} -> {new['cbdir']:.6f}; "
                  f"BATC {old['batc']:.6f} -> {new['batc']:.6f}", flush=True)
        np.savez_compressed(projection_root / f"{dataset}.npz", **arrays)

    write_scores(args.output_root / "scores.csv", scores)
    with (args.output_root / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparisons)
    summary = []
    for variant in VARIANT_NAMES:
        selected = [row for row in comparisons if row["variant"] == variant]
        for metric in ("cbdir", "batc"):
            old = float(np.mean([row[f"scvelo_{metric}"] for row in selected]))
            new = float(np.mean([row[f"graphvelo_projection_{metric}"] for row in selected]))
            delta = np.asarray([row[f"delta_{metric}"] for row in selected])
            summary.append({
                "variant": variant, "metric": metric, "n_datasets": len(selected),
                "scvelo_mean": old, "graphvelo_projection_mean": new, "delta_mean": new - old,
                "improved": int(np.count_nonzero(delta > 1e-12)),
                "reduced": int(np.count_nonzero(delta < -1e-12)),
                "tied": int(np.count_nonzero(np.abs(delta) <= 1e-12)),
            })
    with (args.output_root / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    metadata["outputs_sha256"] = {
        str(path.relative_to(args.output_root)): sha256(path)
        for path in sorted(args.output_root.rglob("*")) if path.is_file()
    }
    (args.output_root / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
