#!/usr/bin/env python3
"""Recalculate gene-level CBDir with reference-30 neighbors and scoring units.

Method-specific moments and predictions are retained. The extra independent
veloVI min-max transforms are reversed for state and velocity before scoring;
native H5ADs and phase-portrait display coordinates remain untouched.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
import scanpy as sc
import scvelo as scv
from sklearn.metrics.pairwise import cosine_similarity
from scipy import sparse


HERE = Path(__file__).resolve()
REPOSITORY_ROOT = HERE.parents[2]
ARTIFACT_ROOT = None
DEFAULT_OUTPUT = REPOSITORY_ROOT / "evaluation" / "gene" / "results"

# Shared with the panel renderer; the evaluation settings are unchanged.
import sys
sys.path.insert(0, str(REPOSITORY_ROOT))
from evaluation.gene.config import (
    METHOD_ORDER, METHOD_KEYS, SIM_METHOD_KEYS, VELOCITY_KEYS,
    REAL_CONFIGS, SIM_CONFIGS, SIM_PATHS,
)
from evaluation.gene.coordinate_units import restore_velovi_coordinates, velocity_cosine


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_artifact_root(value):
    raw_value = value or os.environ.get("VELOCITY_ARTIFACT_ROOT")
    if not raw_value:
        raise SystemExit(
            "Provide --artifact-root or set VELOCITY_ARTIFACT_ROOT to the "
            "external velocity artifact directory."
        )
    return Path(raw_value).expanduser().resolve()


def public_locator(path):
    """Return a repository-relative or logical external-artifact locator."""
    path = Path(path)
    if not path.is_absolute():
        return path.as_posix()
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPOSITORY_ROOT.resolve()).as_posix()
    except ValueError:
        pass
    if ARTIFACT_ROOT is not None:
        try:
            relative = resolved.relative_to(ARTIFACT_ROOT.resolve())
        except ValueError:
            relative = None
        if relative is not None:
            return f"artifact://velocity/{relative.as_posix()}"
    return f"artifact://external/{resolved.name}"


def dense_values(values):
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values)


def read_reference_asset(path):
    with np.load(path, allow_pickle=False) as data:
        return np.asarray(data["obs_names"]).astype(str), np.asarray(data["indices"], dtype=np.int64)


def copy_reference_asset(dataset, output_path):
    """Reuse the committed graph; ordinary reevaluation never rebuilds it."""
    source = HERE.parent / "results" / "reference_neighbors" / f"{dataset}.npz"
    if not source.is_file():
        raise FileNotFoundError(source)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not output_path.is_file() or sha256(output_path) != sha256(source):
        shutil.copy2(source, output_path)
    names, indices = read_reference_asset(output_path)
    if indices.ndim != 2 or indices.shape != (len(names), 30):
        raise ValueError(f"{dataset}: expected the fixed reference-30 graph")
    return {
        "source": public_locator(source),
        "path": public_locator(output_path),
        "source_sha256": sha256(source),
        "n_obs": len(names),
        "n_neighbors": int(indices.shape[1]),
        "reused": True,
    }


def align_reference_neighbors(adata, reference_names, reference_indices):
    names = np.asarray(adata.obs_names).astype(str)
    if len(names) != len(reference_names) or set(names) != set(reference_names):
        raise ValueError("method cell identifiers do not match reference neighbor asset")
    reference_position = {name: position for position, name in enumerate(reference_names)}
    row_order = np.asarray([reference_position[name] for name in names], dtype=np.int64)
    method_position = np.empty(len(row_order), dtype=np.int64)
    method_position[row_order] = np.arange(len(row_order), dtype=np.int64)
    return method_position[reference_indices[row_order]]


def native_neighbor_jaccard(adata, reference_names, reference_indices):
    if "neighbors" not in adata.uns:
        return None
    native = adata.uns["neighbors"].get("indices")
    if native is None:
        return None
    native = np.asarray(native, dtype=np.int64)
    if native.ndim != 2 or native.shape[0] != adata.n_obs:
        return None
    expected = align_reference_neighbors(adata, reference_names, reference_indices)
    scores = []
    for left, right in zip(native, expected):
        a, b = set(left), set(right)
        scores.append(len(a & b) / len(a | b))
    return float(np.mean(scores))


def native_neighbor_indices(adata, k=30):
    if "neighbors" in adata.uns and "indices" in adata.uns["neighbors"]:
        indices = np.asarray(adata.uns["neighbors"]["indices"], dtype=np.int64)
        if indices.ndim == 2:
            return indices
    if "connectivities" in adata.obsp:
        conn = adata.obsp["connectivities"]
        if sparse.issparse(conn):
            conn = conn.tocsr()
            result = []
            for position in range(conn.shape[0]):
                start, end = conn.indptr[position], conn.indptr[position + 1]
                indices = conn.indices[start:end]
                values = conn.data[start:end]
                keep = indices != position
                indices = indices[keep]
                values = values[keep]
                result.append(indices[np.argsort(values)[::-1]][:k])
            return result
    raise KeyError("native H5AD does not contain neighbors indices or connectivities")


def label_vector(adata, names, reference_label_path=None, preferred_key=None):
    if reference_label_path is not None:
        reference = sc.read_h5ad(reference_label_path, backed="r")
        try:
            labels = np.asarray(reference.obs["clusters"].reindex(names).astype(str))
        finally:
            reference.file.close()
        if np.any(labels == "nan"):
            raise ValueError("reference cluster labels do not cover all method cells")
        return labels, "reference_clusters"
    if preferred_key is not None and preferred_key in adata.obs:
        key = preferred_key
    else:
        key = "clusters" if "clusters" in adata.obs else "bin" if "bin" in adata.obs else None
    if key is None:
        raise KeyError("simulation result lacks clusters/bin labels")
    return np.asarray(adata.obs[key].astype(str)), key


def velocity_matrix(adata, method, gene_indices, unit_multipliers=None):
    velocity_s_key, velocity_u_key = VELOCITY_KEYS[method]
    if velocity_s_key not in adata.layers:
        raise KeyError(f"missing layers/{velocity_s_key}")
    spliced_velocity = dense_values(adata.layers[velocity_s_key][:, gene_indices]).astype(np.float64, copy=False)
    if velocity_u_key not in adata.layers:
        raise KeyError(f"missing layers/{velocity_u_key}")
    unspliced_velocity = dense_values(adata.layers[velocity_u_key][:, gene_indices]).astype(np.float64, copy=False)
    if unit_multipliers is not None:
        spliced_velocity = spliced_velocity * unit_multipliers[0]
        unspliced_velocity = unspliced_velocity * unit_multipliers[1]
    invalid = int(np.size(spliced_velocity) - np.count_nonzero(np.isfinite(spliced_velocity)))
    invalid += int(np.size(unspliced_velocity) - np.count_nonzero(np.isfinite(unspliced_velocity)))
    # This is the old evaluator's in-memory cleanup; source H5AD values remain untouched.
    spliced_velocity = np.clip(np.nan_to_num(spliced_velocity, nan=0.0, posinf=0.0, neginf=0.0), -1e10, 1e10)
    unspliced_velocity = np.clip(np.nan_to_num(unspliced_velocity, nan=0.0, posinf=0.0, neginf=0.0), -1e10, 1e10)
    return spliced_velocity, unspliced_velocity, invalid


def parse_old_score_file(path, method_map):
    records = {}
    if not path.is_file():
        return records
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#") or "ERROR" in line:
                continue
            try:
                parts = line.split("|")
                prefix = parts[0].strip().split(":", 2)
                method = method_map.get(prefix[0].strip())
                if method is None:
                    continue
                if len(prefix) != 3:
                    raise ValueError("expected method: gene: score")
                raw_key, gene, score = prefix
                old_edges = {}
                raw_text = next((part.split("Raw_Scores:", 1)[1].strip() for part in parts if "Raw_Scores:" in part), "{}")
                parsed = ast.literal_eval(raw_text)
                for edge, value in parsed.items():
                    old_edges[f"{edge[0]}->{edge[1]}"] = float(value)
                cossim = None
                for part in parts:
                    if "CosSim:" in part:
                        cossim = float(part.split("CosSim:", 1)[1].strip())
                records[(method, gene.strip())] = {"old_cbdir": float(score.strip()), "old_edges": old_edges, "old_cossim": cossim}
            except (ValueError, SyntaxError, TypeError, AttributeError, IndexError) as error:
                raise ValueError(f"{path}:{line_number}: invalid historical score record") from error
    return records


def score_gene(ms, mu, velocity_s, velocity_u, labels, neighbors, edges):
    phase = np.stack([ms, mu], axis=1)
    velocity = np.stack([velocity_s, velocity_u], axis=1)
    edge_rows = []
    edge_means = []
    evaluable_source_cells = 0
    for source, target in edges:
        source_scores = []
        for cell_index in np.flatnonzero(labels == source):
            target_neighbors = neighbors[cell_index]
            target_neighbors = target_neighbors[labels[target_neighbors] == target]
            if target_neighbors.size == 0:
                continue
            displacement = phase[target_neighbors] - phase[cell_index]
            current_velocity = velocity[cell_index]
            if np.all(current_velocity == 0):
                mean_score = 0.0
            else:
                mean_score = float(np.mean(cosine_similarity(displacement, current_velocity.reshape(1, -1)).ravel()))
            source_scores.append(mean_score)
        if source_scores:
            edge_score = float(np.mean(source_scores))
            edge_means.append(edge_score)
            evaluable_source_cells += len(source_scores)
        else:
            edge_score = float("nan")
        edge_rows.append({"source": source, "target": target, "score": edge_score, "evaluable_source_cells": len(source_scores)})
    global_score = float(np.mean(edge_means)) if edge_means else float("nan")
    return global_score, edge_rows, len(edge_means), evaluable_source_cells


def evaluate_dataset(dataset, config, paths, old_scores, output_root, reference_asset_path, reference_label_path=None, simulation=False, neighbor_mode="reference", methods=METHOD_ORDER):
    reference_names, reference_indices = read_reference_asset(reference_asset_path)
    neighbor_rows = []
    score_rows = []
    edge_rows = []
    unit_rows = []
    cossim_rows = []
    for method in methods:
        path = ARTIFACT_ROOT / paths[method]
        if not path.is_file():
            raise FileNotFoundError(path)
        print(f"[{dataset}] reading {method}: {path}", flush=True)
        adata = sc.read_h5ad(path)
        if neighbor_mode == "native":
            neighbors = native_neighbor_indices(adata)
        else:
            neighbors = align_reference_neighbors(adata, reference_names, reference_indices)
        names = np.asarray(adata.obs_names).astype(str)
        labels, label_source = label_vector(
            adata,
            names,
            reference_label_path=reference_label_path,
            preferred_key="bin" if simulation else None,
        )
        native_jaccard = native_neighbor_jaccard(adata, reference_names, reference_indices)
        neighbor_rows.append(
            {
                "dataset": dataset,
                "method": method,
                "source_h5ad": public_locator(path),
                "source_h5ad_sha256": sha256(path),
                "reference_neighbor_asset": public_locator(reference_asset_path),
                "native_neighbor_jaccard": native_jaccard,
                "label_source": label_source,
                "n_cells": int(adata.n_obs),
                "n_reference_neighbors": int(reference_indices.shape[1]),
            }
        )
        gene_positions = {str(name): position for position, name in enumerate(adata.var_names.astype(str))}
        selected = [gene for gene in config["genes"] if gene in gene_positions]
        if not selected:
            for gene in config["genes"]:
                score_rows.append({"dataset": dataset, "method": method, "gene": gene, "new_cbdir": float("nan"), "old_cbdir": old_scores.get((method, gene), {}).get("old_cbdir", float("nan")), "delta": float("nan"), "status": "gene_missing"})
            continue
        positions = np.asarray([gene_positions[gene] for gene in selected], dtype=np.int64)
        ms_matrix = dense_values(adata.layers["Ms"][:, positions]).astype(np.float64, copy=False) if "Ms" in adata.layers else None
        mu_matrix = dense_values(adata.layers["Mu"][:, positions]).astype(np.float64, copy=False) if "Mu" in adata.layers else None
        if ms_matrix is None or mu_matrix is None:
            raise KeyError(f"{method}/{dataset}: native Ms/Mu are required")
        unit_multipliers = None
        if method == "veloVI":
            moment_source = ARTIFACT_ROOT / "data/capo_pre/capo_pre.h5ad" if dataset == "capo_pre" else path
            preprocessed = sc.read_h5ad(moment_source) if dataset == "capo_pre" else None
            ms_matrix, mu_matrix, inverse_s, inverse_u, factors = restore_velovi_coordinates(
                adata, positions, preprocessed=preprocessed)
            unit_multipliers = (inverse_s, inverse_u)
            for factor in factors:
                unit_rows.append(dict(dataset=dataset, method=method,
                                      moment_source=public_locator(moment_source), **factor))
        velocity_s_matrix, velocity_u_matrix, invalid_velocity = velocity_matrix(
            adata, method, positions, unit_multipliers=unit_multipliers)
        if simulation:
            # CosSim historically uses raw predictions and a norm mask, without
            # CBDir's NaN/Inf replacement or clipping. Truth stays in its units.
            s_key, u_key = VELOCITY_KEYS[method]
            raw_s = dense_values(adata.layers[s_key][:, positions]).astype(np.float64)
            raw_u = dense_values(adata.layers[u_key][:, positions]).astype(np.float64)
            truth_s = dense_values(adata.layers["velocity_s_true"][:, positions]).astype(np.float64)
            truth_u = dense_values(adata.layers["velocity_u_true"][:, positions]).astype(np.float64)
            for column, gene in enumerate(selected):
                native = np.column_stack([raw_s[:, column], raw_u[:, column]])
                restored = (native * np.array([unit_multipliers[0][column], unit_multipliers[1][column]])
                            if unit_multipliers is not None else native)
                truth = np.column_stack([truth_s[:, column], truth_u[:, column]])
                native_cosine, native_n = velocity_cosine(native, truth)
                cosine, n = velocity_cosine(restored, truth)
                cossim_rows.append(dict(dataset=dataset, method=method, gene=gene,
                                        native_cossim=native_cosine, cossim=cosine,
                                        native_evaluable_cells=native_n, evaluable_cells=n,
                                        n_cells=int(adata.n_obs), source_h5ad=public_locator(path)))
        for column, gene in enumerate(selected):
            new_score, gene_edges, n_edges, n_cells = score_gene(ms_matrix[:, column], mu_matrix[:, column], velocity_s_matrix[:, column], velocity_u_matrix[:, column], labels, neighbors, config["edges"])
            old = old_scores.get((method, gene), {})
            old_score = old.get("old_cbdir", float("nan"))
            delta = float(new_score - old_score) if np.isfinite(new_score) and np.isfinite(old_score) else float("nan")
            status = "ok" if np.isfinite(new_score) else "no_evaluable_edges"
            score_rows.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "gene": gene,
                    "new_cbdir": new_score,
                    "old_cbdir": old_score,
                    "delta": delta,
                    "absolute_delta": abs(delta)
                    if np.isfinite(delta)
                    else float("nan"),
                    "n_evaluable_edges": n_edges,
                    "n_evaluable_source_cells": n_cells,
                    "invalid_velocity_values_cleaned_in_memory": invalid_velocity,
                    "status": status,
                    "source_h5ad": public_locator(path),
                    "velocity_s_key": VELOCITY_KEYS[method][0],
                    "velocity_u_key": VELOCITY_KEYS[method][1],
                }
            )
            for edge in gene_edges:
                old_edge = old.get("old_edges", {}).get(f"{edge['source']}->{edge['target']}", float("nan"))
                edge_rows.append({"dataset": dataset, "method": method, "gene": gene, "source": edge["source"], "target": edge["target"], "new_edge_score": edge["score"], "old_edge_score": old_edge, "edge_delta": edge["score"] - old_edge if np.isfinite(edge["score"]) and np.isfinite(old_edge) else float("nan"), "evaluable_source_cells": edge["evaluable_source_cells"]})
        for gene in config["genes"]:
            if gene not in selected:
                old = old_scores.get((method, gene), {})
                score_rows.append({"dataset": dataset, "method": method, "gene": gene, "new_cbdir": float("nan"), "old_cbdir": old.get("old_cbdir", float("nan")), "delta": float("nan"), "absolute_delta": float("nan"), "status": "gene_missing"})
    return score_rows, edge_rows, neighbor_rows, unit_rows, cossim_rows


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fields = list(rows[0])
    for row in rows[1:]:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_change_summaries(output_root, rows):
    grouped = {}
    for row in rows:
        try:
            delta = float(row.get("delta", "nan"))
            old_score = float(row.get("old_cbdir", "nan"))
            new_score = float(row.get("new_cbdir", "nan"))
        except (TypeError, ValueError):
            continue
        if not (np.isfinite(delta) and np.isfinite(old_score) and np.isfinite(new_score)):
            continue
        grouped.setdefault((row["dataset"], row["method"]), []).append((old_score, new_score, delta))

    by_dataset_method = []
    for (dataset, method), values in sorted(grouped.items()):
        old_values, new_values, deltas = map(np.asarray, zip(*values))
        by_dataset_method.append({
            "dataset": dataset,
            "method": method,
            "n_genes": int(len(deltas)),
            "old_mean": float(np.mean(old_values)),
            "new_mean": float(np.mean(new_values)),
            "mean_delta": float(np.mean(deltas)),
            "mean_absolute_delta": float(np.mean(np.abs(deltas))),
            "max_absolute_delta": float(np.max(np.abs(deltas))),
            "n_material_changes_abs_delta_gt_1e-4": int(np.sum(np.abs(deltas) > 1e-4)),
        })
    write_csv(output_root / "gene_cbdir_summary_by_dataset_method.csv", by_dataset_method)

    method_groups = {}
    for row in by_dataset_method:
        method_groups.setdefault(row["method"], []).append(row)
    by_method = []
    for method, entries in sorted(method_groups.items()):
        weights = np.asarray([entry["n_genes"] for entry in entries], dtype=float)
        mean_deltas = np.asarray([entry["mean_delta"] for entry in entries], dtype=float)
        mean_abs = np.asarray([entry["mean_absolute_delta"] for entry in entries], dtype=float)
        by_method.append({
            "method": method,
            "n_dataset_or_scenario": int(len(entries)),
            "n_genes": int(np.sum(weights)),
            "weighted_mean_delta": float(np.average(mean_deltas, weights=weights)),
            "weighted_mean_absolute_delta": float(np.average(mean_abs, weights=weights)),
            "max_absolute_delta": float(max(entry["max_absolute_delta"] for entry in entries)),
            "n_material_changes_abs_delta_gt_1e-4": int(sum(entry["n_material_changes_abs_delta_gt_1e-4"] for entry in entries)),
        })
    write_csv(output_root / "gene_cbdir_summary_by_method.csv", by_method)

    ranked = sorted(
        [row for row in rows if row.get("absolute_delta") not in (None, "", "nan")],
        key=lambda row: float(row["absolute_delta"]),
        reverse=True,
    )
    write_csv(output_root / "gene_cbdir_top_changes.csv", ranked[:100])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=None,
        help="External velocity artifact directory (or VELOCITY_ARTIFACT_ROOT).",
    )
    parser.add_argument("--neighbor-mode", choices=["reference", "native"], default="reference")
    parser.add_argument("--datasets", nargs="+", default=["retina", "DentateGyrus_lamanno", "branching", "murk", "all"])
    parser.add_argument("--methods", nargs="+", choices=METHOD_ORDER, default=METHOD_ORDER)
    args = parser.parse_args()
    global ARTIFACT_ROOT
    ARTIFACT_ROOT = resolve_artifact_root(args.artifact_root)
    output_root = args.output_root if args.output_root.is_absolute() else REPOSITORY_ROOT / args.output_root
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Preserve existing outputs; choose an empty output root: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    neighbor_root = output_root / "reference_neighbors"
    neighbor_root.mkdir(parents=True, exist_ok=True)
    all_scores, all_edges, all_neighbors, all_units, all_cossim = [], [], [], [], []
    provenance = {
        "artifact_root": "artifact://velocity",
        "script": "evaluation/gene/run_reference30_gene_cbdir.py",
        "script_sha256": sha256(HERE),
        "scvelo_version": scv.__version__,
        "scanpy_version": sc.__version__,
        "neighbor_mode": args.neighbor_mode,
        "methods": args.methods,
        "only_neighbor_replaced": False,
        "scoring_units": "method-specific moments with veloVI axis min-max transforms reversed for state and velocity",
        "unit_restoration_script": "evaluation/gene/coordinate_units.py",
        "unit_restoration_script_sha256": sha256(HERE.with_name("coordinate_units.py")),
        "historical_comparison": "current units/reference neighbors versus historical native units/neighbors",
        "native_layers_preserved": [
            "Ms", "Mu", "method_spliced_velocity", "method_unspliced_velocity"
        ],
        "reference_neighbor_spec": {
            "n_neighbors": 30,
            "n_pcs": 30,
            "use_rep": "X_pca",
            "metric": "euclidean",
            "random_state": 0,
        },
        "datasets": {},
    }
    for dataset in args.datasets:
        if dataset in REAL_CONFIGS:
            config = REAL_CONFIGS[dataset]
            asset_path = neighbor_root / dataset
            asset_path = asset_path.with_suffix(".npz")
            reference_info = copy_reference_asset(dataset, asset_path)
            old_path = ARTIFACT_ROOT / "plot" / "final_genes" / dataset / "score.txt"
            old_scores = parse_old_score_file(
                old_path, {value: key for key, value in METHOD_KEYS.items()}
            )
            score_rows, edge_rows, neighbor_rows, unit_rows, cossim_rows = evaluate_dataset(dataset, config, config["paths"], old_scores, output_root, asset_path, reference_label_path=ARTIFACT_ROOT / config["label_reference"], simulation=False, neighbor_mode=args.neighbor_mode, methods=args.methods)
            scope_dir = output_root / "real" / dataset
        elif dataset in SIM_CONFIGS:
            config = SIM_CONFIGS[dataset]
            asset_path = neighbor_root / f"{dataset}.npz"
            reference_info = copy_reference_asset(dataset, asset_path)
            old_path = ARTIFACT_ROOT / "indicators" / "gene" / "simulation" / ("sim_all_analyze" if dataset == "all" else f"{dataset}_analyze") / "results_new" / "score_rerun_results.txt"
            old_scores = parse_old_score_file(
                old_path, {value: key for key, value in SIM_METHOD_KEYS.items()}
            )
            score_rows, edge_rows, neighbor_rows, unit_rows, cossim_rows = evaluate_dataset(dataset, config, SIM_PATHS[dataset], old_scores, output_root, asset_path, simulation=True, neighbor_mode=args.neighbor_mode, methods=args.methods)
            scope_dir = output_root / "simulation" / dataset
        else:
            raise ValueError(f"unknown dataset/scenario: {dataset}")
        write_csv(scope_dir / "gene_cbdir_reference30.csv", score_rows)
        write_csv(scope_dir / "gene_cbdir_edge_details_reference30.csv", edge_rows)
        write_csv(scope_dir / "neighbor_audit.csv", neighbor_rows)
        all_scores.extend(score_rows)
        all_edges.extend(edge_rows)
        all_neighbors.extend(neighbor_rows)
        all_units.extend(unit_rows)
        all_cossim.extend(cossim_rows)
        provenance["datasets"][dataset] = {
            "reference": reference_info,
            "old_score_path": public_locator(old_path),
            "old_score_sha256": sha256(old_path) if old_path.is_file() else None,
            "n_score_rows": len(score_rows),
            "n_edge_rows": len(edge_rows),
            "method_artifact_ids": {
                method: f"gene_h5ad_{dataset}_{method}"
                for method in (REAL_CONFIGS[dataset]["paths"] if dataset in REAL_CONFIGS else SIM_PATHS[dataset])
            },
            "method_paths": (REAL_CONFIGS[dataset]["paths"] if dataset in REAL_CONFIGS else SIM_PATHS[dataset]),
            "label_reference": public_locator(ARTIFACT_ROOT / config["label_reference"]) if dataset in REAL_CONFIGS else None,
            "label_reference_sha256": sha256(ARTIFACT_ROOT / config["label_reference"]) if dataset in REAL_CONFIGS else None,
        }
    write_csv(output_root / "gene_cbdir_reference30_all.csv", all_scores)
    write_csv(output_root / "gene_cbdir_edge_details_reference30_all.csv", all_edges)
    write_csv(output_root / "neighbor_audit_all.csv", all_neighbors)
    write_csv(output_root / "unit_restoration" / "velovi_scale_factors.csv", all_units)
    write_csv(output_root / "simulation" / "gene_cossim.csv", all_cossim)
    write_csv(output_root / "simulation" / "velovi_cossim.csv",
              [row for row in all_cossim if row["method"] == "veloVI"])
    provenance["simulation_cossim"] = {
        "path": "simulation/gene_cossim.csv",
        "methods": args.methods,
        "truth_layers": ["velocity_s_true", "velocity_u_true"],
        "norm_mask": "prediction norm > 1e-6 and truth norm > 1e-6",
        "prediction_cleanup": "none",
        "n_rows": len(all_cossim),
    }
    write_change_summaries(output_root, all_scores)
    provenance["counts"] = {"score_rows": int(len(all_scores)), "edge_rows": int(len(all_edges)), "neighbor_rows": int(len(all_neighbors)), "finite_new_scores": int(sum(np.isfinite(float(row["new_cbdir"])) for row in all_scores if row.get("new_cbdir") not in (None, ""))), "changed_scores": int(sum(np.isfinite(float(row["delta"])) and abs(float(row["delta"])) > 1e-12 for row in all_scores if row.get("delta") not in (None, "")))}
    (output_root / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps(provenance["counts"], indent=2), flush=True)


if __name__ == "__main__":
    main()
