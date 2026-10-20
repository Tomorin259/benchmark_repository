#!/usr/bin/env python3
"""Fixed-neighbor CBDir and BATC for the signed-T ensemble outputs."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.spatial import cKDTree
from sklearn.metrics.pairwise import cosine_similarity


BATC_USE_BINS = True
BATC_MIN_PER_BIN = 5
BATC_N_SAMPLES = 800
BATC_DERIVATIVE_STEP = 1e-3


def _as_float_array(value, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} contains NaN or infinity")
    return result


def _validate_inputs(adata, dataset: str, dataset_config: dict) -> tuple:
    label_key = dataset_config["label"]
    if label_key not in adata.obs:
        raise KeyError(f"{dataset}: missing obs[{label_key!r}]")
    if "X_umap" not in adata.obsm:
        raise KeyError(f"{dataset}: missing obsm['X_umap']")
    if "ensemble_velocity_umap" not in adata.obsm:
        raise KeyError(f"{dataset}: missing obsm['ensemble_velocity_umap']")
    neighbors = adata.uns.get("neighbors", {}).get("indices")
    if neighbors is None:
        raise KeyError(f"{dataset}: missing fixed neighbor indices")
    neighbors = np.asarray(neighbors, dtype=np.int64)
    if neighbors.ndim != 2 or neighbors.shape[0] != adata.n_obs:
        raise ValueError(f"{dataset}: invalid fixed neighbor indices shape")
    if neighbors.shape[1] < 2:
        raise ValueError(f"{dataset}: fixed neighbors need self plus one neighbor")
    if np.any(neighbors < 0) or np.any(neighbors >= adata.n_obs):
        raise ValueError(f"{dataset}: fixed neighbors contain an invalid index")
    if not np.array_equal(neighbors[:, 0], np.arange(adata.n_obs)):
        raise ValueError(f"{dataset}: fixed neighbors must store self first")
    labels = np.asarray(adata.obs[label_key].astype(str))
    required = {str(label) for edge in dataset_config["edges"] for label in edge}
    missing = sorted(required.difference(set(labels)))
    if missing:
        raise ValueError(f"{dataset}: missing configured labels {missing}")
    return labels, neighbors


def compute_cbdir(
    labels: np.ndarray,
    neighbors: np.ndarray,
    embedding: np.ndarray,
    velocity: np.ndarray,
    edges,
) -> tuple[float, dict]:
    """Compute equal-edge CBDir using the audited fixed neighbor rows."""
    embedding = _as_float_array(embedding, "CBDir embedding")
    velocity = _as_float_array(velocity, "CBDir velocity")
    if embedding.ndim != 2 or velocity.shape != embedding.shape:
        raise ValueError("CBDir embedding and velocity shapes do not match")
    edge_scores = {}
    edge_counts = {}
    valid_edge_values = []
    for source, target in edges:
        source = str(source)
        target = str(target)
        source_ids = np.flatnonzero(labels == source)
        cell_scores = []
        for cell_id in source_ids:
            target_ids = neighbors[cell_id][labels[neighbors[cell_id]] == target]
            if target_ids.size == 0:
                continue
            displacement = embedding[target_ids] - embedding[cell_id]
            scores = cosine_similarity(
                displacement, velocity[cell_id].reshape(1, -1)
            ).ravel()
            cell_scores.append(float(np.mean(scores)))
        key = f"{source} -> {target}"
        edge_counts[key] = len(cell_scores)
        if cell_scores:
            value = float(np.mean(cell_scores))
            edge_scores[key] = value
            valid_edge_values.append(value)
    score = float(np.mean(valid_edge_values)) if valid_edge_values else float("nan")
    return score, {
        "edge_scores": edge_scores,
        "edge_evaluable_source_cells": edge_counts,
        "evaluable_edges": len(valid_edge_values),
    }


class _ParametricCurve2D:
    def __init__(self, parameters: np.ndarray, points: np.ndarray):
        self._x = PchipInterpolator(parameters, points[:, 0])
        self._y = PchipInterpolator(parameters, points[:, 1])

    def evaluate(self, values) -> np.ndarray:
        values = np.clip(np.asarray(values, dtype=np.float64), 0.0, 1.0)
        return np.column_stack((self._x(values), self._y(values)))

    def derivative(self, values) -> np.ndarray:
        values = np.asarray(values, dtype=np.float64)
        lower = np.clip(values - BATC_DERIVATIVE_STEP, 0.0, 1.0)
        upper = np.clip(values + BATC_DERIVATIVE_STEP, 0.0, 1.0)
        denominator = (upper - lower)[:, None]
        denominator[denominator == 0.0] = 1e-12
        return (self.evaluate(upper) - self.evaluate(lower)) / denominator


def _remove_consecutive_duplicates(points: np.ndarray) -> np.ndarray:
    if len(points) < 2:
        return points
    keep = np.r_[True, np.linalg.norm(np.diff(points, axis=0), axis=1) > 0]
    return points[keep]


def _fit_batc_curve(
    points: np.ndarray,
    source_center: np.ndarray,
    target_center: np.ndarray,
    *,
    n_bins: int,
    min_per_bin: int,
    n_samples: int,
) -> tuple[_ParametricCurve2D, np.ndarray, np.ndarray, float]:
    direction = target_center - source_center
    direction_norm = np.linalg.norm(direction)
    unit_direction = direction / (direction_norm + 1e-12)
    projection = (points - source_center) @ unit_direction
    order = np.argsort(projection, kind="stable")
    ordered_points = points[order]
    ordered_projection = projection[order]

    if BATC_USE_BINS and len(ordered_points) >= n_bins:
        edges = np.linspace(
            ordered_projection.min(), ordered_projection.max(), n_bins + 1
        )
        anchors = []
        for index in range(n_bins):
            if index < n_bins - 1:
                selected = (ordered_projection >= edges[index]) & (
                    ordered_projection < edges[index + 1]
                )
            else:
                selected = (ordered_projection >= edges[index]) & (
                    ordered_projection <= edges[index + 1]
                )
            if np.count_nonzero(selected) >= min_per_bin:
                anchors.append(ordered_points[selected].mean(axis=0))
        skeleton = np.vstack(anchors) if len(anchors) >= 3 else ordered_points
    else:
        skeleton = ordered_points

    skeleton = _remove_consecutive_duplicates(skeleton)
    if len(skeleton) < 2:
        skeleton = np.vstack((ordered_points[0], ordered_points[-1]))
    lengths = np.linalg.norm(np.diff(skeleton, axis=0), axis=1)
    parameters = np.r_[0.0, np.cumsum(lengths)]
    if parameters[-1] == 0.0:
        parameters[-1] = 1.0
    parameters /= parameters[-1]
    curve = _ParametricCurve2D(parameters, skeleton)

    grid = np.linspace(0.0, 1.0, int(n_samples))
    _, nearest = cKDTree(curve.evaluate(grid)).query(points, k=1)
    point_parameters = grid[nearest]
    tangent = curve.evaluate([1.0])[0] - curve.evaluate([0.0])[0]
    orientation = (
        -1.0
        if direction_norm > 0 and np.dot(tangent, direction) < 0
        else 1.0
    )
    return curve, point_parameters[order], order, orientation


def compute_batc(
    labels: np.ndarray,
    embedding: np.ndarray,
    velocity: np.ndarray,
    edges,
    *,
    n_bins: int,
    min_per_bin: int = BATC_MIN_PER_BIN,
    n_samples: int = BATC_N_SAMPLES,
) -> tuple[float, dict]:
    """Compute the branch-aware tangent consistency used by BATC."""
    embedding = _as_float_array(embedding, "BATC embedding")
    velocity = _as_float_array(velocity, "BATC velocity")
    if embedding.ndim != 2 or embedding.shape[1] != 2:
        raise ValueError("BATC expects a two-dimensional embedding")
    if velocity.shape != embedding.shape:
        raise ValueError("BATC embedding and velocity shapes do not match")
    if n_bins < 1 or min_per_bin < 1 or n_samples < 2:
        raise ValueError("BATC parameters must be positive")

    unique_labels = np.unique(labels)
    indices = {
        label: np.flatnonzero(labels == label) for label in unique_labels
    }
    centers = {
        label: embedding[cell_ids].mean(axis=0)
        for label, cell_ids in indices.items()
        if cell_ids.size
    }
    edge_cosines = {}
    edge_list = [(str(a), str(b)) for a, b in edges]
    for source, target in edge_list:
        source_ids = indices[source]
        target_ids = indices[target]
        combined = np.concatenate((source_ids, target_ids))
        curve, parameters, order, orientation = _fit_batc_curve(
            embedding[combined],
            centers[source],
            centers[target],
            n_bins=n_bins,
            min_per_bin=min_per_bin,
            n_samples=n_samples,
        )
        ordered_ids = combined[order]
        mapping = {}
        for position, cell_id in enumerate(ordered_ids):
            tangent = orientation * curve.derivative([parameters[position]])[0]
            cell_velocity = velocity[cell_id]
            denominator = np.linalg.norm(tangent) * np.linalg.norm(cell_velocity)
            mapping[int(cell_id)] = (
                float(np.dot(tangent, cell_velocity) / denominator)
                if denominator
                else np.nan
            )
        edge_cosines[(source, target)] = mapping

    targets_by_source = {}
    for source, target in edge_list:
        targets_by_source.setdefault(source, []).append(target)
    best = np.full(embedding.shape[0], np.nan, dtype=np.float64)
    for source, targets in targets_by_source.items():
        source_ids = indices.get(source, np.empty(0, dtype=np.int64))
        if source_ids.size == 0:
            continue
        scores = np.asarray(
            [
                [
                    edge_cosines.get((source, target), {}).get(int(cell_id), np.nan)
                    for cell_id in source_ids
                ]
                for target in targets
            ],
            dtype=np.float64,
        )
        best[source_ids] = np.nan_to_num(scores, nan=0.0).max(axis=0)
    valid = np.isfinite(best)
    score = float(np.nanmean(best[valid])) if valid.any() else float("nan")
    return score, {
        "evaluable_source_cells": int(valid.sum()),
        "configured_edges": len(edge_list),
    }


def score_global_metrics(
    adata,
    dataset: str,
    dataset_config: dict,
) -> dict:
    """Return CBDir, BATC, and compact evaluator diagnostics."""
    labels, neighbors = _validate_inputs(adata, dataset, dataset_config)
    cbdir_basis = dataset_config.get("cbdir_basis", "umap")
    embedding = _as_float_array(adata.obsm[f"X_{cbdir_basis}"], "CBDir embedding")
    velocity = _as_float_array(
        adata.obsm[f"ensemble_velocity_{cbdir_basis}"], "embedded CBDir velocity"
    )
    if embedding.ndim != 2 or velocity.shape != embedding.shape:
        raise ValueError(f"{dataset}: embedding and velocity shapes do not match")
    expected_dimensions = dataset_config.get("cbdir_dimensions", 2)
    if embedding.shape[1] != expected_dimensions:
        raise ValueError(f"{dataset}: CBDir requires {expected_dimensions} {cbdir_basis} dimensions")
    cbdir, cbdir_details = compute_cbdir(
        labels, neighbors, embedding, velocity, dataset_config["edges"]
    )
    batc_embedding = _as_float_array(adata.obsm["X_umap"], "BATC UMAP embedding")
    batc_velocity = _as_float_array(
        adata.obsm["ensemble_velocity_umap"], "embedded BATC velocity"
    )
    batc, batc_details = compute_batc(
        labels,
        batc_embedding,
        batc_velocity,
        dataset_config["edges"],
        n_bins=int(dataset_config["batc_n_bins"]),
    )
    return {
        "cbdir": cbdir,
        "batc": batc,
        "cbdir_edge_scores": cbdir_details["edge_scores"],
        "cbdir_edge_evaluable_source_cells": cbdir_details[
            "edge_evaluable_source_cells"
        ],
        "cbdir_evaluable_edges": cbdir_details["evaluable_edges"],
        "batc_evaluable_source_cells": batc_details["evaluable_source_cells"],
        "batc_configured_edges": batc_details["configured_edges"],
    }


RANKING_FIELDS = ["dataset", "variant", "label", "metric", "score", "rank", "out_of", "tie_count", "methods_above"]

def _read_csv(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))

def build_rankings(scores, cell_metrics, dataset_map):
    """Purely compute independent 21-way competition ranks from row mappings."""
    methods = sorted({row["method"] for row in cell_metrics})
    if len(methods) != 20:
        raise ValueError("Cell metrics must contain exactly 20 methods")
    score_by_key = {(row["dataset"], row["variant"]): row for row in scores}
    expected = {(d, v) for d in dataset_map for v in ("signed_t_mean4", "signed_t_row_geometric_median4")}
    if len(score_by_key) != len(scores) or set(score_by_key) != expected:
        raise ValueError("scores must contain exactly two unique variants per dataset")
    cell_by_key = {(row["dataset"], row["method"]): row for row in cell_metrics}
    expected_cell = {(dataset_map[d], m) for d in dataset_map for m in methods}
    if len(cell_by_key) != len(cell_metrics) or set(cell_by_key) != expected_cell:
        raise ValueError("Cell metrics must contain a complete mapped 14x20 matrix")
    labels = {
        "signed_t_mean4": "Signed-T Mean-4",
        "signed_t_row_geometric_median4": "Signed-T RowGeometricMedian-4",
    }
    rows = []
    for dataset in dataset_map:
        for variant, label in labels.items():
            source = score_by_key[(dataset, variant)]
            for metric in ("cbdir", "batc"):
                score = float(source[metric])
                values = [float(cell_by_key[(dataset_map[dataset], method)][metric]) for method in methods]
                above = sum(value > score for value in values)
                ties = sum(bool(np.isclose(value, score, rtol=1e-12, atol=1e-12)) for value in values)
                rows.append({"dataset": dataset, "variant": variant, "label": label, "metric": metric, "score": source[metric], "rank": 1 + above, "out_of": 21, "tie_count": ties, "methods_above": above})
    return rows


def write_rankings(scores_path: Path, cell_metrics_path: Path, config_path: Path, rankings_path: Path, summary_path: Path) -> None:
    import json
    config = json.loads(Path(config_path).read_text(encoding="utf-8"))
    scores, cell = _read_csv(scores_path), _read_csv(cell_metrics_path)
    dataset_map = {"bonemarrow": "bonemarrow", "capofibro": "capofibro", "dentategyrus": "dentategyrus", "dentate_lamanno": "DentateGyrus_lamanno", "forebrain": "forebrain", "gastrulation": "gastrulation", "hematopoiesis": "hematopoiesis", "mousebrain": "mousebrain", "neuron": "neuron", "organoid": "organoid", "pancreas": "pancreas", "retina": "retina", "rpe1": "rpe1", "rpe1_gen": "rpe1_gen"}
    if set(dataset_map) != set(config["datasets"]):
        raise ValueError("ranking dataset map is out of sync with ensemble config")
    cell_datasets = {row["dataset"] for row in cell}
    print(f"dataset mapping validation: config={sorted(dataset_map)}, cell={sorted(cell_datasets)}, mapped={sorted(set(dataset_map.values()))}")
    if set(dataset_map.values()) != cell_datasets or len(dataset_map) != len(set(dataset_map.values())):
        raise ValueError("config and Cell dataset sets do not map one-to-one")
    rows = build_rankings(scores, cell, dataset_map)
    with Path(rankings_path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RANKING_FIELDS, lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    lines = ["# Ensemble ranking versus 20 methods", "", "Each score is ranked independently within its dataset and metric against the 20 Cell methods.", "The ensemble is the 21st participant; equal scores use competition rank.", ""]
    for variant, label in (
        ("signed_t_mean4", "Signed-T Mean-4"),
        ("signed_t_row_geometric_median4", "Signed-T RowGeometricMedian-4"),
    ):
        selected = [r for r in rows if r["variant"] == variant]
        average = sum(r["rank"] for r in selected) / len(selected)
        lines += [f"## {label}", "", f"First: {sum(r['rank'] == 1 for r in selected)} | Top 3: {sum(r['rank'] <= 3 for r in selected)} | Top 5: {sum(r['rank'] <= 5 for r in selected)} | Average rank: {average:.6f} (28 dataset-metric records)", "", "| Dataset | CBDir rank | BATC rank |", "| --- | ---: | ---: |"]
        for dataset in dataset_map:
            get = lambda metric: next(r["rank"] for r in selected if r["dataset"] == dataset and r["metric"] == metric)
            lines.append(f"| {dataset} | {get('cbdir')}/21 | {get('batc')}/21 |")
        lines.append("")
    Path(summary_path).write_text("\n".join(lines), encoding="utf-8")
