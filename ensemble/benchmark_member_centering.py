#!/usr/bin/env python3
"""Compare member-before-median centering with frozen projection controls."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import sparse

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ensemble.core import (
    DEFAULT_CONFIG, METHOD_NAMES, align_method_values, load_config,
    offdiagonal, row_l1_norm, sha256, validate_scvelo_runtime,
)
from ensemble.geometric_median import row_geometric_median
from ensemble.pareto_members import read_source, score_fields
from ensemble.reproject import project_graphvelo

ROOT = Path(__file__).resolve().parents[1]
CONDITIONS = ("current", "raw_centered", "member_centered")
METRICS = ("cbdir", "batc")
SCORE_ATOL = 5e-7
FIELD_ATOL = 1e-9
ROW_SUM_ATOL = 1e-10


def center_neighbor_weights(operator):
    """Center each row on its own nonzero, off-diagonal support; no rescaling."""
    centered = offdiagonal(operator)
    centered.sum_duplicates()
    centered.eliminate_zeros()
    counts = np.diff(centered.indptr)
    total = np.asarray(centered.sum(axis=1)).ravel()
    means = np.divide(total, counts, out=np.zeros_like(total), where=counts > 0)
    centered.data -= np.repeat(means, counts)
    centered.eliminate_zeros()
    return centered


def aggregate_rows(operators):
    """Reuse the accepted Euclidean row median without diagonal completion."""
    template, values = align_method_values(operators)
    aggregate, diagnostics, cells = row_geometric_median(
        template, values, return_cell_diagnostics=True
    )
    result = sparse.csr_matrix(
        (aggregate, template.indices.copy(), template.indptr.copy()),
        shape=template.shape,
    )
    result.eliminate_zeros()
    return result, diagnostics, cells["coefficients"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--operator-root", type=Path, required=True)
    parser.add_argument("--ablation-root", type=Path, required=True)
    parser.add_argument("--selected-root", type=Path,
                        default=ROOT / "ensemble/results/legacy_geomedian4_graphvelo",
                        help="frozen four-member GraphVelo-projected baseline bundle")
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args(argv)
    validate_scvelo_runtime()
    config = load_config()
    source_path = args.operator_root / "metadata.json"
    source = json.loads(source_path.read_text())
    if tuple(source["protocol"]["methods"][:4]) != METHOD_NAMES:
        raise ValueError("cached original member identities or order differ")
    if source["completed_datasets"] != list(config["datasets"]):
        raise ValueError("cached datasets differ from the frozen fourteen")
    if source["protocol"]["config_sha256"] != sha256(DEFAULT_CONFIG):
        raise ValueError("cached configuration changed")
    for name, expected in source["protocol"]["code_sha256"].items():
        if sha256(ROOT / name) != expected:
            raise ValueError("cached computation code changed: " + name)

    selected_meta_path = args.selected_root / "metadata.json"
    selected_meta = json.loads(selected_meta_path.read_text())
    selected_scores_path = args.selected_root / "scores.csv"
    if sha256(selected_scores_path) != selected_meta["outputs_sha256"]["scores.csv"]:
        raise ValueError("selected baseline table changed")
    selected_scores = pd.read_csv(selected_scores_path).set_index("dataset")
    ablation_meta_path = args.ablation_root / "metadata.json"
    ablation_meta = json.loads(ablation_meta_path.read_text())
    ablation_scores_path = args.ablation_root / "scores.csv"
    if sha256(ablation_scores_path) != ablation_meta["scores_sha256"]:
        raise ValueError("raw-centered baseline table changed")
    ablation_scores = pd.read_csv(ablation_scores_path)
    ablation_scores = ablation_scores[
        (ablation_scores.variant == "signed_t_row_geometric_median4")
        & (ablation_scores.projection == "raw_centered")
    ].set_index("dataset")
    if set(ablation_scores.index) != set(config["datasets"]):
        raise ValueError("raw-centered baseline coverage differs")

    code_names = [
        "ensemble/benchmark_member_centering.py", "ensemble/core.py",
        "ensemble/geometric_median.py", "ensemble/pareto_members.py",
        "ensemble/metrics.py", "ensemble/reproject.py",
    ]
    protocol = {
        "status": "exploratory",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "members": list(METHOD_NAMES),
        "conditions": {
            "current": "GM(T1,T2,T3,T4), then raw displacement projection",
            "raw_centered": "center GM(T1,T2,T3,T4), then raw displacement projection",
            "member_centered": "GM(center(T1),...,center(T4)), then raw displacement projection",
        },
        "centering": {
            "axis": "neighbor-edge weights of each source cell",
            "support": "each operator's own nonzero off-diagonal entries before zero-padding",
            "empty_or_single_edge_rows": "zero",
            "post_centering_normalization": "none",
            "post_aggregation_normalization": "none",
            "diagonal_completion": "none; diagonal displacement is zero",
        },
        "projection": "sum_j W_ij (Y_j-Y_i); raw displacements; no additional centering",
        "evaluation": "unchanged fixed CBDir/BATC, reference neighbors and coordinate frames",
        "randomness": "none; reuse frozen member operators and coordinates",
        "operator_source_metadata_sha256": sha256(source_path),
        "selected_metadata_sha256": sha256(selected_meta_path),
        "selected_scores_sha256": sha256(selected_scores_path),
        "ablation_metadata_sha256": sha256(ablation_meta_path),
        "ablation_scores_sha256": sha256(ablation_scores_path),
        "config_sha256": sha256(DEFAULT_CONFIG),
        "code_sha256": {name: sha256(ROOT / name) for name in code_names},
        "repository_head": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True
        ).strip(),
        "versions": {name: version(name) for name in
                     ("numpy", "scipy", "pandas", "anndata", "scvelo", "scikit-learn")},
        "score_atol": SCORE_ATOL, "field_atol": FIELD_ATOL,
    }
    output = args.output_root
    output.mkdir(parents=True, exist_ok=False)
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")

    rows, records = [], {}
    for dataset, spec in config["datasets"].items():
        print(dataset + ": validating frozen inputs", flush=True)
        plot, arrays, operators, source_detail = read_source(
            args.operator_root, source, dataset, spec
        )
        operators = operators[:4]
        bases = ("pca", "umap") if dataset == "forebrain" else ("umap",)
        selected_path = args.selected_root / "projections" / (dataset + ".npz")
        if sha256(selected_path) != selected_meta["outputs_sha256"]["projections/" + dataset + ".npz"]:
            raise ValueError("selected projection archive changed: " + dataset)
        with np.load(selected_path, allow_pickle=False) as selected:
            for key in ("cell_ids", "labels", "neighbors", *["X_" + b for b in bases]):
                if not np.array_equal(arrays[key].astype(str) if key == "labels" else arrays[key],
                                      selected[key].astype(str) if key == "labels" else selected[key]):
                    raise ValueError("frozen selected axes or geometry differ: " + dataset + "/" + key)

        print(dataset + ": aggregating original and centered members", flush=True)
        current, current_diagnostics, _ = aggregate_rows(operators)
        centered_members = [center_neighbor_weights(matrix) for matrix in operators]
        member_centered, new_diagnostics, coefficients = aggregate_rows(centered_members)
        conditions = {
            "current": current,
            "raw_centered": center_neighbor_weights(current),
            "member_centered": member_centered,
        }
        zero_sum_error = max(
            float(np.abs(np.asarray(matrix.sum(axis=1)).ravel()).max(initial=0))
            for matrix in [*centered_members, member_centered, conditions["raw_centered"]]
        )
        if zero_sum_error > ROW_SUM_ATOL:
            raise ValueError("centered rows do not have zero signed sum")

        saved = {
            "cell_ids": arrays["cell_ids"].astype(str),
            "labels": arrays["labels"].astype(str),
            "neighbors": arrays["neighbors"],
            "member_centered_coefficients": coefficients,
            **{"X_" + b: arrays["X_" + b] for b in bases},
        }
        condition_details = {}
        for condition, operator in conditions.items():
            fields = {b: project_graphvelo(operator, arrays["X_" + b]) for b in bases}
            if not all(np.isfinite(field).all() for field in fields.values()):
                raise ValueError("nonfinite projected field")
            field_error = 0.0
            if condition == "current":
                field_error = max(
                    float(np.abs(fields[b] - arrays["geomedian4_" + b]).max(initial=0))
                    for b in bases
                )
                if field_error > FIELD_ATOL:
                    raise ValueError("current projection control changed")
            score = score_fields(plot, fields, dataset, spec)
            score_error = 0.0
            if condition != "member_centered":
                expected = selected_scores if condition == "current" else ablation_scores
                score_error = max(abs(score[m] - float(expected.loc[dataset, m])) for m in METRICS)
                if score_error > SCORE_ATOL:
                    raise ValueError("control scores changed: " + dataset + "/" + condition)
            for basis, field in fields.items():
                saved[condition + "_" + basis] = field
            rows.append({
                "dataset": dataset, "condition": condition,
                **{m: score[m] for m in METRICS},
                "control_score_max_abs_error": score_error,
            })
            condition_details[condition] = {
                "metrics": score, "field_control_max_abs_error": field_error,
                "score_control_max_abs_error": score_error,
                "operator_max_abs_row_sum": float(np.abs(np.asarray(operator.sum(axis=1))).max(initial=0)),
                "offdiagonal_l1_mean": float(row_l1_norm(operator).mean()),
            }
            print(json.dumps(rows[-1]), flush=True)
        coverage = condition_details["current"]["metrics"]
        for condition in CONDITIONS[1:]:
            candidate = condition_details[condition]["metrics"]
            for key in ("cbdir_edge_evaluable_source_cells", "batc_evaluable_source_cells"):
                if candidate[key] != coverage[key]:
                    raise ValueError("scored cell coverage changed")
        folder = output / dataset
        folder.mkdir()
        np.savez_compressed(folder / "projections.npz", **saved)
        records[dataset] = {
            "source_dataset_metadata_sha256": source["dataset_metadata_sha256"][dataset],
            "operator_input_sha256": {
                m: source_detail["outputs_sha256"]["operator_input_" + m + ".npz"]
                for m in METHOD_NAMES
            },
            "selected_projection_sha256": sha256(selected_path),
            "n_cells": plot.n_obs, "centered_max_abs_row_sum": zero_sum_error,
            "member_centering": {
                method: {
                    "nonzero_rows_before": int(np.count_nonzero(np.diff(before.indptr))),
                    "nonzero_rows_after": int(np.count_nonzero(np.diff(after.indptr))),
                    "l1_mean_before": float(row_l1_norm(before).mean()),
                    "l1_mean_after": float(row_l1_norm(after).mean()),
                }
                for method, before, after in zip(METHOD_NAMES, operators, centered_members)
            },
            "aggregation": {"current": current_diagnostics, "member_centered": new_diagnostics},
            "conditions": condition_details,
        }
        (folder / "metadata.json").write_text(json.dumps(records[dataset], indent=2) + "\n")
        pd.DataFrame(rows).to_csv(output / "scores.csv", index=False)

    frame = pd.DataFrame(rows)
    current = frame[frame.condition == "current"].set_index("dataset")
    post = frame[frame.condition == "raw_centered"].set_index("dataset")
    summaries, comparisons = [], []
    for condition in CONDITIONS:
        candidate = frame[frame.condition == condition].set_index("dataset")
        summary = {"condition": condition, "n_datasets": len(candidate)}
        for metric in METRICS:
            delta = candidate[metric] - current[metric]
            delta_post = candidate[metric] - post[metric]
            summary.update({
                metric: float(candidate[metric].mean()),
                "delta_vs_current_" + metric: float(delta.mean()),
                "improved_vs_current_" + metric: int((delta > 1e-12).sum()),
                "reduced_vs_current_" + metric: int((delta < -1e-12).sum()),
                "delta_vs_raw_centered_" + metric: float(delta_post.mean()),
                "improved_vs_raw_centered_" + metric: int((delta_post > 1e-12).sum()),
                "reduced_vs_raw_centered_" + metric: int((delta_post < -1e-12).sum()),
            })
        summaries.append(summary)
    for dataset in config["datasets"]:
        row = {"dataset": dataset}
        for condition in CONDITIONS:
            selected = frame[(frame.dataset == dataset) & (frame.condition == condition)].iloc[0]
            for metric in METRICS:
                row[condition + "_" + metric] = float(selected[metric])
        for metric in METRICS:
            row["member_minus_current_" + metric] = row["member_centered_" + metric] - row["current_" + metric]
            row["member_minus_raw_centered_" + metric] = row["member_centered_" + metric] - row["raw_centered_" + metric]
        comparisons.append(row)
    pd.DataFrame(summaries).to_csv(output / "summary.csv", index=False)
    pd.DataFrame(comparisons).to_csv(output / "comparison.csv", index=False)
    manifest = {
        "status": "complete", "finished_utc": datetime.now(timezone.utc).isoformat(),
        "datasets": list(config["datasets"]), "n_score_rows": len(rows),
        "max_control_score_error": float(frame.control_score_max_abs_error.max()),
        "max_centered_row_sum_error": max(x["centered_max_abs_row_sum"] for x in records.values()),
        "outputs_sha256": {
            str(p.relative_to(output)): sha256(p)
            for p in sorted(output.rglob("*")) if p.is_file()
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
