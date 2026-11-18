#!/usr/bin/env python3
"""Compare four/five-member geometric medians before/after row-L1 normalization."""
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
METHODS = (*METHOD_NAMES, "GraphVelo_scv")
METRICS = ("cbdir", "batc")
ZERO_TOL = 1e-12
FIELD_ATOL = 1e-9
SCORE_ATOL = 5e-7


def normalize_offdiagonal(operator):
    """Remove self weights and normalize absolute row mass; zero rows stay zero."""
    matrix = offdiagonal(operator)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    if matrix.shape[0] != matrix.shape[1] or not np.isfinite(matrix.data).all():
        raise ValueError("operator must be square and finite")
    mass = row_l1_norm(matrix)
    scale = np.divide(1., mass, out=np.zeros_like(mass), where=mass > ZERO_TOL)
    normalized = matrix.multiply(scale[:, None]).tocsr()
    normalized.eliminate_zeros()
    return normalized


def aggregate_rows(operators):
    template, values = align_method_values(operators, n_methods=len(operators))
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
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args(argv)
    validate_scvelo_runtime()
    config = load_config()
    source_path = args.operator_root / "metadata.json"
    source = json.loads(source_path.read_text())
    if tuple(source["protocol"]["methods"][:5]) != METHODS:
        raise ValueError("cached member identities/order differ")
    if source["completed_datasets"] != list(config["datasets"]):
        raise ValueError("cached datasets differ")
    if source["protocol"]["config_sha256"] != sha256(DEFAULT_CONFIG):
        raise ValueError("cached configuration changed")
    for name, expected in source["protocol"]["code_sha256"].items():
        if sha256(ROOT / name) != expected:
            raise ValueError("cached computation code changed: " + name)
    protected = {
        str(p.relative_to(ROOT)): sha256(p)
        for p in (ROOT / "ensemble/results/selected").rglob("*") if p.is_file()
    }
    code_names = [
        "ensemble/benchmark_normalization.py", "ensemble/core.py",
        "ensemble/geometric_median.py", "ensemble/pareto_members.py",
        "ensemble/metrics.py", "ensemble/reproject.py",
    ]
    output = args.output_root
    output.mkdir(parents=True, exist_ok=False)
    protocol = {
        "status": "exploratory_not_selected",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "members": list(METHODS), "member_counts": [4, 5],
        "baseline": "cached offdiagonal T native mass; cached Phi unit L1",
        "normalized": "remove diagonal, divide each member row by sum(abs(offdiagonal weights))",
        "zero_mass_tolerance": ZERO_TOL, "near_zero_rows": "remain zero",
        "aggregation": "row Euclidean geometric median, equal objective terms",
        "support": "aligned cell axes; deterministic union, absent weights zero",
        "post_aggregation_normalization": "none",
        "centering": "none",
        "diagonal_completion": "none; diagonal displacement is zero",
        "projection": "sum_j W_ij (Y_j-Y_i); raw displacements",
        "evaluation": "unchanged fixed CBDir/BATC; Forebrain PCA30/UMAP2; others UMAP2",
        "summary": "equal-dataset macro means; paired changes",
        "randomness": "none; frozen operators and coordinates",
        "source_metadata_sha256": sha256(source_path),
        "source_scores_sha256": sha256(args.operator_root / "scores.csv"),
        "config_sha256": sha256(DEFAULT_CONFIG),
        "code_sha256": {name: sha256(ROOT / name) for name in code_names},
        "repository_head": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "versions": {name: version(name) for name in
                     ("numpy", "scipy", "pandas", "anndata", "scvelo", "scikit-learn")},
        "protected_selected_sha256": protected,
        "field_atol": FIELD_ATOL, "score_atol": SCORE_ATOL,
    }
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    rows, records, diagnostics_rows = [], {}, []
    for dataset, spec in config["datasets"].items():
        print(dataset + ": validating frozen input", flush=True)
        plot, arrays, operators, detail = read_source(
            args.operator_root, source, dataset, spec)
        operators = operators[:5]
        normalized = [normalize_offdiagonal(x) for x in operators]
        bases = ("pca", "umap") if dataset == "forebrain" else ("umap",)
        saved = {
            "cell_ids": arrays["cell_ids"].astype(str),
            "labels": arrays["labels"].astype(str), "neighbors": arrays["neighbors"],
            **{"X_" + b: arrays["X_" + b] for b in bases},
        }
        norm_checks = {}
        for method, old, new in zip(METHODS, operators, normalized):
            before, after = row_l1_norm(old), row_l1_norm(new)
            active = before > ZERO_TOL
            error = float(np.max(np.abs(after - active.astype(float)), initial=0))
            expected_sign = old.multiply(active[:, None]).tocsr()
            expected_sign.eliminate_zeros()
            same_support = (np.array_equal(expected_sign.indptr, new.indptr)
                            and np.array_equal(expected_sign.indices, new.indices))
            same_sign = same_support and np.array_equal(
                np.sign(expected_sign.data), np.sign(new.data))
            if error > 1e-12 or not same_sign or np.any(new.diagonal()):
                raise ValueError("normalization invariant failed")
            change = new - old
            norm_checks[method] = {
                "before_mean": float(before.mean()), "before_min": float(before.min()),
                "before_median": float(np.median(before)), "after_mean": float(after.mean()),
                "near_zero_rows": int((~active).sum()), "max_l1_error": error,
                "sign_and_support_preserved": bool(same_sign),
                "max_entry_change": float(np.max(np.abs(change.data), initial=0)),
            }
            diagnostics_rows.append({"dataset": dataset, "method": method, **norm_checks[method]})
        if norm_checks["GraphVelo_scv"]["max_entry_change"] > 1e-12:
            raise ValueError("cached GraphVelo unit normalization unexpectedly changed")
        conditions = {}
        for count in (4, 5):
            for condition, members in (("original", operators), ("normalized", normalized)):
                key = condition + str(count)
                print(dataset + ": " + key, flush=True)
                operator, diagnostics, coefficients = aggregate_rows(members[:count])
                fields = {b: project_graphvelo(operator, arrays["X_" + b]) for b in bases}
                independent_error = 0.
                reconstruction_error = 0.
                for b in bases:
                    coords = arrays["X_" + b]
                    independent = operator @ coords - np.asarray(operator.sum(axis=1)) * coords
                    independent_error = max(independent_error, float(np.max(
                        np.abs(fields[b] - independent), initial=0)))
                    mixture = sum(
                        coefficients[:, m, None] * project_graphvelo(members[m], coords)
                        for m in range(count))
                    reconstruction_error = max(reconstruction_error, float(np.max(
                        np.abs(fields[b] - mixture), initial=0)))
                if max(independent_error, reconstruction_error) > FIELD_ATOL:
                    raise ValueError("independent projection/convex reconstruction failed")
                if np.any(operator.diagonal()) or np.any(row_l1_norm(operator) > 1 + 1e-10):
                    raise ValueError("aggregate diagonal or mass failed")
                score = score_fields(plot, fields, dataset, spec)
                if not all(np.isfinite(score[m]) for m in METRICS):
                    raise ValueError("nonfinite score")
                field_error, score_error = 0., 0.
                if condition == "original":
                    field_error = max(float(np.max(
                        np.abs(fields[b] - arrays[f"geomedian{count}_" + b]), initial=0))
                        for b in bases)
                    expected = next(r for r in detail["rows"] if r["variant"] == f"geomedian{count}")
                    score_error = max(abs(score[m] - float(expected[m])) for m in METRICS)
                    if field_error > FIELD_ATOL or score_error > SCORE_ATOL:
                        raise ValueError("frozen baseline replay failed: " + dataset + "/" + key)
                else:
                    old_score = conditions[f"original{count}"]["metrics"]
                    for coverage in ("cbdir_edge_evaluable_source_cells", "batc_evaluable_source_cells"):
                        if score[coverage] != old_score[coverage]:
                            raise ValueError("evaluation coverage changed")
                    saved[key + "_coefficients"] = coefficients
                for b in bases:
                    saved[key + "_" + b] = fields[b]
                rows.append({"dataset": dataset, "n_cells": plot.n_obs,
                             "n_members": count, "condition": condition,
                             **{m: score[m] for m in METRICS}})
                conditions[key] = {
                    "metrics": score, "aggregation": diagnostics,
                    "control_field_error": field_error, "control_score_error": score_error,
                    "independent_projection_error": independent_error,
                    "convex_projection_error": reconstruction_error,
                    "aggregate_l1_mean": float(row_l1_norm(operator).mean()),
                }
                print(json.dumps(rows[-1]), flush=True)
        folder = output / dataset
        folder.mkdir()
        np.savez_compressed(folder / "projections.npz", **saved)
        records[dataset] = {
            "n_cells": plot.n_obs, "source_metadata_sha256": source["dataset_metadata_sha256"][dataset],
            "input_sha256": {m: detail["outputs_sha256"][f"operator_input_{m}.npz"] for m in METHODS},
            "normalization": norm_checks, "conditions": conditions,
        }
        (folder / "metadata.json").write_text(json.dumps(records[dataset], indent=2) + "\n")
        pd.DataFrame(rows).to_csv(output / "scores.csv", index=False)

    frame = pd.DataFrame(rows)
    frame.groupby(["n_members", "condition"], sort=False).agg(
        n_datasets=("dataset", "size"), n_cells=("n_cells", "sum"),
        cbdir=("cbdir", "mean"), batc=("batc", "mean"),
    ).reset_index().to_csv(output / "summary.csv", index=False)
    comparisons, effects, membership = [], [], []
    for count in (4, 5):
        old = frame[(frame.n_members == count) & (frame.condition == "original")].set_index("dataset")
        new = frame[(frame.n_members == count) & (frame.condition == "normalized")].set_index("dataset")
        for dataset in config["datasets"]:
            row = {"dataset": dataset, "n_members": count}
            for metric in METRICS:
                row.update({metric + "_before": float(old.loc[dataset, metric]),
                            metric + "_after": float(new.loc[dataset, metric]),
                            metric + "_delta": float(new.loc[dataset, metric] - old.loc[dataset, metric])})
            comparisons.append(row)
        for metric in METRICS:
            delta = new[metric] - old[metric]
            effects.append({
                "n_members": count, "metric": metric, "before": float(old[metric].mean()),
                "after": float(new[metric].mean()), "delta_mean": float(delta.mean()),
                "delta_median": float(delta.median()), "delta_min": float(delta.min()),
                "delta_max": float(delta.max()), "improved": int((delta > 1e-12).sum()),
                "reduced": int((delta < -1e-12).sum()), "unchanged": int((abs(delta) <= 1e-12).sum()),
            })
    for condition in ("original", "normalized"):
        four = frame[(frame.n_members == 4) & (frame.condition == condition)].set_index("dataset")
        five = frame[(frame.n_members == 5) & (frame.condition == condition)].set_index("dataset")
        for metric in METRICS:
            delta = five[metric] - four[metric]
            membership.append({"condition": condition, "metric": metric,
                               "delta_5_minus_4": float(delta.mean()),
                               "median_delta": float(delta.median()),
                               "improved": int((delta > 1e-12).sum()),
                               "reduced": int((delta < -1e-12).sum())})
    for name, records_out in (("comparison", comparisons), ("normalization_effects", effects),
                              ("membership_effects", membership), ("input_diagnostics", diagnostics_rows)):
        pd.DataFrame(records_out).to_csv(output / (name + ".csv"), index=False)
    if any(sha256(ROOT / name) != value for name, value in protected.items()):
        raise ValueError("protected selected results changed during experiment")
    manifest = {
        "status": "complete", "finished_utc": datetime.now(timezone.utc).isoformat(),
        "datasets": list(config["datasets"]), "n_cells": sum(r["n_cells"] for r in records.values()),
        "n_score_rows": len(rows), "selected_unchanged": True,
        "max_control_field_error": max(c["control_field_error"] for r in records.values() for c in r["conditions"].values()),
        "max_control_score_error": max(c["control_score_error"] for r in records.values() for c in r["conditions"].values()),
        "max_independent_projection_error": max(c["independent_projection_error"] for r in records.values() for c in r["conditions"].values()),
        "max_convex_projection_error": max(c["convex_projection_error"] for r in records.values() for c in r["conditions"].values()),
        "outputs_sha256": {str(p.relative_to(output)): sha256(p)
                           for p in sorted(output.rglob("*")) if p.is_file()},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(effects).to_string(index=False), flush=True)
    print(pd.DataFrame(membership).to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
