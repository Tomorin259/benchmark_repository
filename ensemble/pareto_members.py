#!/usr/bin/env python3
"""Select and test Pareto-cover ensembles from six frozen signed operators."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ensemble.benchmark_graphvelo import METHODS as GRAPHVELO_METHODS
from ensemble.core import (
    DEFAULT_CONFIG, align_method_values, complete_signed_operator,
    load_config, row_l1_norm, sha256,
)
from ensemble.geometric_median import row_geometric_median
from ensemble.metrics import score_global_metrics
from ensemble.reproject import project_graphvelo

ROOT = Path(__file__).resolve().parents[1]
METHODS = (*GRAPHVELO_METHODS, "TSvelo")
METRICS = ("cbdir", "batc")


def pareto_front(points):
    """Maximize both metrics; ties stay, with no practical-equivalence margin."""
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2 or not len(points):
        raise ValueError("Pareto input must have shape (n >= 1, 2)")
    if not np.isfinite(points).all():
        raise ValueError("Pareto input contains missing or nonfinite scores")
    differences = points[:, None, :] - points[None, :, :]
    dominates = (differences >= 0).all(axis=2) & (differences > 0).any(axis=2)
    return ~dominates.any(axis=0), dominates


def minimum_front_covers(fronts):
    """Return all smallest subsets intersecting every dataset's front."""
    fronts = np.asarray(fronts, dtype=bool)
    if fronts.ndim != 2 or not len(fronts) or not fronts.any(axis=1).all():
        raise ValueError("each dataset must have at least one front member")
    for size in range(1, fronts.shape[1] + 1):
        covers = [
            indices for indices in combinations(range(fronts.shape[1]), size)
            if fronts[:, indices].any(axis=1).all()
        ]
        if covers:
            return covers
    raise AssertionError("the complete method set must cover every front")


def score_fields(plot, fields, dataset, spec):
    for basis in ("pca", "umap") if dataset == "forebrain" else ("umap",):
        plot.obsm["ensemble_velocity_" + basis] = fields[basis]
    return score_global_metrics(plot, dataset, spec)


def read_source(source, metadata, dataset, spec):
    folder = source / dataset
    expected = metadata["dataset_metadata_sha256"][dataset]
    if sha256(folder / "metadata.json") != expected:
        raise ValueError("source metadata changed: " + dataset)
    detail = json.loads((folder / "metadata.json").read_text())
    path = folder / "projections.npz"
    if sha256(path) != detail["outputs_sha256"][path.name]:
        raise ValueError("source projections changed: " + dataset)
    # Object labels are read only after checking the accepted archive hash.
    with np.load(path, allow_pickle=True) as loaded:
        arrays = {key: loaded[key] for key in loaded.files}
    plot = ad.AnnData(obs=pd.DataFrame(
        {spec["label"]: arrays["labels"].astype(str)}, index=arrays["cell_ids"]
    ))
    plot.uns["neighbors"] = {"indices": arrays["neighbors"]}
    for basis in ("pca", "umap") if dataset == "forebrain" else ("umap",):
        plot.obsm["X_" + basis] = arrays["X_" + basis]
    operators = []
    for method in METHODS:
        path = folder / ("operator_input_" + method + ".npz")
        if sha256(path) != detail["outputs_sha256"][path.name]:
            raise ValueError("source operator changed: " + dataset + "/" + method)
        operator = sparse.load_npz(path).tocsr()
        if (operator.shape != (plot.n_obs, plot.n_obs)
                or not np.isfinite(operator.data).all()
                or np.any(operator.diagonal())
                or np.any(row_l1_norm(operator) > 1 + 1e-10)):
            raise ValueError("source operator axes or signed mass are invalid")
        operators.append(operator)
    return plot, arrays, operators, detail


def analyze_members(rows, datasets, output):
    scores = np.asarray([
        [[row[metric] for metric in METRICS] for method in METHODS
         for row in rows if row["dataset"] == dataset and row["method"] == method]
        for dataset in datasets
    ])
    if scores.shape != (len(datasets), len(METHODS), 2):
        raise ValueError("single-method scores must be a complete dataset-method matrix")
    fronts = np.asarray([pareto_front(values)[0] for values in scores])
    macro_front = pareto_front(scores.mean(axis=0))[0]
    for dataset_index, dataset in enumerate(datasets):
        _, dominates = pareto_front(scores[dataset_index])
        for method_index, method in enumerate(METHODS):
            row = next(row for row in rows
                       if row["dataset"] == dataset and row["method"] == method)
            row["pareto"] = bool(fronts[dataset_index, method_index])
            row["dominated_by"] = "|".join(
                METHODS[index] for index in np.flatnonzero(dominates[:, method_index])
            )
    pd.DataFrame(rows).to_csv(output / "single_scores.csv", index=False)
    pd.DataFrame([
        {"method": method, "cbdir": float(scores[:, index, 0].mean()),
         "batc": float(scores[:, index, 1].mean()),
         "macro_pareto": bool(macro_front[index]),
         "front_datasets": int(fronts[:, index].sum()),
         "sole_front_datasets": int((fronts[:, index] & (fronts.sum(axis=1) == 1)).sum()),
         "front_dataset_names": "|".join(
             dataset for dataset, present in zip(datasets, fronts[:, index]) if present
         )}
        for index, method in enumerate(METHODS)
    ]).to_csv(output / "single_macro.csv", index=False)

    minimum = minimum_front_covers(fronts)
    union = tuple(np.flatnonzero(fronts.any(axis=0)).tolist())
    macro = tuple(np.flatnonzero(macro_front).tolist())
    candidates = {}
    def add(indices, role):
        candidates.setdefault(tuple(indices), []).append(role)
    add(tuple(range(4)), "original_four")
    add(tuple(range(5)), "original_four_plus_graphvelo")
    add((0, 1, 2, 3, 5), "original_four_plus_tsvelo")
    add(tuple(range(6)), "all_six")
    for indices in minimum:
        add(indices, "minimum_front_cover")
    add(union, "front_union")
    add(macro, "macro_front_only")
    for removed in range(6):
        add(tuple(index for index in range(6) if index != removed),
            "leave_out_" + METHODS[removed])
    records, covers = [], []
    for size in range(1, 7):
        for indices in combinations(range(6), size):
            key = "m" + "".join(str(index + 1) for index in indices)
            covered = fronts[:, indices].any(axis=1)
            covers.append({
                "subset": key, "n_methods": size,
                "members": "|".join(METHODS[index] for index in indices),
                "covered_datasets": int(covered.sum()),
                "uncovered_datasets": "|".join(
                    dataset for dataset, ok in zip(datasets, covered) if not ok),
                "minimum_complete_cover": indices in minimum,
            })
            if indices in candidates:
                records.append({**covers[-1], "indices": list(indices),
                                "roles": candidates[indices]})
    pd.DataFrame(covers).to_csv(output / "front_coverage.csv", index=False)
    stability = []
    for left_out, dataset in enumerate(datasets):
        training_fronts = np.delete(fronts, left_out, axis=0)
        for indices in minimum_front_covers(training_fronts):
            stability.append({
                "left_out_dataset": dataset,
                "subset": "m" + "".join(str(index + 1) for index in indices),
                "members": "|".join(METHODS[index] for index in indices),
                "n_methods": len(indices),
                "left_out_front_covered": bool(fronts[left_out, indices].any()),
            })
    pd.DataFrame(stability).to_csv(output / "selection_stability.csv", index=False)
    selection = {
        "methods": list(METHODS), "datasets": datasets,
        "definition": "maximize CBDir and BATC; weakly >= both, strictly > at least one",
        "coverage_rule": "at least one selected method belongs to each dataset's front",
        "minimum_cover_size": len(minimum[0]),
        "minimum_covers": [[METHODS[index] for index in indices] for indices in minimum],
        "dataset_fronts": {
            dataset: [METHODS[index] for index in np.flatnonzero(fronts[row])]
            for row, dataset in enumerate(datasets)
        },
        "candidates": records,
        "selection_stability_scope": "descriptive leave-one-dataset-out front coverage, not external validation",
    }
    (output / "selection.json").write_text(json.dumps(selection, indent=2, sort_keys=True) + "\n")
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--stage", choices=("singles", "all"), default="all")
    args = parser.parse_args(argv)
    config = load_config(DEFAULT_CONFIG)
    source_meta = json.loads((args.source_root / "metadata.json").read_text())
    datasets = source_meta["completed_datasets"]
    if datasets != list(config["datasets"]) or source_meta["protocol"]["methods"] != list(METHODS):
        raise ValueError("the accepted six-method source must cover all configured datasets")
    code_names = ("ensemble/pareto_members.py", "ensemble/core.py",
                  "ensemble/geometric_median.py", "ensemble/reproject.py", "ensemble/metrics.py")
    for name in code_names[1:]:
        if sha256(ROOT / name) != source_meta["protocol"]["code_sha256"][name]:
            raise ValueError("computation differs from the accepted six-method source: " + name)
    if sha256(DEFAULT_CONFIG) != source_meta["protocol"]["config_sha256"]:
        raise ValueError("evaluation configuration differs from the accepted source")
    protocol = {
        "status": "exploratory", "methods": list(METHODS), "metrics": list(METRICS),
        "source_metadata_sha256": sha256(args.source_root / "metadata.json"),
        "config_sha256": sha256(DEFAULT_CONFIG),
        "code_sha256": {name: sha256(ROOT / name) for name in code_names},
        "repository_head": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "source_protocol": source_meta["protocol"],
        "aggregation": "row Euclidean geometric median; native directional mass; no post-aggregation rescaling",
        "selection": "all minimum front covers, front union, macro front, four controls and six leave-one-member-out sets",
        "selection_uses_all_evaluation_datasets": True,
        "projection": "sum_j W_ij (Y_j-Y_i); frozen reference frames; no fitting",
    }
    args.output_root.mkdir(parents=True, exist_ok=True)
    previous = args.output_root / "metadata.json"
    if previous.exists() and json.loads(previous.read_text())["protocol"] != protocol:
        raise ValueError("existing Pareto result uses different inputs or code")

    singles = []
    for dataset in datasets:
        folder = args.output_root / dataset
        folder.mkdir(exist_ok=True)
        path = folder / "single_metadata.json"
        if path.exists():
            detail = json.loads(path.read_text())
            if detail["protocol"] != protocol or sha256(folder / "singles.npz") != detail["projections_sha256"]:
                raise ValueError("existing single-method cache changed: " + dataset)
        else:
            spec = config["datasets"][dataset]
            plot, arrays, operators, source = read_source(args.source_root, source_meta, dataset, spec)
            fields = {"cell_ids": arrays["cell_ids"]}
            for basis in ("pca", "umap") if dataset == "forebrain" else ("umap",):
                fields["X_" + basis] = arrays["X_" + basis]
            rows, metrics = [], {}
            for method, operator in zip(METHODS, operators):
                bases = ("pca", "umap") if dataset == "forebrain" else ("umap",)
                projected = {basis: project_graphvelo(operator, plot.obsm["X_" + basis])
                             for basis in bases}
                score = score_fields(plot, projected, dataset, spec)
                if not np.isfinite([score[key] for key in METRICS]).all():
                    raise ValueError("nonfinite single-method score")
                rows.append({"dataset": dataset, "method": method,
                             **{key: score[key] for key in METRICS}})
                metrics[method] = score
                fields.update({method + "_" + basis: field for basis, field in projected.items()})
                print(dataset + "/" + method + ": " + json.dumps(rows[-1]), flush=True)
            np.savez_compressed(folder / "singles.npz", **fields)
            detail = {"protocol": protocol, "n_cells": plot.n_obs, "rows": rows,
                      "metrics": metrics,
                      "source_metadata_sha256": sha256(args.source_root / dataset / "metadata.json"),
                      "projections_sha256": sha256(folder / "singles.npz")}
            path.write_text(json.dumps(detail, indent=2, sort_keys=True) + "\n")
        singles.extend(dict(row) for row in detail["rows"])
    candidates = analyze_members(singles, datasets, args.output_root)
    print("SELECTION " + (args.output_root / "selection.json").read_text(), flush=True)
    if args.stage == "singles":
        previous.write_text(json.dumps({"protocol": protocol, "stage": "singles_complete"}, indent=2, sort_keys=True) + "\n")
        return 0

    ensemble_rows = []
    for dataset in datasets:
        folder = args.output_root / dataset
        spec = config["datasets"][dataset]
        plot, arrays, operators, source = read_source(args.source_root, source_meta, dataset, spec)
        template, values = align_method_values(operators, n_methods=6)
        cached = {tuple(METHODS.index(method) for method in row["members"]): row
                  for row in source["rows"] if row["variant"].startswith("geomedian")}
        fields = {"cell_ids": arrays["cell_ids"]}
        diagnostics, rows = {}, []
        for candidate in candidates:
            indices = tuple(candidate["indices"]); key = candidate["subset"]
            fields_by_basis = {}
            if indices in cached:
                cached_row = cached[indices]
                score = source["diagnostics"][cached_row["variant"]]["metrics"]
                for basis in ("pca", "umap") if dataset == "forebrain" else ("umap",):
                    fields_by_basis[basis] = arrays[cached_row["variant"] + "_" + basis]
                diagnostics[key] = {"cached_variant": cached_row["variant"]}
            else:
                if len(indices) == 1:
                    aggregate = values[indices[0]]
                    diagnostics[key] = {"aggregation": "single member"}
                else:
                    aggregate, detail, cell = row_geometric_median(
                        template, values[list(indices)], return_cell_diagnostics=True)
                    fields[key + "_coefficients"] = cell["coefficients"]
                    diagnostics[key] = {"aggregation": detail}
                operator, completion = complete_signed_operator(template, aggregate)
                diagnostics[key]["completion"] = completion
                for basis in ("pca", "umap") if dataset == "forebrain" else ("umap",):
                    fields_by_basis[basis] = project_graphvelo(operator, plot.obsm["X_" + basis])
                score = score_fields(plot, fields_by_basis, dataset, spec)
            fields.update({key + "_" + basis: field for basis, field in fields_by_basis.items()})
            diagnostics[key]["metrics"] = score
            rows.append({
                "dataset": dataset, "subset": key, "n_methods": len(indices),
                "members": candidate["members"], "roles": "|".join(candidate["roles"]),
                **{metric: score[metric] for metric in METRICS},
            })
            print(dataset + "/" + key + ": " + json.dumps(rows[-1]), flush=True)
        np.savez_compressed(folder / "ensembles.npz", **fields)
        detail = {"protocol": protocol, "rows": rows, "diagnostics": diagnostics,
                  "selection_sha256": sha256(args.output_root / "selection.json"),
                  "projections_sha256": sha256(folder / "ensembles.npz")}
        (folder / "ensemble_metadata.json").write_text(json.dumps(detail, indent=2, sort_keys=True) + "\n")
        ensemble_rows.extend(rows)
        pd.DataFrame(ensemble_rows).to_csv(args.output_root / "ensemble_scores.csv", index=False)
    frame = pd.DataFrame(ensemble_rows)
    macro = frame.groupby("subset", sort=False)[list(METRICS)].mean()
    front, _ = pareto_front(macro.to_numpy())
    summary = []
    for index, (key, row) in enumerate(macro.iterrows()):
        candidate = next(candidate for candidate in candidates if candidate["subset"] == key)
        baseline = frame[frame.subset == "m1234"].set_index("dataset")
        selected = frame[frame.subset == key].set_index("dataset")
        record = {name: candidate[name] for name in ["subset", "n_methods", "members", "covered_datasets"]}
        record["roles"] = "|".join(candidate["roles"])
        record.update({metric: float(row[metric]) for metric in METRICS})
        record["macro_pareto"] = bool(front[index])
        for metric in METRICS:
            delta = selected[metric] - baseline[metric]
            record.update({"delta_" + metric: float(delta.mean()),
                           "improved_" + metric: int((delta > 1e-12).sum()),
                           "reduced_" + metric: int((delta < -1e-12).sum())})
        summary.append(record)
    pd.DataFrame(summary).to_csv(args.output_root / "ensemble_macro.csv", index=False)
    comparisons = []
    for removed in range(6):
        key = "m" + "".join(str(index + 1) for index in range(6) if index != removed)
        for dataset in datasets:
            full = next(row for row in ensemble_rows if row["dataset"] == dataset and row["subset"] == "m123456")
            reduced = next(row for row in ensemble_rows if row["dataset"] == dataset and row["subset"] == key)
            comparisons.append({"dataset": dataset, "removed": METHODS[removed],
                                **{"delta_on_removal_" + metric: reduced[metric] - full[metric]
                                   for metric in METRICS}})
    pd.DataFrame(comparisons).to_csv(args.output_root / "leave_one_member_out.csv", index=False)
    outputs = {str(path.relative_to(args.output_root)): sha256(path)
               for path in sorted(args.output_root.rglob("*")) if path.is_file() and path != previous}
    previous.write_text(json.dumps({"protocol": protocol, "stage": "complete",
                                   "datasets": datasets, "n_single_scores": len(singles),
                                   "n_ensemble_scores": len(ensemble_rows),
                                   "outputs_sha256": outputs}, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

