#!/usr/bin/env python3
"""Screen all registered Cell configurations for dataset-level Pareto coverage."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from anndata.io import read_elem

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ensemble.core import DEFAULT_CONFIG, load_config, sha256
from ensemble.pareto_members import pareto_front, minimum_front_covers

ROOT = Path(__file__).resolve().parents[1]
METRICS = ("cbdir", "batc")


def registered_score_matrix(registry, table, datasets):
    """Join through registered IDs, preserving all methods and complete pairs."""
    methods = registry.operational_id.str.casefold().tolist()
    if len(set(methods)) != len(methods):
        raise ValueError("registered method IDs must be unique ignoring case")
    table = table.copy()
    table["method_key"] = table.method.str.casefold()
    expected = {(dataset, method) for dataset in datasets for method in methods}
    observed = list(zip(table.dataset, table.method_key))
    if len(observed) != len(set(observed)) or set(observed) != expected:
        raise ValueError("scores must cover every registered method and dataset exactly once")
    indexed = table.set_index(["dataset", "method_key"])
    values = np.asarray([[[indexed.loc[(dataset, method), metric] for metric in METRICS]
                          for method in methods] for dataset in datasets], dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("every registered score must be finite")
    return values


def inspect_native_inputs(inputs, registry, artifact_root, cached_root):
    """Inspect source headers; do not silently convert unavailable expression states."""
    records = []
    cached = {}
    metadata_path = cached_root / "metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        for dataset, digest in metadata["dataset_metadata_sha256"].items():
            path = cached_root / dataset / "metadata.json"
            if sha256(path) != digest:
                raise ValueError("cached operator metadata changed")
            detail = json.loads(path.read_text())
            evaluation = Path(load_config()["datasets"][dataset]["reference_neighbor_file"]).stem
            for method, source in detail["methods"].items():
                operator = path.parent / ("operator_input_" + method + ".npz")
                if sha256(operator) != detail["outputs_sha256"][operator.name]:
                    raise ValueError("cached member operator changed")
                cached[evaluation, method.casefold()] = source["sha256"]
            source = detail["graphvelo"]
            operator = path.parent / "operator_input_GraphVelo_scv.npz"
            if sha256(operator) != detail["outputs_sha256"][operator.name]:
                raise ValueError("cached native-Phi operator changed")
            cached[evaluation, "graphvelo_scv"] = source["source_sha256"]
    names = dict(zip(registry.config_id, registry.manuscript_name))
    for row in inputs.to_dict("records"):
        relative = Path(row["result_path"].split(":", 1)[1])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("unsafe registered result path")
        path = artifact_root / relative
        if path.stat().st_size != int(row["file_size_bytes"]):
            raise ValueError("registered result size differs: " + str(relative))
        with h5py.File(path) as handle:
            layers = list(handle["layers"])
            obs = read_elem(handle["obs"])
            var_columns = list(handle["var"])
            stored_phi = "gv" in handle.get("obsp", {})
            stored_native_neighbors = "neighbors" in handle["uns"] and "indices" in handle["uns/neighbors"]
            layer = row["velocity_layer_key"]
            if layer not in layers:
                raise ValueError("registered velocity layer missing: " + str(relative))
            native_mask = layer + "_genes"
            mask = native_mask if native_mask in var_columns else "velocity_genes" if "velocity_genes" in var_columns else ""
            if layer.startswith("velo_s_norm") and (not mask or layer + "_params" not in handle["uns"]):
                raise ValueError("MultiVelo native mask or parameters missing")
            method = row["operational_id"].casefold()
            cache = cached.get((row["evaluation_dataset"], method)) == row["sha256"]
            if cache:
                route = "verified_cached_signed_operator"
            elif method.startswith("graphvelo"):
                route = "stored_native_phi_requires_replay" if stored_phi else "native_phi_missing"
            elif "Ms" in layers:
                route = "construct_T_from_native_Ms_and_registered_velocity"
            elif method == "dynamo" and "M_s" in layers:
                route = "construct_T_from_native_M_s_and_velocity_S"
            else:
                route = "native_expression_basis_requires_alignment"
            records.append({
                "config_id": row["config_id"], "method": row["operational_id"],
                "display_name": names[row["config_id"]], "dataset": row["evaluation_dataset"],
                "result_path": row["result_path"], "registered_sha256": row["sha256"],
                "size_matches_registry": True,
                "source_hash_scope": "registered identity; only cached operator sources were hash-bound previously",
                "n_cells": len(obs), "cell_ids_unique": obs.index.is_unique,
                "n_reference": int(row["n_reference"]), "n_common": int(row["n_common"]),
                "registered_velocity_layer": layer, "native_velocity_mask": mask,
                "has_Ms": "Ms" in layers, "has_M_s": "M_s" in layers,
                "has_M_total": "M_total" in layers, "has_native_phi": stored_phi,
                "has_native_neighbor_indices": stored_native_neighbors,
                "integration_route": route,
            })
    return pd.DataFrame(records)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True,
                        help="velocity root containing the registered results directory")
    parser.add_argument("--cached-root", type=Path, default=ROOT / "ensemble/results/tsvelo_members")
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output_root.exists():
        raise ValueError("choose a new output directory")
    registry = pd.read_csv(ROOT / "methods/configurations.tsv", sep="\t", dtype=str).fillna("")
    inputs = pd.read_csv(ROOT / "evaluation/cell/inputs.tsv", sep="\t", dtype=str).fillna("")
    table_path = ROOT / "evaluation/cell/results/cell_metrics.csv"
    audit_path = ROOT / "evaluation/cell/results/cell_metrics_audit.csv"
    table = pd.read_csv(table_path)
    audit = pd.read_csv(audit_path, dtype=str).fillna("")
    config = load_config(DEFAULT_CONFIG)
    datasets = [Path(spec["reference_neighbor_file"]).stem for spec in config["datasets"].values()]
    values = registered_score_matrix(registry, table, datasets)
    if values.shape != (14, 20, 2) or len(inputs) != 280:
        raise ValueError("the formal benchmark must contain all 20 configurations on 14 datasets")
    mapped = {(row.evaluation_dataset, row.operational_id.casefold()): row for row in inputs.itertuples()}
    for row in audit.itertuples():
        original = mapped[row.dataset, row.method.casefold()]
        if row.status != "ok" or row.result_sha256 != original.sha256 or row.result_path != original.result_path:
            raise ValueError("formal score audit is not bound to the current input registry")
    if len(audit) != 280 or len(mapped) != 280:
        raise ValueError("formal input and audit coverage is incomplete")
    args.output_root.mkdir(parents=True)
    output = args.output_root
    fronts, rows = [], []
    for dataset, points in zip(datasets, values):
        front, dominates = pareto_front(points)
        fronts.append(front)
        for index, method in enumerate(registry.itertuples()):
            rows.append({
                "dataset": dataset, "config_id": method.config_id,
                "method": method.operational_id, "display_name": method.manuscript_name,
                "cbdir": points[index, 0], "batc": points[index, 1], "pareto": bool(front[index]),
                "dominated_by": "|".join(registry.iloc[np.flatnonzero(dominates[:, index])].config_id),
            })
    fronts = np.asarray(fronts)
    macro_points = values.mean(axis=0)
    macro_front = pareto_front(macro_points)[0]
    macro = pd.DataFrame([
        {"config_id": row.config_id, "method": row.operational_id,
         "display_name": row.manuscript_name,
         "cbdir": macro_points[index, 0], "batc": macro_points[index, 1],
         "macro_pareto": bool(macro_front[index]),
         "front_datasets": int(fronts[:, index].sum()),
         "sole_front_datasets": int((fronts[:, index] & (fronts.sum(axis=1) == 1)).sum()),
         "front_dataset_names": "|".join(dataset for dataset, present in zip(datasets, fronts[:, index]) if present)}
        for index, row in enumerate(registry.itertuples())])
    pd.DataFrame(rows).to_csv(output / "scores_and_fronts.csv", index=False)
    macro.to_csv(output / "method_summary.csv", index=False)
    active = np.flatnonzero(fronts.any(axis=0))
    covers = [tuple(int(active[index]) for index in ids)
              for ids in minimum_front_covers(fronts[:, active])]
    cover_rows = []
    for index, ids in enumerate(covers):
        cover_rows.append({
            "candidate": "cover_" + str(index + 1), "n_methods": len(ids),
            "config_ids": "|".join(registry.iloc[list(ids)].config_id),
            "members": "|".join(registry.iloc[list(ids)].manuscript_name),
            "covered_datasets": int(fronts[:, ids].any(axis=1).sum())})
    pd.DataFrame(cover_rows).to_csv(output / "minimum_covers.csv", index=False)
    budgets = []
    for size in range(1, len(active) + 1):
        counts = [int(fronts[:, ids].any(axis=1).sum()) for ids in combinations(active, size)]
        best = max(counts)
        budgets.append({"n_methods": size, "maximum_front_coverage": best,
                        "n_subsets_attaining_maximum": counts.count(best)})
    pd.DataFrame(budgets).to_csv(output / "coverage_by_size.csv", index=False)
    selections = {
        "original_four": ("C05", "C08", "C12", "C16"),
        "original_four_plus_graphvelo": ("C05", "C08", "C12", "C16", "C19"),
        "previous_six": ("C05", "C08", "C12", "C14", "C16", "C19"),
        "macro_front_only": tuple(registry.loc[macro_front, "config_id"]),
    }
    for row in cover_rows:
        selections[row["candidate"]] = tuple(row["config_ids"].split("|"))
    selected_rows = []
    for name, ids in selections.items():
        mask = registry.config_id.isin(ids).to_numpy()
        covered = fronts[:, mask].any(axis=1)
        selected_rows.append({
            "selection": name, "n_methods": int(mask.sum()),
            "members": "|".join(registry.loc[mask, "manuscript_name"]),
            "covered_datasets": int(covered.sum()),
            "uncovered_datasets": "|".join(dataset for dataset, present in zip(datasets, covered) if not present)})
    pd.DataFrame(selected_rows).to_csv(output / "selection_coverage.csv", index=False)
    compatibility = inspect_native_inputs(inputs, registry, args.artifact_root, args.cached_root)
    compatibility.to_csv(output / "input_compatibility.csv", index=False)
    compatibility.groupby(["config_id", "display_name", "integration_route"], sort=False).size().rename(
        "n_datasets").reset_index().to_csv(output / "compatibility_summary.csv", index=False)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(12.8, 9), layout="constrained")
    ax.imshow(fronts.T.astype(int), cmap=ListedColormap(["#F1F3F5", "#0072B2"]), aspect="auto", vmin=0, vmax=1)
    for dataset in range(14):
        for method in range(20):
            if fronts[dataset, method]:
                ax.text(dataset, method, "F*" if fronts[dataset].sum() == 1 else "F",
                        ha="center", va="center", color="white", fontsize=10)
    ax.set_yticks(range(20), [row.config_id + " " + row.manuscript_name + "  " +
                             str(int(fronts[:, index].sum())) + "/14"
                             for index, row in enumerate(registry.itertuples())])
    ax.set_xticks(range(14), datasets, rotation=45, ha="right", fontsize=9)
    ax.set_xticks(np.arange(-.5, 14, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 20, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="both", length=0)
    ax.set_title("All 20 configurations | native benchmark CBDir and BATC", fontsize=12, pad=12)
    for extension in ("png", "pdf"):
        fig.savefig(output / ("dataset_fronts." + extension), dpi=180, facecolor="white")
    plt.close(fig)
    metadata = {
        "status": "exploratory", "scope": "all 20 registered configurations; 14 formal datasets",
        "score_space": "formal method-native Cell scores; not common-displacement ensemble scores",
        "objectives": list(METRICS), "directions": "maximize both",
        "dominance": "no worse on both, strictly better on at least one; no practical-equivalence margin",
        "coverage": "at least one selected method on each dataset's front",
        "zero_front_methods_included_in_score_analysis": True,
        "cover_search": "exact on the 14 active-front methods; zero-coverage columns cannot improve minimum cardinality",
        "n_methods": 20, "n_datasets": 14, "n_scores": 280,
        "n_active_front_methods": int(len(active)),
        "minimum_cover_size": len(covers[0]), "n_minimum_covers": len(covers),
        "repository_head": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "code_sha256": {name: sha256(ROOT / name) for name in
                       ("ensemble/pareto_all_methods.py", "ensemble/pareto_members.py", "ensemble/core.py")},
        "inputs_sha256": {name: sha256(ROOT / name) for name in
                         ("methods/configurations.tsv", "evaluation/cell/inputs.tsv",
                          "evaluation/cell/results/cell_metrics.csv",
                          "evaluation/cell/results/cell_metrics_audit.csv", "ensemble/config.json")},
        "cached_operators_metadata_sha256": sha256(args.cached_root / "metadata.json"),
        "caption": "Pareto-front membership on each dataset among all 20 registered method configurations. "
                   "F denotes a nondominated method and F* the sole front method. Both native Cell metrics "
                   "are maximized; totals count datasets with front membership. Forebrain CBDir uses PCA30, "
                   "BATC uses UMAP2. Coverage is descriptive and does not imply that a fused operator "
                   "retains a member's score. Member selection is exploratory on these evaluation datasets.",
        "outputs_sha256": {path.name: sha256(path) for path in sorted(output.iterdir()) if path.is_file()},
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: metadata[key] for key in
                      ("n_methods", "n_datasets", "n_active_front_methods", "minimum_cover_size", "n_minimum_covers")}))
    print(pd.DataFrame(selected_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

