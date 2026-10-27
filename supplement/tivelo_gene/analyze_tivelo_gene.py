#!/usr/bin/env python3
"""Diagnose TIvelo gene-level direction on the frozen simulation outputs.

Read the saved method H5ADs and benchmark tables to summarize directionality.
"""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


SCENARIOS = ("murk", "branching", "all")
LABEL = re.compile(r"^Branch([AB])_(\d+)$")
REQUIRED_LAYERS = ("velocity", "velocity_u")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def edge_type(source, target):
    left, right = LABEL.fullmatch(str(source)), LABEL.fullmatch(str(target))
    if left is None or right is None:
        return "unrecognized"
    if left.group(1) != right.group(1):
        return "cross_branch"
    difference = int(right.group(2)) - int(left.group(2))
    if difference == 1:
        return "forward"
    if difference == -1:
        return "reverse"
    return "skipped_or_same"


def dense(values):
    return values.toarray() if hasattr(values, "toarray") else np.asarray(values)


def cosine_rows(estimated_s, estimated_u, true_s, true_u):
    estimated_s, estimated_u = np.asarray(estimated_s, dtype=np.float64), np.asarray(estimated_u, dtype=np.float64)
    dot = estimated_s * true_s + estimated_u * true_u
    estimate_norm = np.hypot(estimated_s, estimated_u)
    truth_norm = np.hypot(true_s, true_u)
    valid = (
        np.isfinite(dot)
        & np.isfinite(estimate_norm)
        & np.isfinite(truth_norm)
        & (estimate_norm > 1e-6)
        & (truth_norm > 1e-6)
    )
    result = np.full(dot.shape, np.nan, dtype=float)
    result[valid] = np.clip(dot[valid] / (estimate_norm[valid] * truth_norm[valid]), -1, 1)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()

    repo, artifact_root, output = args.repo_root, args.artifact_root, args.output_dir
    sys.path.insert(0, str(repo))
    from evaluation.gene.simulation.reference import load_velocity_reference
    canonical_path = repo / "evaluation/gene/results/gene_cbdir_reference30_all.csv"
    cosine_path = repo / "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv"
    canonical = pd.read_csv(canonical_path)
    original = pd.read_csv(cosine_path)
    canonical = canonical[
        canonical.dataset.isin(SCENARIOS)
        & canonical.method.eq("TIvelo")
        & canonical.status.eq("ok")
    ].copy()
    original = original[original.Method.eq("TIvelo")].copy()
    if canonical.duplicated(["dataset", "gene"]).any() or original.duplicated(["Dataset", "Gene"]).any():
        raise ValueError("Duplicate TIvelo dataset-gene rows in an input table")
    scores = canonical.merge(
        original[["Dataset", "Gene", "CBDir", "CosSim"]],
        left_on=["dataset", "gene"],
        right_on=["Dataset", "Gene"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    if not scores._merge.eq("both").all():
        raise ValueError("Canonical reference-30 and original CosSim gene sets differ")
    scores["gene_class"] = scores.gene.str.replace(r"_\d+$", "", regex=True)
    scores = scores.rename(columns={"new_cbdir": "cbdir_reference30", "CBDir": "cbdir_original", "CosSim": "cossim_original"})

    score_rows, bin_rows, path_rows, child_rows, cosine_rows_out, path_summary = [], [], [], [], [], []
    h5ad_paths = []
    for scenario in SCENARIOS:
        path = artifact_root / f"indicators/gene/simulation/results_for_analyze/tivelo/tivelo_{scenario}.h5ad"
        h5ad_paths.append(path)
        if not path.is_file():
            raise FileNotFoundError(path)
        adata = ad.read_h5ad(path)
        missing = set(REQUIRED_LAYERS) - set(adata.layers)
        if missing:
            raise ValueError(f"{scenario}: missing layers {sorted(missing)}")
        if "bin" not in adata.obs or "path_dict" not in adata.uns or "child_dict" not in adata.uns:
            raise ValueError(f"{scenario}: missing bin, path_dict, or child_dict")
        paths = adata.uns["path_dict"]
        child = adata.uns["child_dict"]
        bins = adata.obs["bin"].astype(str)

        for name, nodes in paths.items():
            if name == "start_node":
                continue
            for source, target in zip(nodes[:-1], nodes[1:]):
                path_rows.append({"scenario": scenario, "path": name, "source": source,
                                  "target": target, "type": edge_type(source, target)})
        for source, targets in child.items():
            for target in targets:
                child_rows.append({"scenario": scenario, "source": source,
                                   "target": target, "type": edge_type(source, target)})

        for bin_name, frame in adata.obs.groupby(bins, observed=True):
            bin_rows.append({
                "scenario": scenario, "bin": bin_name, "n_cells": len(frame),
                "mean_root_score": pd.to_numeric(frame.root_cells, errors="coerce").mean(),
                "mean_end_score": pd.to_numeric(frame.end_points, errors="coerce").mean(),
                "mean_dpt": pd.to_numeric(frame.dpt_pseudotime, errors="coerce").mean(),
            })

        subset = scores[scores.dataset.eq(scenario)]
        genes = subset.gene.tolist()
        indices = adata.var_names.get_indexer(genes)
        if np.any(indices < 0):
            raise ValueError(f"{scenario}: H5AD lacks benchmark genes")
        estimated_s = dense(adata.layers["velocity"][:, indices])
        estimated_u = dense(adata.layers["velocity_u"][:, indices])
        true_s, true_u = load_velocity_reference(scenario, adata.obs_names, genes)
        cell_cosines = cosine_rows(estimated_s, estimated_u, true_s, true_u)
        for j, gene in enumerate(genes):
            for branch, mask in (
                ("all_cells", np.ones(adata.n_obs, dtype=bool)),
                ("BranchA", bins.str.startswith("BranchA").to_numpy()),
                ("BranchB", bins.str.startswith("BranchB").to_numpy()),
            ):
                values = cell_cosines[mask, j]
                finite = values[np.isfinite(values)]
                if not len(finite):
                    continue
                cosine_rows_out.append({
                    "scenario": scenario, "gene": gene,
                    "gene_class": re.sub(r"_\d+$", "", gene), "cell_branch": branch,
                    "n_cells": int(mask.sum()), "n_valid_vectors": len(finite),
                    "mean_cosine": float(finite.mean()),
                    "fraction_opposed": float(np.mean(finite < 0)),
                    "sign_flipped_mean_cosine_diagnostic_only": float(-finite.mean()),
                })
        score_main = np.asarray(adata.uns.get("score_main", []), dtype=float)
        child_types = [r["type"] for r in child_rows if r["scenario"] == scenario]
        path_summary.append({
            "scenario": scenario, "n_cells": adata.n_obs, "n_genes": adata.n_vars,
            "start_node": str(paths["start_node"]),
            "main_path_first": str(paths["main_path"][0]),
            "main_path_last": str(paths["main_path"][-1]),
            "mean_score_main": float(np.nanmean(score_main)) if len(score_main) else np.nan,
            "median_score_main": float(np.nanmedian(score_main)) if len(score_main) else np.nan,
            "child_forward": child_types.count("forward"),
            "child_reverse": child_types.count("reverse"),
            "child_cross_branch": child_types.count("cross_branch"),
            "child_other": len(child_types) - sum(child_types.count(k) for k in ("forward", "reverse", "cross_branch")),
        })
        print(f"{scenario}: root={paths['start_node']}, child edges forward/reverse="
              f"{child_types.count('forward')}/{child_types.count('reverse')}", flush=True)

    output.mkdir(parents=True, exist_ok=True)
    scores[["dataset", "gene", "gene_class", "cbdir_reference30", "cbdir_original", "cossim_original"]].to_csv(
        output / "tivelo_simulation_gene_scores.csv", index=False
    )
    summary = scores.groupby(["dataset", "gene_class"], observed=True).agg(
        n_genes=("gene", "size"), mean_cbdir_reference30=("cbdir_reference30", "mean"),
        median_cbdir_reference30=("cbdir_reference30", "median"),
        fraction_negative_cbdir=("cbdir_reference30", lambda x: np.mean(x < 0)),
        mean_cossim=("cossim_original", "mean"),
        fraction_negative_cossim=("cossim_original", lambda x: np.mean(x < 0)),
    ).reset_index()
    summary.to_csv(output / "tivelo_simulation_summary.csv", index=False)
    pd.DataFrame(bin_rows).to_csv(output / "root_end_dpt_by_bin.csv", index=False)
    pd.DataFrame(path_rows).to_csv(output / "inferred_path_edges.csv", index=False)
    pd.DataFrame(child_rows).to_csv(output / "directed_neighbor_cluster_edges.csv", index=False)
    pd.DataFrame(cosine_rows_out).to_csv(output / "velocity_truth_cosine_by_gene_branch.csv", index=False)
    pd.DataFrame(path_summary).to_csv(output / "path_summary.csv", index=False)

    cosine_check = pd.DataFrame(cosine_rows_out).query("cell_branch == 'all_cells'")
    cosine_check = cosine_check.merge(scores[["dataset", "gene", "cossim_original"]],
                                      left_on=["scenario", "gene"], right_on=["dataset", "gene"],
                                      validate="one_to_one")
    cosine_check["absolute_difference"] = (cosine_check.mean_cosine - cosine_check.cossim_original).abs()
    cosine_check[["scenario", "gene", "mean_cosine", "cossim_original", "absolute_difference"]].to_csv(
        output / "cosine_recalculation_check.csv", index=False
    )
    real = pd.read_csv(canonical_path)
    real = real[~real.dataset.isin(SCENARIOS) & real.status.eq("ok")].copy()
    real_summary = real.groupby(["dataset", "method"], observed=True).agg(
        n_genes=("gene", "size"), mean_cbdir_reference30=("new_cbdir", "mean")
    ).reset_index()
    real_summary["rank_within_dataset"] = real_summary.groupby("dataset").mean_cbdir_reference30.rank(
        method="min", ascending=False
    )
    real_summary.to_csv(output / "real_method_context.csv", index=False)

    manifest = {
        "purpose": "Exploratory TIvelo direction diagnosis",
        "expected_bin_direction": "BranchA_0 to BranchA_9; BranchB_0 to BranchB_9",
        "inputs": [
            {"locator": f"artifact://velocity/indicators/gene/simulation/results_for_analyze/tivelo/tivelo_{scenario}.h5ad",
             "bytes": path.stat().st_size}
            for scenario, path in zip(SCENARIOS, h5ad_paths)
        ],
        "tables": [
            {"path": str(path.relative_to(repo)), "sha256": sha256(path)}
            for path in (canonical_path, cosine_path)
        ],
        "max_absolute_cosine_recalculation_difference": float(cosine_check.absolute_difference.max()),
    }
    (output / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)
    print("maximum CosSim recalculation difference:", manifest["max_absolute_cosine_recalculation_difference"], flush=True)


if __name__ == "__main__":
    main()
