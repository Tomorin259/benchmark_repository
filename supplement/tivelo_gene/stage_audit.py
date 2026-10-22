#!/usr/bin/env python3
"""Audit the TIvelo simulation call chain and intermediate direction scores.

This diagnostic exposes preprocessing, process_path, correct_path,
directed-neighbor construction and velocity fitting as separate stages.
See methods/tivelo.py for method execution and methods/source_manifest.json
for retained source identities.
"""

from __future__ import annotations

import argparse
import json
import random
import re
import sys
from copy import deepcopy
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch


SCENARIOS = {
    "branching": "indicators/gene/simulation/branching_analyze/dataset/s1_branch.h5ad",
    "murk": "indicators/gene/simulation/murk_analyze/data/s1_only_murk.h5ad",
    "all": "indicators/gene/simulation/sim_all_analyze/data/sim_all.h5ad",
}
LABEL = re.compile(r"^Branch([AB])_(\d+)$")


def dense(values):
    return values.toarray() if hasattr(values, "toarray") else np.asarray(values)


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def serial_path(path_dict):
    return {str(key): [str(value) for value in values] if hasattr(values, "__iter__") and not isinstance(values, str)
            else str(values) for key, values in path_dict.items()}


def classify_edge(source, target):
    left, right = LABEL.fullmatch(str(source)), LABEL.fullmatch(str(target))
    if left is None or right is None:
        return "unrecognized"
    if left.group(1) != right.group(1):
        return "cross_branch"
    delta = int(right.group(2)) - int(left.group(2))
    if delta > 0:
        return "forward"
    if delta < 0:
        return "reverse"
    return "same_bin"


def path_rows(scenario, stage, path_dict):
    rows = []
    for name, nodes in path_dict.items():
        if name == "start_node":
            continue
        for source, target in zip(nodes[:-1], nodes[1:]):
            rows.append({"scenario": scenario, "stage": stage, "path": name,
                         "source": str(source), "target": str(target),
                         "edge_type": classify_edge(source, target)})
    return rows


def child_rows(scenario, child_dict):
    rows = []
    for source, targets in child_dict.items():
        for target in targets:
            rows.append({"scenario": scenario, "source": str(source), "target": str(target),
                         "edge_type": classify_edge(source, target)})
    return rows


def dnn_summary(scenario, d_nn, labels):
    source, target = np.nonzero(dense(d_nn))
    left = [LABEL.fullmatch(str(labels[i])) for i in source]
    right = [LABEL.fullmatch(str(labels[i])) for i in target]
    forward = reverse = same = cross = 0
    for lmatch, rmatch in zip(left, right):
        if lmatch is None or rmatch is None:
            continue
        if lmatch.group(1) != rmatch.group(1):
            cross += 1
        else:
            delta = int(rmatch.group(2)) - int(lmatch.group(2))
            forward += delta > 0
            reverse += delta < 0
            same += delta == 0
    return {"scenario": scenario, "n_pairs": int(len(source)),
            "forward_bin_pairs": int(forward), "reverse_bin_pairs": int(reverse),
            "same_bin_pairs": int(same), "cross_branch_pairs": int(cross)}


def cosine_rows(estimated_s, estimated_u, true_s, true_u):
    dot = estimated_s * true_s + estimated_u * true_u
    estimate_norm = np.hypot(estimated_s, estimated_u)
    truth_norm = np.hypot(true_s, true_u)
    valid = np.isfinite(dot) & np.isfinite(estimate_norm) & np.isfinite(truth_norm)
    valid &= (estimate_norm > 1e-6) & (truth_norm > 1e-6)
    result = np.full(dot.shape, np.nan, dtype=float)
    result[valid] = np.clip(dot[valid] / (estimate_norm[valid] * truth_norm[valid]), -1, 1)
    return result


def prepare(artifact_root, rel_path):
    input_path = artifact_root / rel_path
    adata = sc.read(input_path)
    adata.var_names_make_unique()
    # Prepare the saved simulation input with the recorded TIvelo preprocessing.
    if "spliced" not in adata.layers:
        if "X_spliced" in adata.layers:
            adata.layers["spliced"] = adata.layers["X_spliced"]
            adata.layers["unspliced"] = adata.layers["X_unspliced"]
        elif "new" in adata.layers:
            adata.layers["unspliced"] = adata.layers["new"].copy()
            adata.layers["spliced"] = adata.layers["total"] - adata.layers["new"]
            if "velocity_N" in adata.layers and "velocity_T" in adata.layers:
                adata.layers["velocity"] = adata.layers["velocity_N"] + adata.layers["velocity_T"]
    elif "clusters" not in adata.obs:
        for candidate in ["celltype", "cell_type", "phase", "Annotation", "time_str",
                          "Clusters", "cell_cycle_phase"]:
            if candidate in adata.obs:
                adata.obs["clusters"] = adata.obs[candidate].astype(str).astype("category")
                adata.obs["bin"] = adata.obs["clusters"].copy()
                break
    if "X_umap" not in adata.obsm:
        for candidate in ("X_RFP_GFP", "X_tsne", "X_cdr", "X_pca"):
            if candidate in adata.obsm:
                adata.obsm["X_umap"] = adata.obsm[candidate].copy()
                break
    raw_shape = list(adata.shape)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    if "X_umap" not in adata.obsm:
        for candidate in ("X_cdr", "X_pca", "X_tsne"):
            if candidate in adata.obsm:
                adata.obsm["X_umap"] = adata.obsm[candidate].copy()
                break
    return adata, raw_shape, input_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "stage_audit")
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--scenarios", nargs="+", choices=sorted(SCENARIOS), default=sorted(SCENARIOS))
    args = parser.parse_args()

    sys.path.insert(0, str(args.repo_root / "evaluation/gene"))
    from run_reference30_gene_cbdir import SIM_CONFIGS, align_reference_neighbors, read_reference_asset, score_gene
    from tivelo.direction.correct import correct_path
    from tivelo.path.process import process_path
    from tivelo.utils.velocity_genes import compute_velocity_genes
    from tivelo.velocity.DTI import get_child_dict, get_d_nn
    from tivelo.velocity.model import get_velocity

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows, root_rows, path_rows_all, child_rows_all = [], [], [], []
    orientation_rows, stage_score_rows, dnn_rows = [], [], []
    manifest = {"purpose": "TIvelo call-chain mechanism audit",
                "seed": args.seed, "settings": {"group_key": "bin", "emb_key": "X_umap",
                "start_mode": "stochastic", "rev_stat": "mean", "t1": 0.1, "t2": 2,
                "loss_fun": "mse", "only_s": False, "constrain": True, "alpha_1": 1,
                "alpha_2": 0.1, "batch_size": 1024, "n_epochs": 100, "filter_genes": False},
                "inputs": []}

    for scenario in args.scenarios:
        set_seed(args.seed)
        adata, raw_shape, input_path = prepare(args.artifact_root, SCENARIOS[scenario])
        manifest["inputs"].append({"scenario": scenario, "locator": f"artifact://velocity/{SCENARIOS[scenario]}",
                                   "raw_shape": raw_shape, "preprocessed_shape": list(adata.shape),
                                   "n_velocity_layers": int("Ms" in adata.layers and "Mu" in adata.layers)})
        print(f"{scenario}: preprocessing complete {raw_shape} -> {adata.shape}", flush=True)

        raw_path, ax = process_path(adata, "bin", "X_umap", njobs=12, start_mode="stochastic")
        plt.close("all")
        for bin_name, frame in adata.obs.groupby(adata.obs["bin"].astype(str), observed=True):
            root_rows.append({"scenario": scenario, "stage": "process_path", "bin": bin_name,
                              "n_cells": len(frame), "mean_root_score": float(pd.to_numeric(frame.root_cells).mean()),
                              "mean_end_score": float(pd.to_numeric(frame.end_points).mean())})
        path_rows_all.extend(path_rows(scenario, "process_path", raw_path))
        raw_main = [str(x) for x in raw_path["main_path"]]

        velocity_mask = compute_velocity_genes(adata, n_top_genes=2000, r2_adjust=True, inplace=False)
        velocity_names = np.asarray(adata.var_names)[velocity_mask]
        before_correct = deepcopy(raw_path)
        corrected_path, ax = correct_path(raw_path, adata, "bin", rev_stat="mean", root_select="connectivities")
        plt.close("all")
        path_rows_all.extend(path_rows(scenario, "correct_path", corrected_path))
        for idx, gene in zip(np.flatnonzero(velocity_mask), velocity_names):
            orientation_rows.append({"scenario": scenario, "gene": str(gene), "selected": True,
                                     "score_main": float(adata.uns["score_main"][list(velocity_names).index(gene)])
                                     if len(adata.uns.get("score_main", [])) == len(velocity_names) else np.nan,
                                     "score_branch_1": float(adata.uns.get("score_branch_1", [np.nan] * len(velocity_names))[list(velocity_names).index(gene)])
                                     if len(adata.uns.get("score_branch_1", [])) == len(velocity_names) else np.nan})
        for gene in np.asarray(adata.var_names)[~velocity_mask]:
            orientation_rows.append({"scenario": scenario, "gene": str(gene), "selected": False,
                                     "score_main": np.nan, "score_branch_1": np.nan})

        selected_before = set(str(x) for x in velocity_names)
        score_main = np.asarray(adata.uns.get("score_main", []), dtype=float)
        score_branch = np.asarray(adata.uns.get("score_branch_1", []), dtype=float)
        child_dict, level_dict, threshold_list = get_child_dict(
            adata, "bin", corrected_path, threshold=0.1, threshold_trans=2, adjust=False
        )
        dti_start = list(level_dict.keys())[0]
        d_nn, knn, child_c, same_c = get_d_nn(adata, "bin", child_dict, dti_start,
                                               emb_key="X_umap", root_select="connectivities")
        plt.close("all")
        child_rows_all.extend(child_rows(scenario, child_dict))
        dnn_rows.append(dnn_summary(scenario, d_nn, adata.obs["bin"].astype(str).to_numpy()))

        reference_names, reference_indices = read_reference_asset(
            args.repo_root / f"evaluation/gene/results/reference_neighbors/{scenario}.npz"
        )
        neighbors = align_reference_neighbors(adata, reference_names, reference_indices)
        labels = adata.obs["bin"].astype(str).to_numpy()
        genes = SIM_CONFIGS[scenario]["genes"]
        v_adata, v_u, v_s = get_velocity(
            adata, dense(d_nn), dense(knn), same_c=dense(same_c), n_epochs=100,
            loss_fun="mse", only_s=False, constrain=True, alpha_1=1, alpha_2=0.1,
            batch_size=1024, filter_genes=False,
        )
        true_s = dense(adata.layers["velocity_s_true"])
        true_u = dense(adata.layers["velocity_u_true"])
        indices = adata.var_names.get_indexer(genes)
        cell_cos = cosine_rows(v_s[:, indices], v_u[:, indices], true_s[:, indices], true_u[:, indices])
        for j, gene in enumerate(genes):
            cbdir, _, _, _ = score_gene(
                dense(adata.layers["Ms"][:, indices[j]]), dense(adata.layers["Mu"][:, indices[j]]),
                v_s[:, indices[j]], v_u[:, indices[j]], labels, neighbors, SIM_CONFIGS[scenario]["edges"]
            )
            finite = cell_cos[:, j][np.isfinite(cell_cos[:, j])]
            stage_score_rows.append({"scenario": scenario, "gene": gene,
                                     "selected_for_orientation": gene in selected_before,
                                     "cbdir_reference30": float(cbdir),
                                     "cossim": float(np.mean(finite)),
                                     "fraction_cell_cosine_negative": float(np.mean(finite < 0))})
        summary_rows.append({
            "scenario": scenario, "raw_n_cells": raw_shape[0], "raw_n_genes": raw_shape[1],
            "preprocessed_n_cells": adata.n_obs, "preprocessed_n_genes": adata.n_vars,
            "n_velocity_genes": len(velocity_names), "velocity_target_genes_selected": int(sum(g in selected_before for g in genes)),
            "process_start_node": str(before_correct["start_node"]), "process_main_first": raw_main[0],
            "process_main_last": raw_main[-1], "correct_start_node": str(corrected_path["start_node"]),
            "correct_main_first": str(corrected_path["main_path"][0]), "correct_main_last": str(corrected_path["main_path"][-1]),
            "dti_start_node": str(dti_start), "mean_score_main": float(np.nanmean(score_main)) if len(score_main) else np.nan,
            "median_score_main": float(np.nanmedian(score_main)) if len(score_main) else np.nan,
            "fraction_orientation_negative": float(np.mean(score_main < 0)) if len(score_main) else np.nan,
            "mean_score_branch_1": float(np.nanmean(score_branch)) if len(score_branch) else np.nan,
            "median_score_branch_1": float(np.nanmedian(score_branch)) if len(score_branch) else np.nan,
            "n_dti_edges": len(child_rows(scenario, child_dict)),
        })
        print(f"{scenario}: process={before_correct['start_node']} -> correct={corrected_path['start_node']} -> dti={dti_start}; model complete", flush=True)

    pd.DataFrame(summary_rows).to_csv(args.output_dir / "stage_summary.csv", index=False)
    pd.DataFrame(root_rows).to_csv(args.output_dir / "root_end_scores_at_process_path.csv", index=False)
    pd.DataFrame(path_rows_all).to_csv(args.output_dir / "path_edges_by_stage.csv", index=False)
    pd.DataFrame(child_rows_all).to_csv(args.output_dir / "child_edges_after_correct_path.csv", index=False)
    pd.DataFrame(orientation_rows).to_csv(args.output_dir / "orientation_scores.csv", index=False)
    pd.DataFrame(dnn_rows).to_csv(args.output_dir / "dnn_direction_summary.csv", index=False)
    pd.DataFrame(stage_score_rows).to_csv(args.output_dir / "stage_final_gene_scores.csv", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summary_rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
