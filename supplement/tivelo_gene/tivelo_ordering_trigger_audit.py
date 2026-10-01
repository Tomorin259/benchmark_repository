#!/usr/bin/env python3
"""Isolate cell ordering in TIvelo's orientation score, without model fitting.

Each paired comparison uses identical cells, Mu/Ms and selected genes. Only
the ordering changes: TIvelo's DPT ordering versus known simulation bin order.
The process-main paths are read from the existing fixed-seed stage audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed, parallel_config
from scipy.stats import spearmanr
from tivelo.direction.model import orient_score
from tivelo.utils.velocity_genes import compute_velocity_genes

from stage_audit import SCENARIOS, prepare, set_seed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--n-jobs", type=int, default=12)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parent
    edges = pd.read_csv(root / "stage_audit/path_edges_by_stage.csv")
    stage = pd.read_csv(root / "stage_audit/stage_summary.csv")
    gene_rows, cell_rows, order_rows, input_rows = [], [], [], []

    for scenario in ("murk", "branching", "all"):
        set_seed(args.seed)
        adata, _, _ = prepare(args.artifact_root, SCENARIOS[scenario])
        selected = compute_velocity_genes(adata, n_top_genes=2000, r2_adjust=True, inplace=False)
        selected_genes = np.asarray(adata.var_names)[selected].astype(str)
        e = edges[(edges.scenario == scenario) & (edges.stage == "process_path") & (edges.path == "main_path")]
        nodes = list(e.source) + [e.iloc[-1].target]
        start = str(stage.loc[stage.scenario == scenario, "process_start_node"].iloc[0])
        variants = {
            "process_main": (nodes, start),
            "single_B": ([f"BranchB_{i}" for i in range(10)], "BranchB_0"),
        }
        input_rows.append({"scenario": scenario, "input": SCENARIOS[scenario],
                           "Mu_sha256": hashlib.sha256(np.asarray(adata.layers["Mu"]).tobytes()).hexdigest(),
                           "Ms_sha256": hashlib.sha256(np.asarray(adata.layers["Ms"]).tobytes()).hexdigest(),
                           "n_selected_genes": int(np.sum(selected))})
        for path_name, (nodes, start) in variants.items():
            sub = adata[adata.obs["bin"].isin(nodes), :].copy()
            sub.uns["iroot"] = int(np.flatnonzero(sub.obs["bin"].astype(str) == start)[0])
            sc.tl.diffmap(sub)
            sc.tl.dpt(sub)
            labels = sub.obs["bin"].astype(str).to_numpy()
            bins = np.array([int(x.rsplit("_", 1)[1]) for x in labels])
            branches = np.array([x.split("_")[0] for x in labels])
            dpt = sub.obs["dpt_pseudotime"].to_numpy()
            if not np.isfinite(dpt).all():
                raise ValueError(f"Nonfinite DPT: {scenario}/{path_name}")
            mu, ms = np.asarray(sub.layers["Mu"])[:, selected], np.asarray(sub.layers["Ms"])[:, selected]
            orders = {"dpt": np.argsort(dpt), "known_bin": np.argsort(bins, kind="stable")}
            for branch in np.unique(branches):
                mask = branches == branch
                order_rows.append({"scenario": scenario, "path": path_name, "root": start,
                                   "branch": branch, "n_cells": int(np.sum(mask)),
                                   "spearman_dpt_bin": float(spearmanr(dpt[mask], bins[mask]).statistic)})
            for i, cell in enumerate(sub.obs_names):
                cell_rows.append({"scenario": scenario, "path": path_name, "cell": str(cell),
                                  "branch": branches[i], "bin": int(bins[i]), "dpt": float(dpt[i])})
            for order_name, order in orders.items():
                with parallel_config(backend="loky", inner_max_num_threads=1):
                    scores = Parallel(n_jobs=args.n_jobs)(
                        delayed(orient_score)(mu[order, j], ms[order, j]) for j in range(len(selected_genes)))
                for j, (gene, result) in enumerate(zip(selected_genes, scores)):
                    gene_rows.append({"scenario": scenario, "path": path_name, "order": order_name,
                                      "gene": gene, "gene_class": gene.rsplit("_", 1)[0],
                                      "score": float(result[0])})
                mean = float(np.mean([x[0] for x in scores]))
                print(f"{scenario}/{path_name}/{order_name}: mean={mean:.6f}", flush=True)
                pd.DataFrame(gene_rows).to_csv(args.output_dir / "gene_orientation_scores.csv", index=False)
                pd.DataFrame(cell_rows).to_csv(args.output_dir / "cell_orders.csv", index=False)

    gene_scores = pd.DataFrame(gene_rows)
    gene_scores.to_csv(args.output_dir / "gene_orientation_scores.csv", index=False)
    summary = gene_scores.groupby(["scenario", "path", "order"], observed=True).score.agg(
        n_genes="size", mean="mean", median="median", fraction_negative=lambda x: np.mean(x < 0)).reset_index()
    summary.to_csv(args.output_dir / "orientation_summary.csv", index=False)
    gene_scores.groupby(["scenario", "path", "order", "gene_class"], observed=True).score.agg(
        n_genes="size", mean="mean").reset_index().to_csv(args.output_dir / "class_summary.csv", index=False)
    pd.DataFrame(cell_rows).to_csv(args.output_dir / "cell_orders.csv", index=False)
    pd.DataFrame(order_rows).to_csv(args.output_dir / "dpt_bin_correlations.csv", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({
        "purpose": "Ordering-only orientation diagnostic; no velocity fit or benchmark changes",
        "seed": args.seed, "inputs": input_rows,
        "held_fixed": ["cells within each paired path", "Mu", "Ms", "selected genes", "orient_score defaults"],
        "changed": "cell ordering only",
        "known_bin_order": "bin number, then original row index; not exact within-bin kinetic time",
        "process_main_paths": "stage_audit/path_edges_by_stage.csv",
    }, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
