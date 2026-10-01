#!/usr/bin/env python3
"""Compare TIvelo orientation scores under true and TIvelo-selected path orders.

This is a read-only mechanism audit. It reuses the saved TIvelo H5AD outputs,
keeps the observed Mu/Ms layers fixed, and compares target synthetic gene
classes under branch-forward orders with the final TIvelo path order.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from tivelo.direction.model import orient_score


TARGET_CLASSES = {
    "murk": ["OneSidedBoost"],
    "branching": ["MultiForward", "BoostALL"],
    "all": ["MultiForward", "BoostALL", "OneSidedBoost", "OneSidedBoost2"],
}


def ordered_branch_indices(labels: np.ndarray):
    bins = np.array([int(x.rsplit("_", 1)[1]) for x in labels])
    idx = np.arange(len(labels))
    a = idx[np.array([x.startswith("BranchA") for x in labels])]
    b = idx[np.array([x.startswith("BranchB") for x in labels])]
    a = a[np.argsort(bins[a], kind="stable")]
    b = b[np.argsort(bins[b], kind="stable")]
    return a, b


def target_indices(names: np.ndarray, gene_class: str):
    return np.flatnonzero(np.array([x.startswith(gene_class + "_") for x in names]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "orientation_order_audit")
    parser.add_argument("--scenarios", nargs="+", choices=sorted(TARGET_CLASSES), default=sorted(TARGET_CLASSES))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for scenario in args.scenarios:
        path = args.artifact_root / f"indicators/gene/simulation/results_for_analyze/tivelo/tivelo_{scenario}.h5ad"
        adata = ad.read_h5ad(path)
        labels = np.asarray(adata.obs["bin"].astype(str).to_numpy(), dtype=str)
        names = np.asarray(adata.var_names).astype(str)
        mu = np.asarray(adata.layers["Mu"])
        ms = np.asarray(adata.layers["Ms"])
        a_forward, b_forward = ordered_branch_indices(labels)
        dpt = np.asarray(adata.obs["dpt_pseudotime"], dtype=float)
        for path_name, nodes in adata.uns["path_dict"].items():
            if path_name != "main_path":
                continue
            mask = np.isin(labels, np.asarray(nodes).astype(str))
            final_order = np.flatnonzero(mask)
            final_order = final_order[np.argsort(dpt[final_order], kind="stable")]
        orders = {"branch_A_forward": a_forward, "branch_B_forward": b_forward,
                  "tivelo_final_main_path": final_order}
        for order_name, order in orders.items():
            if len(order) == 0:
                continue
            for gene_class in TARGET_CLASSES[scenario]:
                indices = target_indices(names, gene_class)
                values = np.array([float(orient_score(mu[order, j], ms[order, j])[0]) for j in indices])
                rows.append({"scenario": scenario, "order": order_name, "gene_class": gene_class,
                             "n_genes": len(values), "n_cells": len(order),
                             "mean_orientation_score": float(values.mean()),
                             "median_orientation_score": float(np.median(values)),
                             "fraction_negative": float(np.mean(values < 0))})
            print(f"completed {scenario} / {order_name}", flush=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "orientation_score_by_order_full_target.csv", index=False)


if __name__ == "__main__":
    main()
