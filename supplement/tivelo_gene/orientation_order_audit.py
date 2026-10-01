#!/usr/bin/env python3
"""Test whether branch mixing and ordering affect TIvelo orientation scores.

The TIvelo orientation score treats an ordered list of cells as one line and
smooths it with a 100-cell convolution. This read-only audit applies the same
`orient_score` to true bin-ordered branch-specific lists and to concatenated
lists, providing a mechanism control for the observed path reversal.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tivelo.direction.model import orient_score
from tivelo.utils.velocity_genes import compute_velocity_genes


def score_gene(mu, ms, order, index):
    return float(orient_score(mu[order, index], ms[order, index])[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "orientation_order_audit")
    parser.add_argument("--n-jobs", type=int, default=12)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for scenario in ("murk", "branching", "all"):
        path = args.artifact_root / f"indicators/gene/simulation/results_for_analyze/tivelo/tivelo_{scenario}.h5ad"
        adata = ad.read_h5ad(path)
        mask = compute_velocity_genes(adata, n_top_genes=2000, r2_adjust=True, inplace=False)
        indices = np.flatnonzero(mask)
        names = np.asarray(adata.var_names)[indices]
        mu, ms = np.asarray(adata.layers["Mu"]), np.asarray(adata.layers["Ms"])
        labels = adata.obs["bin"].astype(str).to_numpy()
        parsed = [(str(x), int(re.search(r"_(\d+)$", str(x)).group(1))) for x in labels]
        branch_a_only = np.array(sorted(np.flatnonzero(np.char.startswith(labels.astype(str), "BranchA")), key=lambda i: (parsed[i][1], i)), dtype=int)
        branch_b_only = np.array(sorted(np.flatnonzero(np.char.startswith(labels.astype(str), "BranchB")), key=lambda i: (parsed[i][1], i)), dtype=int)
        true_concat = np.concatenate([branch_a_only, branch_b_only])
        # This deliberately mixes the two lineages in a single ordered list.
        mixed_by_bin = np.array(sorted(np.arange(len(labels)), key=lambda i: (parsed[i][1], 0 if labels[i].startswith("BranchA") else 1, i)))
        orders = {"branch_A_forward": branch_a_only, "branch_B_forward": branch_b_only,
                  "branches_concatenated": true_concat, "time_interleaved": mixed_by_bin}
        for order_name, order in orders.items():
            if len(order) == 0:
                continue
            values = Parallel(n_jobs=args.n_jobs)(delayed(score_gene)(mu, ms, order, index) for index in indices)
            frame = pd.DataFrame({"gene": names, "score": values})
            frame["gene_class"] = frame.gene.str.replace(r"_\d+$", "", regex=True)
            for gene_class, group in frame.groupby("gene_class", observed=True):
                rows.append({"scenario": scenario, "order": order_name, "gene_class": gene_class,
                             "n_genes": len(group), "mean_orientation_score": float(group.score.mean()),
                             "median_orientation_score": float(group.score.median()),
                             "fraction_negative": float(np.mean(group.score < 0))})
            print(f"completed {scenario} / {order_name}", flush=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "orientation_score_by_order.csv", index=False)


if __name__ == "__main__":
    main()
