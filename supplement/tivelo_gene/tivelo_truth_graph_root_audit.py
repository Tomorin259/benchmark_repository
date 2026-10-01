#!/usr/bin/env python3
"""Diagnostic control: TIvelo/scVelo root scores with kinetic-equation velocity.

Uses the archived TIvelo AnnData only as a read-only input. It replaces the
upstream scVelo velocity layer in memory with the simulator's spliced-velocity
truth, while retaining the archived moments, graph neighbors, and the original
scVelo velocity-gene subset. This is a diagnostic, not a benchmark variant.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scvelo as scv
from scvelo.tools.transition_matrix import transition_matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    score_rows = []
    flux_rows = []
    for scenario in ("murk", "branching", "all"):
        input_path = args.artifact_root / f"indicators/gene/simulation/results_for_analyze/tivelo/tivelo_{scenario}.h5ad"
        adata = ad.read_h5ad(input_path)
        labels = adata.obs["bin"].astype(str).to_numpy()
        selected_genes = adata.var_names[np.asarray(adata.var["velocity_genes"], bool)]
        adata.layers["velocity_truth"] = np.asarray(adata.layers["velocity_s_true"])
        adata.uns["velocity_truth_params"] = {"mode": "stochastic"}
        scv.tl.velocity_graph(
            adata,
            vkey="velocity_truth",
            xkey="Ms",
            gene_subset=selected_genes,
            sqrt_transform=True,
            n_jobs=1,
        )
        scv.tl.terminal_states(adata, vkey="velocity_truth")
        T = transition_matrix(adata, vkey="velocity_truth", backward=False, self_transitions=False).tocsr()
        for source in sorted(np.unique(labels)):
            source_rows = np.flatnonzero(labels == source)
            branch, bin_id = source.rsplit("_", 1)
            bin_id = int(bin_id)
            mass = T[source_rows, :].sum(axis=0).A.ravel() / len(source_rows)
            probs = {"forward": 0.0, "reverse": 0.0, "same_bin": 0.0, "cross_branch": 0.0}
            for target in np.unique(labels):
                target_branch, target_bin = target.rsplit("_", 1)
                p = float(mass[labels == target].sum())
                if target_branch != branch:
                    probs["cross_branch"] += p
                elif int(target_bin) > bin_id:
                    probs["forward"] += p
                elif int(target_bin) < bin_id:
                    probs["reverse"] += p
                else:
                    probs["same_bin"] += p
            score_rows.append(
                {
                    "scenario": scenario,
                    "bin": source,
                    "n_cells": len(source_rows),
                    "mean_root_cells": float(np.mean(adata.obs["root_cells"].to_numpy()[source_rows])),
                    "mean_end_points": float(np.mean(adata.obs["end_points"].to_numpy()[source_rows])),
                    "n_graph_genes": len(selected_genes),
                }
            )
            flux_rows.append({"scenario": scenario, "bin": source, **probs})
        print(f"done {scenario}", flush=True)

    pd.DataFrame(score_rows).to_csv(args.output_dir / "truth_velocity_root_end_by_bin.csv", index=False)
    pd.DataFrame(flux_rows).to_csv(args.output_dir / "truth_velocity_transition_flux_by_bin.csv", index=False)


if __name__ == "__main__":
    main()
