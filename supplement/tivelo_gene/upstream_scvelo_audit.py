#!/usr/bin/env python3
"""Audit the scVelo upstream direction used by TIvelo's path selection.

This follows the same simulation preparation as the TIvelo driver, but stops
after stochastic scVelo velocity, velocity_graph, and terminal_states. It is
read-only and is intended to identify whether the direction is already wrong
before TIvelo's path correction and gene-velocity model.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scvelo as scv


TRUE_EDGES = {
    "murk": [(f"BranchB_{i}", f"BranchB_{i + 1}") for i in range(9)],
    "branching": [(f"BranchA_{i}", f"BranchA_{i + 1}") for i in range(9)]
    + [(f"BranchB_{i}", f"BranchB_{i + 1}") for i in range(9)],
    "all": [(f"BranchA_{i}", f"BranchA_{i + 1}") for i in range(9)]
    + [(f"BranchB_{i}", f"BranchB_{i + 1}") for i in range(9)],
}
EARLY = {"murk": {"BranchB_0"}, "branching": {"BranchA_0", "BranchB_0"},
         "all": {"BranchA_0", "BranchB_0"}}


def dense(values):
    return values.toarray() if hasattr(values, "toarray") else np.asarray(values)


def row_cosine(estimate, truth):
    dot = np.sum(estimate * truth, axis=1)
    denom = np.linalg.norm(estimate, axis=1) * np.linalg.norm(truth, axis=1)
    valid = np.isfinite(dot) & np.isfinite(denom) & (denom > 1e-8)
    out = np.full(len(dot), np.nan)
    out[valid] = np.clip(dot[valid] / denom[valid], -1, 1)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "upstream_scvelo_audit")
    parser.add_argument("--seeds", nargs="+", type=int, default=[20260926, 20260927])
    parser.add_argument("--n-jobs", type=int, default=12)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from stage_audit import SCENARIOS, prepare, set_seed

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows, root_rows, velocity_rows, graph_rows = [], [], [], []
    for seed in args.seeds:
        for scenario in sorted(SCENARIOS):
            set_seed(seed)
            adata, raw_shape, _ = prepare(args.artifact_root, SCENARIOS[scenario])
            scv.tl.velocity(adata, mode="stochastic")
            scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)
            scv.tl.terminal_states(adata)
            labels = adata.obs["bin"].astype(str).to_numpy()
            groups = list(adata.obs["bin"].cat.categories)
            root_means = adata.obs.assign(_bin=labels).groupby("_bin", observed=True)["root_cells"].mean()
            end_means = adata.obs.assign(_bin=labels).groupby("_bin", observed=True)["end_points"].mean()
            root_rank = root_means.sort_values(ascending=False)
            end_rank = end_means.sort_values(ascending=False)
            summary_rows.append({
                "seed": seed, "scenario": scenario, "raw_n_cells": raw_shape[0],
                "raw_n_genes": raw_shape[1], "root_top": root_rank.index[0],
                "root_top_score": float(root_rank.iloc[0]), "end_top": end_rank.index[0],
                "end_top_score": float(end_rank.iloc[0]), "max_root_score": float(root_means.max()),
                "max_end_score": float(end_means.max()), "early_root_score_max": float(root_means[root_means.index.isin(EARLY[scenario])].max()),
                "early_end_score_max": float(end_means[end_means.index.isin(EARLY[scenario])].max()),
            })
            for group in groups:
                mask = labels == group
                root_rows.append({"seed": seed, "scenario": scenario, "bin": group,
                                  "n_cells": int(mask.sum()), "mean_root_cells": float(root_means[group]),
                                  "mean_end_points": float(end_means[group])})

            inferred = dense(adata.layers["velocity"])
            truth = dense(adata.layers["velocity_s_true"])
            cell_cos = row_cosine(inferred, truth)
            for group in groups:
                values = cell_cos[labels == group]
                velocity_rows.append({"seed": seed, "scenario": scenario, "bin": group,
                                      "mean_spliced_velocity_truth_cosine": float(np.nanmean(values)),
                                      "fraction_negative": float(np.nanmean(values < 0))})

            graph = dense(adata.uns["velocity_graph"])
            for source, target in TRUE_EDGES[scenario]:
                source_mask, target_mask = labels == source, labels == target
                forward = graph[np.ix_(source_mask, target_mask)]
                reverse = graph[np.ix_(target_mask, source_mask)]
                graph_rows.append({"seed": seed, "scenario": scenario, "source": source, "target": target,
                                   "mean_forward_graph": float(np.nanmean(forward)),
                                   "mean_reverse_graph": float(np.nanmean(reverse)),
                                   "forward_minus_reverse": float(np.nanmean(forward) - np.nanmean(reverse))})
            print(f"completed seed={seed} scenario={scenario}: root={root_rank.index[0]}", flush=True)

    pd.DataFrame(summary_rows).to_csv(args.output_dir / "scvelo_stage_summary.csv", index=False)
    pd.DataFrame(root_rows).to_csv(args.output_dir / "scvelo_root_end_by_bin.csv", index=False)
    pd.DataFrame(velocity_rows).to_csv(args.output_dir / "scvelo_velocity_truth_by_bin.csv", index=False)
    pd.DataFrame(graph_rows).to_csv(args.output_dir / "scvelo_graph_true_edge_scores.csv", index=False)
    manifest = {"purpose": "Audit upstream stochastic scVelo direction used by TIvelo; not a benchmark replacement",
                "seeds": args.seeds, "n_jobs": args.n_jobs,
                "held_fixed": "TIvelo simulation preprocessing and input H5ADs",
                "stages": ["scv.tl.velocity(mode='stochastic')", "scv.tl.velocity_graph", "scv.tl.terminal_states"]}
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summary_rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
