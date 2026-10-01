#!/usr/bin/env python3
"""Audit the root decision and PAGA topology used before TIvelo path building."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scvelo as scv


def dense(value):
    return value.toarray() if hasattr(value, "toarray") else np.asarray(value)


def edge_type(source: str, target: str):
    left, right = source.split("_"), target.split("_")
    if left[0] != right[0]:
        return "cross_branch"
    delta = int(right[1]) - int(left[1])
    return "forward" if delta > 0 else "reverse" if delta < 0 else "same_bin"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "root_decision_audit")
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from stage_audit import SCENARIOS, prepare, set_seed
    from tivelo.path.prune import graph_prune
    from tivelo.path.start_node import get_start_node

    args.output_dir.mkdir(parents=True, exist_ok=True)
    root_rows, edge_rows, summary_rows = [], [], []
    for scenario in ("branching", "all", "murk"):
        set_seed(args.seed)
        adata, raw_shape, _ = prepare(args.artifact_root, SCENARIOS[scenario])
        select_info, select_weight = graph_prune(adata, "bin", "X_umap")
        # Let TIvelo's own selector run the upstream scVelo calls exactly once.
        # Read the scores only after this call; re-running velocity_graph can
        # change stochastic/parallel outputs and would decouple the table from
        # the reported selected root.
        chosen = get_start_node(adata, select_info, "bin", njobs=12, mode="stochastic")
        labels = adata.obs["bin"].astype(str).to_numpy()
        root = adata.obs.assign(_bin=labels).groupby("_bin", observed=True)["root_cells"].mean()
        end = adata.obs.assign(_bin=labels).groupby("_bin", observed=True)["end_points"].mean()
        for group in select_info.index.astype(str):
            root_rows.append({"scenario": scenario, "bin": group,
                              "mean_root_cells": float(root[group]),
                              "mean_end_points": float(end[group]),
                              "root_minus_end": float(root[group] - end[group]),
                              "root_rank": int(root.rank(ascending=False, method="min")[group]),
                              "end_rank": int(end.rank(ascending=False, method="min")[group])})
        names = [str(x) for x in select_info.index]
        matrix = dense(select_weight)
        for i in range(matrix.shape[0]):
            for j in range(i + 1, matrix.shape[1]):
                if matrix[i, j] > 0:
                    edge_rows.append({"scenario": scenario, "source": names[i], "target": names[j],
                                      "edge_type": edge_type(names[i], names[j]),
                                      "paga_tree_weight": float(matrix[i, j])})
        top_root = root.sort_values(ascending=False).iloc[0]
        top_end = end.sort_values(ascending=False).iloc[0]
        summary_rows.append({"scenario": scenario, "selected_root": str(chosen),
                             "top_root_bin": str(root.idxmax()), "top_root_score": float(top_root),
                             "second_root_score": float(root.sort_values(ascending=False).iloc[1]),
                             "top_end_bin": str(end.idxmax()), "top_end_score": float(top_end),
                             "selected_root_root_minus_end": float(root[chosen] - end[chosen]),
                             "n_paga_nodes": len(names), "n_paga_edges": int(np.sum(matrix > 0) // 2)})
        print(f"completed {scenario}: chosen={chosen}", flush=True)

    pd.DataFrame(root_rows).to_csv(args.output_dir / "root_end_scores.csv", index=False)
    pd.DataFrame(edge_rows).to_csv(args.output_dir / "paga_tree_edges.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(args.output_dir / "root_decision_summary.csv", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({
        "purpose": "Read-only audit of TIvelo root selection before process_path",
        "seed": args.seed,
        "stages": ["mirrored preprocessing", "graph_prune", "scv stochastic velocity", "velocity_graph", "terminal_states", "get_start_node"],
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
