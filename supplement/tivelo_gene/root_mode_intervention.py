#!/usr/bin/env python3
"""Compare TIvelo's stochastic and dynamical root-inference modes."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "root_mode_intervention")
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--scenarios", nargs="+", default=["branching", "all", "murk"])
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from stage_audit import SCENARIOS, prepare, set_seed
    from tivelo.path.prune import graph_prune
    from tivelo.path.start_node import get_start_node

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for scenario in args.scenarios:
        for mode in ("stochastic", "dynamical"):
            set_seed(args.seed)
            adata, raw_shape, _ = prepare(args.artifact_root, SCENARIOS[scenario])
            select_info, _ = graph_prune(adata, "bin", "X_umap")
            chosen = get_start_node(adata, select_info, "bin", njobs=12, mode=mode)
            labels = adata.obs["bin"].astype(str).to_numpy()
            root = adata.obs.assign(_bin=labels).groupby("_bin", observed=True)["root_cells"].mean()
            end = adata.obs.assign(_bin=labels).groupby("_bin", observed=True)["end_points"].mean()
            rank = root.sort_values(ascending=False)
            rows.append({"scenario": scenario, "mode": mode, "selected_root": str(chosen),
                         "top_root": str(rank.index[0]), "top_root_score": float(rank.iloc[0]),
                         "second_root": str(rank.index[1]), "second_root_score": float(rank.iloc[1]),
                         "top_end": str(end.idxmax()), "top_end_score": float(end.max()),
                         "selected_root_score": float(root[chosen]),
                         "selected_root_end_score": float(end[chosen])})
            print(f"completed {scenario} / {mode}: {chosen}", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(args.output_dir / "root_mode_summary.csv", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({
        "purpose": "Read-only comparison of TIvelo root-inference modes",
        "seed": args.seed,
        "modes": ["stochastic", "dynamical"],
        "not_a_benchmark_variant": True,
    }, indent=2) + "\n")
    print(frame.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
