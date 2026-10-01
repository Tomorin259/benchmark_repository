#!/usr/bin/env python3
"""Compare mean- and median-based TIvelo path orientation decisions."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "orientation_stat_intervention")
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from stage_audit import SCENARIOS, prepare, set_seed
    from tivelo.direction.correct import correct_path
    from tivelo.path.process import process_path

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for scenario in ("branching", "all", "murk"):
        set_seed(args.seed)
        base, _, _ = prepare(args.artifact_root, SCENARIOS[scenario])
        observed, _ = process_path(base, "bin", "X_umap", njobs=12, start_mode="stochastic")
        for rev_stat in ("mean", "median"):
            adata = base.copy()
            corrected, _ = correct_path(copy.deepcopy(observed), adata, "bin", rev_stat=rev_stat,
                                        root_select="connectivities")
            scores = np.asarray(adata.uns["score_main"], dtype=float)
            rows.append({"scenario": scenario, "rev_stat": rev_stat,
                         "mean_score": float(scores.mean()), "median_score": float(np.median(scores)),
                         "fraction_positive": float(np.mean(scores > 0)),
                         "returned_start_node": str(corrected["start_node"]),
                         "returned_main_first": str(corrected["main_path"][0]),
                         "returned_main_last": str(corrected["main_path"][-1]),
                         "n_velocity_genes": len(scores)})
            print(f"completed {scenario} / {rev_stat}", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(args.output_dir / "orientation_stat_summary.csv", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({
        "purpose": "Read-only comparison of TIvelo orientation reversal statistics",
        "seed": args.seed,
        "held_fixed": "mirrored preprocessing and process_path output",
        "changed": "correct_path rev_stat only",
    }, indent=2) + "\n")
    print(frame.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
