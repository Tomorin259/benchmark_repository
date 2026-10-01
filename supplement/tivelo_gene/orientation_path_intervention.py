#!/usr/bin/env python3
"""Audit TIvelo orientation scores under controlled path definitions.

The input layers and preprocessing are held fixed.  Only the ordered
cluster path supplied to ``correct_path`` is changed, allowing us to separate
path-order effects from the upstream velocity and from the later optimizer.
This is exploratory and never writes benchmark artifacts.
"""

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
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "orientation_path_intervention")
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from stage_audit import SCENARIOS, prepare, set_seed
    from tivelo.direction.correct import correct_path
    from tivelo.path.process import process_path
    from tivelo.utils.velocity_genes import compute_velocity_genes

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    path_rows = []
    for scenario in ("branching", "all", "murk"):
        set_seed(args.seed)
        base, _, _ = prepare(args.artifact_root, SCENARIOS[scenario])
        observed, _ = process_path(base, "bin", "X_umap", njobs=12, start_mode="stochastic")
        variants = {"process_path": observed}
        if scenario in {"branching", "all"}:
            variants["true_A_forward"] = {"start_node": "BranchA_0",
                                           "main_path": [f"BranchA_{i}" for i in range(10)]}
            variants["true_B_forward"] = {"start_node": "BranchB_0",
                                           "main_path": [f"BranchB_{i}" for i in range(10)]}
        else:
            variants["true_B_forward"] = {"start_node": "BranchB_0",
                                           "main_path": [f"BranchB_{i}" for i in range(10)]}

        for variant, path_dict in variants.items():
            adata = base.copy()
            corrected, _ = correct_path(copy.deepcopy(path_dict), adata, "bin", rev_stat="mean",
                                        root_select="connectivities")
            scores = np.asarray(adata.uns["score_main"], dtype=float)
            velocity_mask = compute_velocity_genes(adata, n_top_genes=2000, r2_adjust=True, inplace=False)
            names = np.asarray(adata.var_names)[velocity_mask].astype(str)
            classes = sorted(set(x.rsplit("_", 1)[0] for x in names
                                 if "_" in x and x.rsplit("_", 1)[1].isdigit()))
            for gene_class in classes:
                indices = np.flatnonzero(np.array([x.startswith(gene_class + "_") for x in names]))
                if len(indices) == 0:
                    continue
                values = scores[indices]
                rows.append({"scenario": scenario, "variant": variant,
                             "gene_class": gene_class, "n_genes": len(values),
                             "mean_orientation_score": float(values.mean()),
                             "median_orientation_score": float(np.median(values)),
                             "fraction_negative": float(np.mean(values < 0)),
                             "returned_start_node": str(corrected["start_node"]),
                             "returned_main_first": str(corrected["main_path"][0]),
                             "returned_main_last": str(corrected["main_path"][-1])})
            path_rows.append({"scenario": scenario, "variant": variant,
                              "input_start_node": str(path_dict["start_node"]),
                              "input_main_path": "|".join(map(str, path_dict["main_path"])),
                              "returned_start_node": str(corrected["start_node"]),
                              "returned_main_path": "|".join(map(str, corrected["main_path"]))})
            print(f"completed {scenario} / {variant}", flush=True)

    pd.DataFrame(rows).to_csv(args.output_dir / "orientation_by_controlled_path.csv", index=False)
    pd.DataFrame(path_rows).to_csv(args.output_dir / "controlled_paths.csv", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({
        "purpose": "Path-only TIvelo orientation audit; not a benchmark replacement",
        "seed": args.seed,
        "held_fixed": "mirrored simulation preprocessing and saved input layers",
        "changed": "ordered path supplied to tivelo.direction.correct.correct_path",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
