#!/usr/bin/env python3
"""Compare TIvelo's process_path root, corrected root, and oracle roots.

This is a mechanism experiment, not a benchmark variant. It follows the
simulation preprocessing and path calls, then fits the same TIvelo velocity
model under each root while holding the prepared data, graph, optimizer, and
reference-30 scoring fixed.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import torch


LABEL = re.compile(r"^Branch([AB])_(\d+)$")


def dense(values):
    return values.toarray() if hasattr(values, "toarray") else np.asarray(values)


def edge_counts(d_nn, labels):
    source, target = np.nonzero(dense(d_nn))
    values = {"forward": 0, "reverse": 0, "same_bin": 0, "cross_branch": 0}
    for left, right in zip(labels[source], labels[target]):
        lmatch, rmatch = LABEL.fullmatch(str(left)), LABEL.fullmatch(str(right))
        if lmatch is None or rmatch is None:
            continue
        if lmatch.group(1) != rmatch.group(1):
            values["cross_branch"] += 1
        else:
            delta = int(rmatch.group(2)) - int(lmatch.group(2))
            values["forward" if delta > 0 else "reverse" if delta < 0 else "same_bin"] += 1
    return {"n_dnn_pairs": int(len(source)), **values}


def cosine_rows(v_s, v_u, t_s, t_u):
    dot = v_s * t_s + v_u * t_u
    vnorm, tnorm = np.hypot(v_s, v_u), np.hypot(t_s, t_u)
    valid = np.isfinite(dot) & np.isfinite(vnorm) & np.isfinite(tnorm) & (vnorm > 1e-6) & (tnorm > 1e-6)
    out = np.full(dot.shape, np.nan)
    out[valid] = np.clip(dot[valid] / (vnorm[valid] * tnorm[valid]), -1, 1)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "root_stage_intervention")
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from stage_audit import SCENARIOS, dense as stage_dense, prepare, set_seed
    sys.path.insert(0, str(args.repo_root / "evaluation/gene"))
    from run_reference30_gene_cbdir import SIM_CONFIGS, align_reference_neighbors, read_reference_asset, score_gene
    from tivelo.direction.correct import correct_path
    from tivelo.path.process import process_path
    from tivelo.velocity.DTI import get_child_dict, get_d_nn
    from tivelo.velocity.model import get_velocity

    args.output_dir.mkdir(parents=True, exist_ok=True)
    score_rows, edge_rows, path_rows = [], [], []
    manifest = {"type": "root-stage intervention; not a benchmark correction", "seed": args.seed,
                "held_fixed": "simulation input, mirrored preprocessing, reference-30 neighbors, gene metric, TIvelo optimizer",
                "variants": ["process_path_root", "correct_path_root", "oracle_early_root"]}

    for scenario in sorted(SCENARIOS):
        set_seed(args.seed)
        adata, raw_shape, input_path = prepare(args.artifact_root, SCENARIOS[scenario])
        process_dict, _ = process_path(adata, "bin", "X_umap", njobs=12, start_mode="stochastic")
        process_root = str(process_dict["start_node"])
        corrected_dict, _ = correct_path(deepcopy(process_dict), adata, "bin", rev_stat="mean",
                                         root_select="connectivities")
        corrected_root = str(corrected_dict["start_node"])
        labels = adata.obs["bin"].astype(str).to_numpy()
        genes = SIM_CONFIGS[scenario]["genes"]
        ref_names, ref_indices = read_reference_asset(args.repo_root / f"evaluation/gene/results/reference_neighbors/{scenario}.npz")
        neighbors = align_reference_neighbors(adata, ref_names, ref_indices)
        truth_s = stage_dense(adata.layers["velocity_s_true"])
        truth_u = stage_dense(adata.layers["velocity_u_true"])
        roots = [("process_path_root", process_root), ("correct_path_root", corrected_root)]
        early = "BranchB_0" if scenario == "murk" else "BranchA_0"
        roots.append(("oracle_early_root", early))

        for variant, root in roots:
            work = adata.copy()
            child, level, thresholds = get_child_dict(work, "bin", {"start_node": root}, threshold=0.1,
                                                      threshold_trans=2, adjust=False)
            dti_root = str(list(level.keys())[0])
            d_nn, knn, child_c, same_c = get_d_nn(work, "bin", child, dti_root, emb_key="X_umap",
                                                  root_select="connectivities")
            set_seed(args.seed)
            _, v_u, v_s = get_velocity(work, dense(d_nn), dense(knn), same_c=dense(same_c), n_epochs=100,
                                       loss_fun="mse", only_s=False, constrain=True, alpha_1=1, alpha_2=0.1,
                                       batch_size=1024, filter_genes=False)
            indices = work.var_names.get_indexer(genes)
            cossim = cosine_rows(v_s[:, indices], v_u[:, indices], truth_s[:, indices], truth_u[:, indices])
            edge_rows.append({"scenario": scenario, "variant": variant, "requested_root": root,
                              "dti_root": dti_root, **edge_counts(d_nn, labels)})
            for j, gene in enumerate(genes):
                cbdir, _, _, _ = score_gene(stage_dense(work.layers["Ms"][:, indices[j]]),
                                             stage_dense(work.layers["Mu"][:, indices[j]]), v_s[:, indices[j]],
                                             v_u[:, indices[j]], labels, neighbors, SIM_CONFIGS[scenario]["edges"])
                valid = cossim[:, j][np.isfinite(cossim[:, j])]
                score_rows.append({"scenario": scenario, "variant": variant, "requested_root": root,
                                   "dti_root": dti_root, "gene": gene, "cbdir_reference30": float(cbdir),
                                   "cossim": float(np.mean(valid)), "fraction_cell_cossim_negative": float(np.mean(valid < 0))})
            print(f"completed {scenario} / {variant} ({root})", flush=True)

    scores = pd.DataFrame(score_rows)
    scores.to_csv(args.output_dir / "root_stage_gene_scores.csv", index=False)
    edges = pd.DataFrame(edge_rows)
    edges.to_csv(args.output_dir / "root_stage_dnn_summary.csv", index=False)
    summary = scores.groupby(["scenario", "variant"], observed=True).agg(
        n_genes=("gene", "size"), mean_cbdir_reference30=("cbdir_reference30", "mean"),
        mean_cossim=("cossim", "mean"), fraction_negative_cbdir=("cbdir_reference30", lambda x: np.mean(x < 0)),
        fraction_negative_cossim=("cossim", lambda x: np.mean(x < 0)),
        mean_fraction_cell_cossim_negative=("fraction_cell_cossim_negative", "mean"),
    ).reset_index()
    summary.to_csv(args.output_dir / "root_stage_summary.csv", index=False)
    manifest["scenarios"] = sorted(SCENARIOS)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)
    print(edges.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
