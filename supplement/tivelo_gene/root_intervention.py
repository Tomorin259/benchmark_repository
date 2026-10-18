#!/usr/bin/env python3
"""Test whether TIvelo's inferred root/directed neighbors explain simulation errors.

Use the simulator's known early bins to set the root while retaining the
preprocessed cells, expression layers, TIvelo optimizer and evaluation metric.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from tivelo.velocity.DTI import get_child_dict, get_d_nn
from tivelo.velocity.model import get_velocity

from analyze_tivelo_gene import cosine_rows, dense


SCENARIOS = ("branching", "all")
LABEL = re.compile(r"^Branch([AB])_(\d+)$")


def directed_neighbor_types(d_nn, labels):
    source, target = np.nonzero(dense(d_nn))
    parsed = [LABEL.fullmatch(str(label)) for label in labels]
    branch = np.array([match.group(1) for match in parsed])
    time = np.array([int(match.group(2)) for match in parsed])
    same_branch = branch[source] == branch[target]
    delta = time[target] - time[source]
    return {
        "n_directed_neighbor_pairs": len(source),
        "n_forward_bin_pairs": int(np.sum(same_branch & (delta > 0))),
        "n_reverse_bin_pairs": int(np.sum(same_branch & (delta < 0))),
        "n_same_bin_pairs": int(np.sum(same_branch & (delta == 0))),
        "n_cross_branch_pairs": int(np.sum(~same_branch)),
    }


def score_velocity(adata, genes, velocity_s, velocity_u, neighbors, edges, score_gene):
    indices = adata.var_names.get_indexer(genes)
    if np.any(indices < 0):
        raise ValueError("A benchmark gene is absent from the baseline H5AD")
    cosines = cosine_rows(
        velocity_s[:, indices], velocity_u[:, indices],
        dense(adata.layers["velocity_s_true"][:, indices]),
        dense(adata.layers["velocity_u_true"][:, indices]),
    )
    bins = adata.obs["bin"].astype(str)
    rows = []
    for col, gene in enumerate(genes):
        index = indices[col]
        cbdir, _, _, _ = score_gene(
            dense(adata.layers["Ms"][:, index]), dense(adata.layers["Mu"][:, index]),
            velocity_s[:, index], velocity_u[:, index], bins.to_numpy(), neighbors, edges,
        )
        for branch, mask in (
            ("all_cells", np.ones(adata.n_obs, dtype=bool)),
            ("BranchA", bins.str.startswith("BranchA").to_numpy()),
            ("BranchB", bins.str.startswith("BranchB").to_numpy()),
        ):
            values = cosines[mask, col]
            valid = values[np.isfinite(values)]
            if len(valid):
                rows.append({
                    "gene": gene, "gene_class": re.sub(r"_\d+$", "", gene),
                    "cell_branch": branch, "n_valid_vectors": len(valid),
                    "mean_cossim": float(valid.mean()),
                    "cbdir_reference30": float(cbdir),
                    "fraction_opposed": float(np.mean(valid < 0)),
                })
    return rows


def fit_velocity(adata, d_nn, knn, same_c, seed):
    np.random.seed(seed)
    torch.manual_seed(seed)
    _, v_u, v_s = get_velocity(
        adata, d_nn, knn, same_c=same_c, n_epochs=100,
        loss_fun="mse", only_s=False, constrain=True,
        alpha_1=1, alpha_2=0.1, batch_size=1024,
        filter_genes=False,
    )
    return v_s, v_u


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    sys.path.insert(0, str(args.repo_root / "evaluation/gene"))
    from run_reference30_gene_cbdir import (
        SIM_CONFIGS, align_reference_neighbors, read_reference_asset, score_gene,
    )

    source = args.repo_root / "evaluation/gene/results/gene_cbdir_reference30_all.csv"
    table = pd.read_csv(source)
    table = table[table.method.eq("TIvelo") & table.dataset.isin(SCENARIOS) & table.status.eq("ok")]
    if table.duplicated(["dataset", "gene"]).any():
        raise ValueError("Duplicate TIvelo gene in benchmark table")

    all_scores, all_neighbors = [], []
    for scenario in SCENARIOS:
        h5ad = args.artifact_root / f"indicators/gene/simulation/results_for_analyze/tivelo/tivelo_{scenario}.h5ad"
        original = ad.read_h5ad(h5ad)
        genes = table.loc[table.dataset.eq(scenario), "gene"].tolist()
        labels = original.obs["bin"].astype(str).to_numpy()
        reference_names, reference_indices = read_reference_asset(
            args.repo_root / f"evaluation/gene/results/reference_neighbors/{scenario}.npz"
        )
        neighbors = align_reference_neighbors(original, reference_names, reference_indices)
        edges = SIM_CONFIGS[scenario]["edges"]
        knn = (dense(original.obsp["distances"]) != 0).astype(float)
        baseline_dnn = dense(original.obsp["d_nn"])
        all_neighbors.append({"scenario": scenario, "variant": "baseline_archived_graph",
                              **directed_neighbor_types(baseline_dnn, labels)})
        archived_s = dense(original.layers["velocity"])
        archived_u = dense(original.layers["velocity_u"])
        for row in score_velocity(original, genes, archived_s, archived_u, neighbors, edges, score_gene):
            all_scores.append({"scenario": scenario, "variant": "baseline_archived_velocity", **row})

        for variant in ("baseline_refit", "oracle_root_BranchA_0", "oracle_root_BranchB_0"):
            work = original.copy()
            if variant == "baseline_refit":
                d_nn = baseline_dnn
            else:
                root = variant.removeprefix("oracle_root_")
                child, _, _ = get_child_dict(work, "bin", {"start_node": root}, adjust=True)
                d_nn, _, _, _ = get_d_nn(work, "bin", child, root, root_select="connectivities")
            all_neighbors.append({"scenario": scenario, "variant": variant,
                                  **directed_neighbor_types(d_nn, labels)})
            same_c = (labels[:, None] == labels[None, :]).astype(float)
            v_s, v_u = fit_velocity(work, d_nn, knn, same_c, args.seed)
            for row in score_velocity(original, genes, v_s, v_u, neighbors, edges, score_gene):
                all_scores.append({"scenario": scenario, "variant": variant, **row})
            print(f"completed {scenario} / {variant}", flush=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.DataFrame(all_scores)
    scores.to_csv(args.output_dir / "root_intervention_gene_scores.csv", index=False)
    neighbors = pd.DataFrame(all_neighbors)
    neighbors.to_csv(args.output_dir / "root_intervention_neighbor_pairs.csv", index=False)
    summary = scores[scores.cell_branch.eq("all_cells")].groupby(
        ["scenario", "variant", "gene_class"], observed=True
    ).agg(n_genes=("gene", "size"), mean_cossim=("mean_cossim", "mean"),
          median_cossim=("mean_cossim", "median"),
          mean_cbdir_reference30=("cbdir_reference30", "mean"),
          median_cbdir_reference30=("cbdir_reference30", "median"),
          fraction_genes_negative=("mean_cossim", lambda x: np.mean(x < 0)),
          fraction_genes_negative_cbdir=("cbdir_reference30", lambda x: np.mean(x < 0)),
          mean_fraction_opposed=("fraction_opposed", "mean")).reset_index()
    summary.to_csv(args.output_dir / "root_intervention_summary.csv", index=False)
    manifest = {
        "type": "exploratory oracle-root intervention, not a benchmark correction",
        "method_settings": {
            "n_epochs": 100, "loss_fun": "mse", "only_s": False, "constrain": True,
            "alpha_1": 1, "alpha_2": 0.1, "batch_size": 1024, "filter_genes": False,
            "seed_per_variant": args.seed,
        },
        "held_fixed": "archived TIvelo-preprocessed H5AD, cell and gene set, Mu/Ms, KNN, optimizer, reference-30 neighbors and metrics",
        "changed": "root used to construct child_dict, directed neighbors and DPT; oracle roots known from simulation",
        "inputs": [
            f"artifact://velocity/indicators/gene/simulation/results_for_analyze/tivelo/tivelo_{scenario}.h5ad"
            for scenario in SCENARIOS
        ],
    }
    (args.output_dir / "root_intervention_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)
    print(neighbors.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
