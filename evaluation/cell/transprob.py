#!/usr/bin/env python3
"""Calculate cell-level TransProb/TransCosine from explicit result inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd

from metrics import (
    DEFAULT_CONFIG,
    DEFAULT_INPUTS,
    ensure_neighbors,
    ensure_velocity_graph,
    load_config,
    prepare_labels,
    prepare_velocity,
    resolve_artifact,
)


def load_dependencies():
    import scvelo as scv
    from tivelo.utils.metrics import cross_boundary_scvelo_probs

    return scv, cross_boundary_scvelo_probs


def evaluate_one(artifact_root, row, config, scv, cross_boundary_scvelo_probs):
    dataset = row["dataset"]
    path = resolve_artifact(artifact_root, row["result_path"])
    adata = ad.read_h5ad(path)
    prepare_velocity(adata, row)
    label = prepare_labels(adata, dataset, config)
    ensure_neighbors(adata, scv, config["metrics"]["n_neighbors"])
    ensure_velocity_graph(
        adata,
        scv,
        n_jobs=config["metrics"]["velocity_graph_n_jobs"],
    )
    _, score = cross_boundary_scvelo_probs(
        adata,
        cluster_key=label,
        cluster_edges=config["datasets"][dataset]["cluster_edges"],
        trans_g_key="velocity_graph",
    )
    return {
        "dataset": dataset,
        "method": row["operational_id"],
        "transprob": float(score),
        "result_path": row["result_path"],
        "velocity_layer_key": row["velocity_layer_key"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--inputs", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--methods", nargs="*")
    parser.add_argument("--datasets", nargs="*")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.check_only == (args.output is not None):
        parser.error("use exactly one of --output or --check-only")

    config = load_config(args.config)
    from metrics import read_result_inputs

    rows = read_result_inputs(args.inputs)
    if args.methods:
        rows = rows[rows["operational_id"].isin(args.methods)]
    if args.datasets:
        rows = rows[rows["dataset"].isin(args.datasets)]
    if rows.empty:
        raise ValueError("No result inputs selected")
    artifact_root = args.artifact_root.resolve()
    if args.check_only:
        for _, row in rows.iterrows():
            resolve_artifact(artifact_root, row["result_path"])
        print(f"Validated {len(rows)} result inputs")
        return

    scv, cross_boundary_scvelo_probs = load_dependencies()
    results = []
    for _, row in rows.iterrows():
        print(f"[{row['dataset']}] {row['operational_id']}", flush=True)
        results.append(
            evaluate_one(
                artifact_root,
                row,
                config,
                scv,
                cross_boundary_scvelo_probs,
            )
        )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = ["dataset", "method", "transprob"]
    pd.DataFrame(results).loc[:, columns].to_csv(output, index=False)
    print(f"Wrote {len(results)} TransProb rows to {output}")


if __name__ == "__main__":
    main()
