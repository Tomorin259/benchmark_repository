#!/usr/bin/env python3
"""Score a prepared AnnData object with the official TiVelo 0.1.4 functions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
from tivelo.utils.metrics import (
    cross_boundary_correctness2,
    cross_boundary_scvelo_probs,
)


def serialise_raw(raw: dict) -> dict:
    return {
        f"{source}->{target}": [float(value) for value in values]
        for (source, target), values in raw.items()
    }


def aggregate(raw: dict) -> float:
    values = [np.mean(item) for item in raw.values() if len(item)]
    values = [value for value in values if np.isfinite(value)]
    return float(np.mean(values)) if values else float("nan")


def normalize_cluster_edges(adata: ad.AnnData, cluster_key: str, edges: list) -> list:
    """Match manifest edge labels to the serialized AnnData label dtype."""
    observed = list(adata.obs[cluster_key].dropna().unique())
    by_text = {str(value): value for value in observed}
    normalized = []
    for source, target in edges:
        try:
            normalized.append([by_text[str(source)], by_text[str(target)]])
        except KeyError as exc:
            raise ValueError(
                f"{adata.uns['evaluation_metadata']['dataset']}: cluster edge "
                f"{source!r}->{target!r} is absent from {cluster_key!r}"
            ) from exc
    return normalized


def json_ready(value):
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-input", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    if args.dataset not in config["datasets"]:
        raise KeyError(f"unknown dataset {args.dataset!r}")
    adata = ad.read_h5ad(args.prepared_input)
    cluster_key = adata.uns["evaluation_metadata"]["label_key"]
    edges = normalize_cluster_edges(
        adata, cluster_key, config["datasets"][args.dataset]["cluster_edges"]
    )
    highdim = cross_boundary_correctness2(
        adata,
        cluster_key=cluster_key,
        velocity_key="velocity",
        cluster_edges=edges,
        return_raw=True,
    )
    transcosine = cross_boundary_scvelo_probs(
        adata,
        cluster_key=cluster_key,
        cluster_edges=edges,
        trans_g_key="velocity_graph",
        return_raw=True,
    )
    payload = {
        "highdim_cbdir": aggregate(highdim),
        "transcosine": aggregate(transcosine),
        "highdim_cbdir_raw": serialise_raw(highdim),
        "transcosine_raw": serialise_raw(transcosine),
        "tivelo_version": "0.1.4",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "highdim_cbdir": payload["highdim_cbdir"],
                "transcosine": payload["transcosine"],
                "tivelo_version": payload["tivelo_version"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
