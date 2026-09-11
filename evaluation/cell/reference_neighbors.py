#!/usr/bin/env python3
"""Load and apply the fixed reference neighbor graphs used by cell metrics."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import sparse


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MODULE_DIR / "config.json"


def load_reference_params(config_path: Path = DEFAULT_CONFIG) -> dict:
    """Load the single fixed-neighbor contract used by Cell evaluation."""
    with config_path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    reference = config.get("reference")
    if not isinstance(reference, dict):
        raise ValueError(f"{config_path}: missing reference contract")
    required = {
        "n_neighbors",
        "n_pcs",
        "use_rep",
        "metric",
        "random_state",
        "n_recurse_neighbors",
        "mode_neighbors",
    }
    missing = required.difference(reference)
    if missing:
        raise ValueError(f"{config_path}: missing reference fields {sorted(missing)}")
    return {
        "min_shared_counts": 20,
        "n_top_genes_requested": 2000,
        "n_pcs": int(reference["n_pcs"]),
        "n_neighbors": int(reference["n_neighbors"]),
        "use_rep": str(reference["use_rep"]),
        "metric": str(reference["metric"]),
        "random_state": int(reference["random_state"]),
        "n_recurse_neighbors": int(reference["n_recurse_neighbors"]),
        "mode_neighbors": str(reference["mode_neighbors"]),
    }


REFERENCE_PARAMS = load_reference_params()
DEFAULT_REFERENCE_NEIGHBORS_ROOT = MODULE_DIR / "reference_neighbors"


def _load_sparse(data, name: str):
    shape = tuple(int(value) for value in data[f"{name}_shape"])
    return sparse.csr_matrix(
        (
            data[f"{name}_data"],
            data[f"{name}_indices"],
            data[f"{name}_indptr"],
        ),
        shape=shape,
    )


def load_reference_neighbors(
    dataset: str, output_root: Path = DEFAULT_REFERENCE_NEIGHBORS_ROOT
):
    output_root = Path(output_root)
    path = output_root / f"{dataset}.npz"
    if not path.is_file():
        raise FileNotFoundError(f"missing reference neighbor asset: {path}")
    with np.load(path, allow_pickle=False) as data:
        return {
            "path": path,
            "obs_names": np.asarray(data["obs_names"]).astype(str),
            "indices": np.asarray(data["indices"], dtype=np.int64),
            "distances": _load_sparse(data, "distances"),
            "connectivities": _load_sparse(data, "connectivities"),
            "metadata": json.loads(str(data["metadata"])),
        }


def apply_reference_neighbors(
    adata, dataset: str, output_root: Path = DEFAULT_REFERENCE_NEIGHBORS_ROOT
):
    output_root = Path(output_root)
    reference = load_reference_neighbors(dataset, output_root)
    current_names = np.asarray(adata.obs_names).astype(str)
    reference_names = reference["obs_names"]

    if len(current_names) != len(reference_names):
        return False, (
            f"cell count mismatch: method={len(current_names)}, "
            f"canonical={len(reference_names)}"
        )
    if len(np.unique(current_names)) != len(current_names):
        return False, "duplicate method cell identifiers"
    if set(current_names) != set(reference_names):
        return False, "cell identifier set mismatch"

    reference_position = {name: index for index, name in enumerate(reference_names)}
    row_order = np.asarray(
        [reference_position[name] for name in current_names], dtype=np.int64
    )
    current_position = np.empty(len(row_order), dtype=np.int64)
    current_position[row_order] = np.arange(len(row_order), dtype=np.int64)

    reference_indices = reference["indices"][row_order]
    mapped_indices = current_position[reference_indices]
    adata.obsp["distances"] = reference["distances"][row_order][:, row_order].tocsr()
    adata.obsp["connectivities"] = (
        reference["connectivities"][row_order][:, row_order].tocsr()
    )
    adata.uns["neighbors"] = {
        "indices": mapped_indices,
        "distances_key": "distances",
        "connectivities_key": "connectivities",
        "params": {
            "n_neighbors": REFERENCE_PARAMS["n_neighbors"],
            "n_pcs": REFERENCE_PARAMS["n_pcs"],
            "use_rep": REFERENCE_PARAMS["use_rep"],
            "metric": REFERENCE_PARAMS["metric"],
            "random_state": REFERENCE_PARAMS["random_state"],
            "reference_dataset": dataset,
        },
    }
    return True, (
        f"applied canonical graph: n_obs={len(current_names)}, "
        f"indices_shape={mapped_indices.shape}"
    )
