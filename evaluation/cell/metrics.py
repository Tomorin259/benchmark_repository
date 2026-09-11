#!/usr/bin/env python3
"""Pure helpers shared by the canonical cell-level evaluation workers."""

from __future__ import annotations

import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MODULE_DIR / "config.json"


def load_config(path: Path = DEFAULT_CONFIG) -> dict:
    """Read the dataset and metric contract used by the canonical workflow."""
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def resolve_artifact(artifact_root: Path, relative_path: str) -> Path:
    """Resolve a file below an explicit artifact root."""
    path = Path(relative_path)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"result paths must be relative to artifact root: {relative_path}")
    candidate = artifact_root / path
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    return candidate.resolve()


def dense(values):
    if sparse.issparse(values):
        return values.toarray()
    return np.asarray(values)


def prepare_labels(adata: ad.AnnData, dataset: str, config: dict) -> str:
    dataset_config = config["datasets"][dataset]
    target = dataset_config["label"]
    label_config = config["labels"]
    source = label_config.get("dataset_specific_sources", {}).get(dataset)
    if source and source in adata.obs:
        adata.obs[target] = adata.obs[source].copy()
    elif target not in adata.obs:
        for candidate in label_config.get("fallback_sources", []):
            if candidate in adata.obs:
                adata.obs[target] = adata.obs[candidate].copy()
                break
    if target not in adata.obs:
        raise KeyError(f"{dataset}: label column {target!r} is missing")
    return target


def consistency_scores(
    velocity: np.ndarray, neighbor_indices, n_neighbors: int
) -> np.ndarray:
    """Return each cell's mean cosine similarity to its nearest neighbors."""
    velocity = np.asarray(velocity, dtype=np.float64)
    scores = np.zeros(velocity.shape[0], dtype=float)
    for index in range(velocity.shape[0]):
        # Reference neighbor assets store the source cell in slot zero.
        neighbors = np.asarray(neighbor_indices[index])[1 : n_neighbors + 1]
        if neighbors.size == 0:
            continue
        current = velocity[index]
        neighbor_values = velocity[neighbors]
        valid = np.isfinite(current) & np.isfinite(neighbor_values).all(axis=0)
        if not valid.any():
            continue
        current = current[valid]
        neighbor_values = neighbor_values[:, valid]
        current_norm = np.linalg.norm(current)
        neighbor_norms = np.linalg.norm(neighbor_values, axis=1)
        denominators = current_norm * neighbor_norms
        valid_norms = denominators > 0
        if valid_norms.any():
            scores[index] = np.mean(
                neighbor_values[valid_norms].dot(current)
                / denominators[valid_norms]
            )
    return scores


def parse_unitvelo_means(output: str) -> tuple[float, float]:
    cbdir = re.search(
        r"# Cross-Boundary Direction Correctness \(A->B\).*?Total Mean:\s*([^\n\r]+)",
        output,
        flags=re.DOTALL,
    )
    iccoh = re.search(
        r"# In-cluster Coherence.*?Total Mean:\s*([^\n\r]+)",
        output,
        flags=re.DOTALL,
    )
    if cbdir is None or iccoh is None:
        raise ValueError("unitvelo.evaluate output did not contain both Total Mean values")
    return float(cbdir.group(1).strip()), float(iccoh.group(1).strip())
