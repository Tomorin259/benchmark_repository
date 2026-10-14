"""The two frozen Forebrain coordinate frames shared by Cell and Ensemble.

PCA30 is the historical reference H5AD's X_pca, retained in the neighbor
asset. UMAP2 is the already evaluated BATC frame; neither frame is rebuilt.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_REFERENCE_EMBEDDINGS_ROOT = HERE / "reference_embeddings"
PCA_COORDINATE_SHA256 = "b1f98d84bd0f7fbd84d2c4e36bfb3b36d1e5b0bfd4dc0585e6867c0c470baaeb"
UMAP_COORDINATE_SHA256 = "0d8448c8e6cb2d7a635f773874257a1f40e57a17720df597138364f7491fec7d"


def array_sha256(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


def load_forebrain_frames(
    reference_neighbors_root: Path = None,
    reference_embeddings_root: Path = None,
) -> dict:
    """Validate the frozen bytes, dimensions, and common cell order."""
    neighbors = Path(reference_neighbors_root) if reference_neighbors_root else HERE / "reference_neighbors"
    embeddings = Path(reference_embeddings_root) if reference_embeddings_root else DEFAULT_REFERENCE_EMBEDDINGS_ROOT
    manifest = json.loads((embeddings / "manifest.json").read_text(encoding="utf-8"))
    entry = manifest["datasets"]["forebrain"]
    umap_path = embeddings / entry["umap_file"]
    if entry["umap_file"] != "forebrain_umap.npz":
        raise ValueError("Forebrain UMAP must use the module-local frozen asset")
    if hashlib.sha256(umap_path.read_bytes()).hexdigest() != entry["umap_file_sha256"]:
        raise ValueError("Forebrain UMAP asset SHA256 mismatch")
    with np.load(neighbors / "forebrain.npz", allow_pickle=False) as data:
        names = np.asarray(data["obs_names"]).astype(str)
        pca = np.asarray(data["reference_pca"]).copy()
    with np.load(umap_path, allow_pickle=False) as data:
        umap_names = np.asarray(data["obs_names"]).astype(str)
        umap = np.asarray(data["X_umap"]).copy()
    if len(names) != 1720 or len(np.unique(names)) != len(names):
        raise ValueError("Forebrain reference must contain 1720 unique cell identifiers")
    if not np.array_equal(names, umap_names):
        raise ValueError("Forebrain PCA and UMAP cell orders differ")
    if pca.shape != (1720, 30) or pca.dtype != np.dtype("float32"):
        raise ValueError("Forebrain PCA must be the frozen float32 1720x30 frame")
    if umap.shape != (1720, 2) or not np.isfinite(umap).all():
        raise ValueError("Forebrain UMAP must be a finite 1720x2 frame")
    pca_sha, umap_sha = array_sha256(pca), array_sha256(umap)
    if pca_sha != PCA_COORDINATE_SHA256 or entry["pca_coordinate_sha256"] != pca_sha:
        raise ValueError("Forebrain PCA coordinate SHA256 mismatch")
    if umap_sha != UMAP_COORDINATE_SHA256 or entry["umap_coordinate_sha256"] != umap_sha:
        raise ValueError("Forebrain UMAP coordinate SHA256 mismatch")
    metadata = {
        "cbdir_space": "pca", "pca_dimensions": 30, "batc_space": "umap",
        "umap_dimensions": 2, "pca_coordinate_sha256": pca_sha,
        "umap_coordinate_sha256": umap_sha,
        "pca_source": "reference_neighbors/forebrain.npz:reference_pca",
        "umap_source": "reference_embeddings/forebrain_umap.npz:X_umap",
    }
    return {"obs_names": names, "X_pca": pca, "X_umap": umap, "metadata": metadata}
