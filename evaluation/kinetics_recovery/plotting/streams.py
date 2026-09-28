#!/usr/bin/env python3
"""Render Figure 6 streams from validated projections in frozen result H5ADs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

PLOTTING_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PLOTTING_ROOT / "stream_config.json"
DEFAULT_ARTIFACT_ROOT = PLOTTING_ROOT.parents[3]
DEFAULT_OUTPUT_ROOT = PLOTTING_ROOT.parent / "figures" / "fig6"


def resolve(path, root):
    value = Path(path)
    if value.is_absolute() or ".." in value.parts:
        raise ValueError(f"stream paths must be relative to artifact root: {path}")
    resolved = root / value
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    return resolved


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, default=DEFAULT_ARTIFACT_ROOT)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--panel", choices=("organoid", "organoid_manual", "rpe1", "rpe1_manual", "rpe1_deepkinet"),
        default="organoid",
    )
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    root = args.artifact_root.expanduser().resolve()
    config = json.loads(args.config.read_text())
    panel = "rpe1" if args.panel == "rpe1_deepkinet" else args.panel
    specs = config[panel]
    if args.panel == "rpe1_deepkinet":
        specs = [spec for spec in specs if spec["method"] == "DeepKINET"]
    if panel.endswith("_manual"):
        for spec in specs:
            resolve(spec["asset"], root)
        if not args.check_only:
            raise ValueError("Dynamo is a recorded labelled-data manual panel; no renderer")
        print(f"{panel} manual asset exists; H5AD rendering provenance remains unverified")
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import scanpy as sc
    import scvelo as scv

    output = args.output or DEFAULT_OUTPUT_ROOT / (
        "figure6A_rpe1" if panel == "rpe1" else "figure6E_organoid"
    )
    if not args.check_only:
        if output.exists() and any(output.iterdir()):
            raise FileExistsError(f"Output directory must be empty: {output}")
        output.mkdir(parents=True, exist_ok=True)
    records = []
    for spec in specs:
        path = resolve(spec["result"], root)
        source_hash = sha256(path)
        adata = sc.read_h5ad(path)
        if adata.obs_names.has_duplicates:
            raise ValueError(f"{path}: duplicate cell IDs")
        vkey = spec["velocity_layer"]
        projection_key = f"{vkey}_umap"
        if vkey not in adata.layers or projection_key not in adata.obsm:
            raise KeyError(f"{path}: missing {vkey} or {projection_key}")
        X = np.asarray(adata.obsm["X_umap"])
        V = np.asarray(adata.obsm[projection_key])
        if X.shape != (adata.n_obs, 2) or V.shape != X.shape:
            raise ValueError(f"{path}: invalid embedding/projection shape")
        if not np.isfinite(X).all() or not np.isfinite(V).all():
            raise ValueError(f"{path}: non-finite embedding/projection values")
        record = {
            "method": spec["method"], "input": spec["result"], "input_sha256": source_hash,
            "n_cells": adata.n_obs, "n_genes": adata.n_vars,
            "embedding_key": "X_umap", "projection_key": projection_key,
            "projection_recomputed": False,
        }
        if spec.get("reference"):
            refpath = resolve(spec["reference"], root)
            reference = sc.read_h5ad(refpath)
            if not adata.obs_names.isin(reference.obs_names).all():
                raise ValueError(f"{refpath}: missing result cells")
            reference = reference[adata.obs_names]
            reference_X = np.asarray(reference.obsm["X_umap"])
            exact_match = np.array_equal(X, reference_X)
            dtype_match = (X.dtype != reference_X.dtype and
                           np.array_equal(X, reference_X.astype(X.dtype)))
            if not (exact_match or dtype_match):
                raise ValueError(f"{path}: saved UMAP differs from configured reference")
            label = "clusters" if "clusters" in reference.obs else "cell_type"
            adata.obs["clusters"] = reference.obs[label].copy()
            color_key = f"{label}_colors"
            if color_key in reference.uns:
                adata.uns["clusters_colors"] = np.asarray(reference.uns[color_key])
            record["display_reference"] = spec["reference"]
            record["display_reference_sha256"] = sha256(refpath)
            record["reference_embedding_exact_match"] = exact_match
            record["reference_embedding_equal_at_saved_dtype"] = bool(exact_match or dtype_match)
            record["reference_embedding_max_abs_difference"] = float(np.max(np.abs(X - reference_X)))
        elif "clusters" not in adata.obs:
            label = next((key for key in ("cell_type", "celltype") if key in adata.obs), None)
            if label is None:
                raise KeyError(f"{path}: no cluster/cell-type label")
            adata.obs["clusters"] = adata.obs[label].copy()
            if f"{label}_colors" in adata.uns:
                adata.uns["clusters_colors"] = adata.uns[f"{label}_colors"].copy()
        provenance = adata.uns.get(f"{vkey}_projection", {})
        record["projection_provenance_status"] = provenance.get(
            "status", "saved_H5AD_projection_without_independent_provenance_record"
        )
        print(f"Validated {spec['method']} {panel}: {projection_key}", flush=True)
        if not args.check_only:
            axis = scv.pl.velocity_embedding_stream(
                adata, basis="umap", vkey=vkey, X=X, V=V, recompute=False,
                color="clusters", title="", legend_fontsize=9 if panel == "rpe1" else None,
                show=False,
            )
            stem = spec["output_stem"]
            filename = f"{stem}.png" if panel == "rpe1" else f"{stem}_organoid_stream_title_removed.png"
            target = output / filename
            if panel == "rpe1":
                axis.figure.savefig(target, dpi=args.dpi, bbox_inches="tight", pad_inches=0.1)
            else:
                axis.figure.savefig(target, dpi=args.dpi)
            plt.close(axis.figure)
            record["output"] = filename
            record["output_sha256"] = sha256(target)
        if sha256(path) != source_hash:
            raise RuntimeError(f"Source H5AD changed during rendering: {path}")
        records.append(record)
    if not args.check_only:
        manifest = {
            "panel": panel, "projection_policy": "reuse_saved_X_and_V_without_recomputation",
            "versions": {"scanpy": sc.__version__, "scvelo": scv.__version__,
                         "numpy": np.__version__, "matplotlib": matplotlib.__version__},
            "dpi": args.dpi, "records": records,
        }
        if panel == "rpe1":
            manifest["savefig"] = {"bbox_inches": "tight", "pad_inches": 0.1}
        (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
