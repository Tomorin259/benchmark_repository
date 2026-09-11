#!/usr/bin/env python3
"""Reproduce the stream panels used by Figure 5A and Figure 5E."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


PLOTTING_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PLOTTING_ROOT / "stream_config.json"
VELOCITY_ROOT = PLOTTING_ROOT.parents[3]
DEFAULT_ARTIFACT_ROOT = VELOCITY_ROOT
DEFAULT_OUTPUTS = {
    "organoid": PLOTTING_ROOT.parent / "figures" / "fig5" / "figure5E_organoid",
    "rpe1_deepkinet": PLOTTING_ROOT.parent / "figures" / "fig5" / "figure5A_rpe1",
}


def load_config(path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def resolve(path, root):
    value = Path(path)
    if value.is_absolute() or ".." in value.parts:
        raise ValueError(f"stream paths must be relative to artifact root: {path}")
    return root / value


def validate_config(config, artifact_root, panel):
    specs = config.get(panel, [])
    if not specs:
        raise ValueError(f"no stream specs configured for panel {panel!r}")
    for spec in specs:
        for key in ("result", "asset", "reference"):
            value = spec.get(key)
            if value:
                path = resolve(value, artifact_root)
                if not path.is_file():
                    raise FileNotFoundError(path)
    return specs


def plot_organoid(specs, artifact_root, output, dpi):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc
    import scvelo as scv

    output.mkdir(parents=True, exist_ok=True)
    scv.settings.verbosity = 1
    for spec in specs:
        method = spec["method"]
        output_stem = spec["output_stem"]
        path = resolve(spec["result"], artifact_root)
        velocity_layer = spec["velocity_layer"]
        print(f"[{method}] {path}", flush=True)
        adata = sc.read_h5ad(path)
        if "clusters" not in adata.obs:
            for source in ("cell_type", "celltype"):
                if source in adata.obs:
                    adata.obs["clusters"] = adata.obs[source].copy()
                    break
        if "clusters" not in adata.obs:
            raise KeyError(f"{path}: no clusters/cell_type/celltype label")
        if velocity_layer not in adata.layers:
            raise KeyError(f"{path}: missing velocity layer {velocity_layer!r}")
        if method == "DeepKINET" and "velocity_graph" not in adata.uns:
            if "neighbors" not in adata.uns:
                scv.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
            scv.tl.velocity_graph(adata, vkey=velocity_layer, n_jobs=8)
        axis = scv.pl.velocity_embedding_stream(
            adata,
            basis="umap",
            color="clusters",
            vkey=velocity_layer,
            title="",
            show=False,
        )
        axis.figure.savefig(
            output / f"{output_stem}_organoid_stream_title_removed.png",
            dpi=dpi,
        )
        plt.close(axis.figure)


def plot_rpe1_deepkinet(specs, artifact_root, output, dpi):
    """Render the current DeepKINET RPE1 result in the Figure 5A style."""
    import matplotlib

    matplotlib.use("Agg")
    import scanpy as sc
    import scvelo as scv
    import numpy as np

    output.mkdir(parents=True, exist_ok=True)
    scv.settings.verbosity = 1
    for spec in specs:
        method = spec["method"]
        result_path = resolve(spec["result"], artifact_root)
        reference_path = resolve(spec["reference"], artifact_root)
        velocity_layer = spec["velocity_layer"]
        print(f"[{method}] {result_path}", flush=True)

        adata = sc.read_h5ad(result_path)
        reference = sc.read_h5ad(reference_path)
        if not np.array_equal(np.asarray(adata.obs_names), np.asarray(reference.obs_names)):
            reference = reference[adata.obs_names].copy()
        adata.obs["clusters"] = reference.obs["clusters"].copy()
        adata.uns["clusters_colors"] = np.asarray(reference.uns["clusters_colors"])
        adata.obsm["X_umap"] = np.asarray(reference.obsm["X_umap"]).copy()

        if velocity_layer not in adata.layers:
            raise KeyError(f"{result_path}: missing velocity layer {velocity_layer!r}")

        # The seeded result does not store scVelo's neighbor graph or embedding.
        # Recreate both from the saved gene-level velocity without modifying H5AD.
        adata.uns.pop("neighbors", None)
        for key in ("connectivities", "distances"):
            adata.obsp.pop(key, None)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30, use_highly_variable=True)
        scv.tl.velocity_graph(adata, vkey=velocity_layer, xkey="Ms", n_jobs=8)
        scv.tl.velocity_embedding(adata, basis="umap", vkey=velocity_layer)

        scv.pl.velocity_embedding_stream(
            adata,
            basis="umap",
            color="clusters",
            legend_fontsize=9,
            save=str(output / f"{spec['output_stem']}.png"),
            title="",
            vkey=velocity_layer,
            dpi=dpi,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, default=DEFAULT_ARTIFACT_ROOT)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--panel",
        choices=("organoid", "rpe1_manual", "rpe1_deepkinet"),
        default="organoid",
    )
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    artifact_root = args.artifact_root.expanduser().resolve()
    config = load_config(args.config.resolve())
    specs = validate_config(config, artifact_root, args.panel)
    if args.check_only:
        print(f"stream configuration valid: {args.panel}, {len(specs)} inputs")
        return
    if args.panel == "rpe1_manual":
        raise ValueError("rpe1_manual is a recorded manual panel and has no renderer")
    output = args.output or DEFAULT_OUTPUTS[args.panel]
    if args.panel == "organoid":
        plot_organoid(specs, artifact_root, output.resolve(), args.dpi)
    else:
        plot_rpe1_deepkinet(specs, artifact_root, output.resolve(), args.dpi)


if __name__ == "__main__":
    main()
