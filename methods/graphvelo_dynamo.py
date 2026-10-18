#!/usr/bin/env python3
"""Reproduce C20 (GraphVelo on a Dynamo velocity input).

The C20 AnnData output retains the layers and projected embeddings consumed
by benchmark evaluation.
"""

from pathlib import Path
import sys

# Direct execution must import installed packages, not same-named launchers.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)

import scanpy as sc
import dynamo as dyn
from graphvelo.graph_velocity import GraphVelo


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs" / "graphvelo_dynamo"

# Edit this list for a focused rerun.
DATASETS_TO_RUN = [
    "bonemarrow",
    "capo_pre",
    "DentateGyrus_lamanno",
    "ForebrainGlut",
    "Gastrulation",
    "Pancreas",
    "dentategyrus",
    "hematopoiesis_dynamo",
    "mousebrain",
    "neuron",
    "organoid",
    "retina",
    "rpe1",
    "rpe1_gen",
]

DATASETS = {
    "bonemarrow": {
        "input": "bonemarrow/bonemarrow.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_tsne",
    },
    "capo_pre": {"input": "capo_pre/capo_pre.h5ad", "label_key": "phase"},
    "DentateGyrus_lamanno": {
        "input": "DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_tsne",
    },
    "ForebrainGlut": {
        "input": "ForebrainGlut/ForebrainGlut.loom",
        "label_key": "Clusters",
    },
    "Gastrulation": {"input": "Gastrulation/Gastrulation.h5ad", "label_key": "celltype"},
    "Pancreas": {"input": "Pancreas/Pancreas.h5ad", "label_key": "clusters"},
    "dentategyrus": {"input": "dentategyrus/dentategyrus.h5ad", "label_key": "clusters"},
    "hematopoiesis_dynamo": {
        "input": "hematopoiesis_dynamo/hematopoiesis_dynamo.h5ad",
        "label_key": "cell_type",
        "remove_metabolic_layers": True,
    },
    "mousebrain": {"input": "mousebrain/mousebrain.h5ad", "label_key": "celltype"},
    "neuron": {"input": "neuron/neuron.h5ad", "label_key": "time_str"},
    "organoid": {"input": "organoid/organoid_pre.h5ad", "label_key": "cell_type"},
    "retina": {"input": "retina/retina_full.h5ad", "label_key": "Annotation"},
    "rpe1": {
        "input": "rpe1/rpe1.h5ad",
        "label_key": "cell_cycle_phase",
        "embedding_key": "X_RFP_GFP",
        "remove_metabolic_layers": True,
    },
    "rpe1_gen": {"input": "rpe1_gen/rpe1_gen.h5ad", "label_key": "phase"},
}


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / dataset
    output_path = output_dir / f"graphvelo_dynamo_{dataset}.h5ad"

    print(f"\n=== GraphVelo / Dynamo: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    if input_path.suffix.lower() == ".loom":
        adata = sc.read_loom(input_path.as_posix(), sparse=True)
    else:
        adata = sc.read_h5ad(input_path.as_posix())

    missing_layers = {"spliced", "unspliced"} - set(adata.layers)
    if missing_layers:
        raise KeyError(f"{dataset}: missing required count layers {sorted(missing_layers)}")
    if config["label_key"] not in adata.obs:
        raise KeyError(f"{dataset}: missing label key {config['label_key']!r}")
    adata.obs["clusters"] = adata.obs[config["label_key"]].copy()

    embedding_key = config.get("embedding_key")
    if embedding_key is not None:
        if embedding_key not in adata.obsm:
            raise KeyError(f"{dataset}: missing embedding key {embedding_key!r}")
        adata.obsm["X_umap"] = adata.obsm[embedding_key].copy()
    elif "X_umap" not in adata.obsm:
        for candidate in ("X_cdr", "X_tsne", "X_pca"):
            if candidate in adata.obsm:
                adata.obsm["X_umap"] = adata.obsm[candidate].copy()
                break

    # Dynamo's conventional recipe is used for every C20 input. The two
    # metabolic-label inputs must first drop new/total and time metadata, as
    # in the archived C20 run, so dynamics uses the spliced/unspliced model.
    if config.get("remove_metabolic_layers", False):
        for layer in ("new", "total"):
            if layer in adata.layers:
                del adata.layers[layer]
        if "time" in adata.obs:
            del adata.obs["time"]

    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata)
    dyn.tl.neighbors(adata)

    if "M_s" not in adata.layers or "velocity_S" not in adata.layers:
        raise KeyError(f"{dataset}: Dynamo did not produce M_s and velocity_S")
    graphvelo = GraphVelo(adata, xkey="M_s", vkey="velocity_S")
    graphvelo.train()
    adata.layers["velocity_gv"] = graphvelo.project_velocity(adata.layers["M_s"])
    if "X_pca" in adata.obsm:
        adata.obsm["gv_pca"] = graphvelo.project_velocity(adata.obsm["X_pca"])
    if "X_umap" in adata.obsm:
        adata.obsm["gv_umap"] = graphvelo.project_velocity(adata.obsm["X_umap"])
    else:
        raise KeyError(f"{dataset}: no embedding available for GraphVelo projection")

    output_dir.mkdir(parents=True, exist_ok=True)
    dyn.pl.streamline_plot(
        adata,
        color=["clusters"],
        vector="gv",
        save_show_or_return="save",
        save_kwargs={
            # Dynamo's ``path`` includes the file stem rather than a directory.
            # Giving it an explicit stem keeps the PNG beside this dataset's H5AD.
            "path": str(output_dir / "graphvelo_dynamo"),
            "prefix": dataset,
            "ext": "png",
        },
    )

    # Write exactly the compact C20 interchange contract rather than trying
    # to serialize Dynamo's internal kinetic objects.  It retains every field
    # consumed by evaluation and by the GraphVelo visualizations.
    keep_obs = [
        key
        for key in (
            "clusters",
            "celltype",
            "cell_type",
            "Clusters",
            "phase",
            "Annotation",
            "time_str",
            "cell_cycle_phase",
            "true_branch",
            "palantir_pseudotime",
        )
        if key in adata.obs
    ]
    result = sc.AnnData(
        X=None,
        obs=adata.obs[keep_obs].copy(),
        var=adata.var.iloc[:, :0].copy(),
    )
    result.obs_names = adata.obs_names.copy()
    result.var_names = adata.var_names.copy()
    for key in ("X_umap", "X_pca", "X_tsne", "gv_umap", "gv_pca"):
        if key in adata.obsm:
            result.obsm[key] = adata.obsm[key].copy()
    for key in ("spliced", "unspliced", "velocity_S", "velocity_gv"):
        if key in adata.layers:
            result.layers[key] = adata.layers[key].copy()
    if "clusters_colors" in adata.uns:
        result.uns["clusters_colors"] = adata.uns["clusters_colors"]

    required_layers = {"spliced", "unspliced", "velocity_S", "velocity_gv"}
    required_obsm = {"X_umap", "gv_umap"}
    missing_layers = required_layers - set(result.layers)
    missing_obsm = required_obsm - set(result.obsm)
    if missing_layers or missing_obsm:
        raise RuntimeError(
            f"{dataset}: incomplete GraphVelo/Dynamo result; "
            f"missing layers={sorted(missing_layers)}, obsm={sorted(missing_obsm)}"
        )

    result.write_h5ad(output_path.as_posix())
    print(f"wrote: {output_path}", flush=True)


def main():
    unknown_datasets = set(DATASETS_TO_RUN) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets: {sorted(unknown_datasets)}")
    print(f"Dynamo version: {dyn.__version__}", flush=True)
    for dataset in DATASETS_TO_RUN:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
