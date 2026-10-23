#!/usr/bin/env python3
"""Reproduce C07: conventional RNA-seq Dynamo on the 14 benchmark datasets.

C07 uses Dynamo's conventional stochastic model with unspliced/spliced RNA.
Inputs and preprocessing follow the dataset configuration below.
"""

from importlib.metadata import version
from pathlib import Path
import sys

# Direct execution must import installed packages, not same-named launchers.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)

import dynamo as dyn
import scanpy as sc
from dynamo.preprocessing import Preprocessor


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs"

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

# Fixed dataset inputs and annotation mappings. label_key
# supplies the standard clusters column. X_umap_ori is the embedding on which
# Dynamo projects velocity; its name is part of the C07 result contract.
DATASETS = {
    "bonemarrow": {
        "input": "bonemarrow/bonemarrow.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_tsne",
    },
    "capo_pre": {
        "input": "capo_pre/capo_pre.h5ad",
        "label_key": "phase",
        "embedding_key": "X_umap",
    },
    "DentateGyrus_lamanno": {
        "input": "DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_tsne",
    },
    "ForebrainGlut": {
        "input": "ForebrainGlut/ForebrainGlut.loom",
        "label_key": "Clusters",
        "embedding_key": None,
    },
    "Gastrulation": {
        "input": "Gastrulation/Gastrulation.h5ad",
        "label_key": "celltype",
        "embedding_key": "X_umap",
    },
    "Pancreas": {
        "input": "Pancreas/Pancreas.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_umap",
    },
    "dentategyrus": {
        "input": "dentategyrus/dentategyrus.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_umap",
    },
    "hematopoiesis_dynamo": {
        "input": "hematopoiesis_dynamo/hematopoiesis_dynamo.h5ad",
        "label_key": "cell_type",
        "embedding_key": "X_umap",
        "drop_layers": ("new", "total"),
    },
    "mousebrain": {
        "input": "mousebrain/mousebrain.h5ad",
        "label_key": "celltype",
        "embedding_key": "X_umap",
    },
    "neuron": {
        "input": "neuron/neuron.h5ad",
        "label_key": "time_str",
        "embedding_key": "X_umap",
    },
    "organoid": {
        "input": "organoid/organoid_pre.h5ad",
        "label_key": "cell_type",
        "embedding_key": "X_umap",
    },
    "retina": {
        "input": "retina/retina_full.h5ad",
        "label_key": "Annotation",
        "embedding_key": "X_umap",
    },
    "rpe1": {
        "input": "rpe1/rpe1.h5ad",
        "label_key": "cell_cycle_phase",
        "embedding_key": "X_RFP_GFP",
        "drop_layers": ("new", "total"),
        "drop_obs": (
            "Cell_cycle_possition",
            "Cell_cycle_relativePos",
            "exp_type",
            "time",
            "cell_cycle_phase",
        ),
    },
    "rpe1_gen": {
        "input": "rpe1_gen/rpe1_gen.h5ad",
        "label_key": "phase",
        "embedding_key": "X_umap",
    },
}


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / "dynamo"
    output_path = output_dir / f"dynamo_{dataset}.h5ad"

    print(f"\n=== Dynamo conventional RNA-seq: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    adata = sc.read(input_path, cache=False)

    missing_layers = {"spliced", "unspliced"} - set(adata.layers)
    if missing_layers:
        raise KeyError(f"{dataset}: missing required count layers {sorted(missing_layers)}")
    if config["label_key"] not in adata.obs:
        raise KeyError(f"{dataset}: missing label key {config['label_key']!r}")
    adata.obs["clusters"] = adata.obs[config["label_key"]].copy()

    for layer in config.get("drop_layers", ()):
        if layer not in adata.layers:
            raise KeyError(f"{dataset}: expected auxiliary layer {layer!r}")
        del adata.layers[layer]
    for key in config.get("drop_obs", ()):
        if key not in adata.obs:
            raise KeyError(f"{dataset}: expected auxiliary obs key {key!r}")
        del adata.obs[key]

    embedding_key = config["embedding_key"]
    if embedding_key is not None:
        if embedding_key not in adata.obsm:
            raise KeyError(f"{dataset}: missing embedding key {embedding_key!r}")
        adata.obsm["X_umap_ori"] = adata.obsm[embedding_key].copy()

    # This is the established C07 preprocessing and conventional stochastic
    # dynamics sequence. It creates a fresh visualization UMAP in X_umap while
    # retaining X_umap_ori for the benchmark velocity projection.
    preprocessor = Preprocessor()
    preprocessor.preprocess_adata(adata, recipe="monocle")
    dyn.tl.dynamics(adata)
    dyn.tl.moments(adata)
    dyn.tl.reduceDimension(adata, basis="umap", enforce=True)

    if "X_umap_ori" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            raise KeyError(f"{dataset}: Dynamo preprocessing did not provide X_pca")
        adata.obsm["X_umap_ori"] = adata.obsm["X_pca"].copy()
    dyn.tl.cell_velocities(adata, basis="umap_ori")
    dyn.tl.cell_wise_confidence(adata)

    if "velocity_S" not in adata.layers:
        raise RuntimeError(f"{dataset}: Dynamo did not create layers['velocity_S']")
    missing_obsm = {"X_umap", "X_umap_ori", "velocity_umap_ori"} - set(adata.obsm)
    if missing_obsm:
        raise RuntimeError(f"{dataset}: incomplete Dynamo result; missing {sorted(missing_obsm)}")

    output_dir.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    if version("dynamo-release") != "1.4.1":
        raise RuntimeError(
            f"Dynamo 1.4.1 is required, found {version('dynamo-release')}"
        )
    unknown_datasets = set(DATASETS_TO_RUN) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets={sorted(unknown_datasets)}")

    print(f"Dynamo version: {version('dynamo-release')}", flush=True)
    for dataset in DATASETS_TO_RUN:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
