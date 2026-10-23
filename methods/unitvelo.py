#!/usr/bin/env python3
"""Reproduce C05 and C06 with the UniTVelo tutorial workflow.

Inputs and preprocessing follow the dataset configuration below.
"""

import argparse
from pathlib import Path
import sys

import os
# Direct execution must import installed packages, not same-named launchers.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)

import pkg_resources
import scanpy as sc
import scvelo as scv
import unitvelo as utv


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))
from methods.simulation_inputs import simulation_inputs
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs"

# Edit these lists for a focused rerun. FIT_OPTION 1 is UniTVelo's
# unified-time model (C05); FIT_OPTION 2 is its independent-time model (C06).
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
FIT_OPTIONS_TO_RUN = ["1", "2"]
OUTPUT_DIRECTORY = {"1": "unitvelo", "2": "unitvelo_mode2"}
GPU_INDEX = 0

# Dataset inputs, preprocessing settings and annotation mappings.
# label_key supplies clusters; Forebrain uses the PCA coordinates.
DATASETS = {
    "bonemarrow": {
        "input": "bonemarrow/bonemarrow.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_tsne",
        "preprocessed": False,
    },
    "capo_pre": {
        "input": "capo_pre/capo_pre.h5ad",
        "label_key": "phase",
        "embedding_key": None,
        "preprocessed": True,
    },
    "DentateGyrus_lamanno": {
        "input": "DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_tsne",
        "preprocessed": False,
    },
    "ForebrainGlut": {
        "input": "ForebrainGlut/ForebrainGlut.loom",
        "label_key": "Clusters",
        "embedding_key": None,
        "preprocessed": False,
    },
    "Gastrulation": {
        "input": "Gastrulation/Gastrulation.h5ad",
        "label_key": "celltype",
        "embedding_key": None,
        "preprocessed": False,
    },
    "Pancreas": {
        "input": "Pancreas/Pancreas.h5ad",
        "label_key": "clusters",
        "embedding_key": None,
        "preprocessed": False,
    },
    "dentategyrus": {
        "input": "dentategyrus/dentategyrus.h5ad",
        "label_key": "clusters",
        "embedding_key": None,
        "preprocessed": False,
    },
    "hematopoiesis_dynamo": {
        "input": "hematopoiesis_dynamo/hematopoiesis_dynamo.h5ad",
        "label_key": "cell_type",
        "embedding_key": None,
        "preprocessed": False,
    },
    "mousebrain": {
        "input": "mousebrain/mousebrain.h5ad",
        "label_key": "celltype",
        "embedding_key": None,
        "preprocessed": True,
    },
    "neuron": {
        "input": "neuron/neuron.h5ad",
        "label_key": "time_str",
        "embedding_key": None,
        "preprocessed": True,
    },
    "organoid": {
        "input": "organoid/organoid_pre.h5ad",
        "label_key": "cell_type",
        "embedding_key": None,
        "preprocessed": False,
    },
    "retina": {
        "input": "retina/retina_full.h5ad",
        "label_key": "Annotation",
        "embedding_key": None,
        "preprocessed": False,
    },
    "rpe1": {
        "input": "rpe1/rpe1.h5ad",
        "label_key": "cell_cycle_phase",
        "embedding_key": "X_RFP_GFP",
        "preprocessed": False,
    },
    "rpe1_gen": {
        "input": "rpe1_gen/rpe1_gen.h5ad",
        "label_key": "phase",
        "embedding_key": None,
        "preprocessed": True,
    },
}

DATASETS.update(simulation_inputs())


def run_dataset(dataset, fit_option, data_root, output_root):
    config = DATASETS[dataset]
    input_path = data_root / config["input"]
    output_dir = output_root / OUTPUT_DIRECTORY[fit_option]
    output_path = output_dir / f"unitvelo_{dataset}.h5ad"
    if output_path.exists():
        raise FileExistsError(f"output exists: {output_path}; choose a new output root")

    print(f"\n=== UniTVelo FIT_OPTION={fit_option}: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    adata = sc.read(input_path, cache=False)

    missing_layers = {"spliced", "unspliced"} - set(adata.layers)
    if missing_layers:
        raise KeyError(f"{dataset}: missing required count layers {sorted(missing_layers)}")
    if config["label_key"] not in adata.obs:
        raise KeyError(f"{dataset}: missing label key {config['label_key']!r}")
    adata.obs["clusters"] = adata.obs[config["label_key"]].copy()

    embedding_key = config["embedding_key"]
    if embedding_key is not None:
        if embedding_key not in adata.obsm:
            raise KeyError(f"{dataset}: missing embedding key {embedding_key!r}")
        adata.obsm["X_umap"] = adata.obsm[embedding_key].copy()

    if config["preprocessed"]:
        missing_moments = {"Ms", "Mu"} - set(adata.layers)
        if missing_moments:
            raise KeyError(f"{dataset}: missing fixed preprocessing layers {sorted(missing_moments)}")
        print("input moments: Ms, Mu", flush=True)

    # UniTVelo writes its own temporary res directory relative to the current
    # directory when it receives an AnnData object. Keep that package side
    # effect inside this dataset's ignored output work directory.
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / f"{dataset}_unitvelo_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    # UniTVelo's tutorial uses Configuration and run_model directly. All
    # remaining values are the package defaults; only the benchmarked mode and
    # visible GPU selection differ between C05 and C06.
    velo_config = utv.config.Configuration()
    velo_config.FIT_OPTION = fit_option
    velo_config.GPU = GPU_INDEX
    if config.get("simulation"):
        velo_config.R2_ADJUST = True
        velo_config.IROOT = None
        velo_config.AGENES_R2 = 1
    previous_cwd = Path.cwd()
    os.chdir(work_dir)
    try:
        adata = utv.run_model(
            adata,
            "clusters",
            config_file=velo_config,
            normalize=not config["preprocessed"],
        )
    finally:
        os.chdir(previous_cwd)

    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            raise KeyError(f"{dataset}: UniTVelo preprocessing did not provide X_pca")
        adata.obsm["X_umap"] = adata.obsm["X_pca"].copy()
    if "velocity" not in adata.layers:
        raise RuntimeError(f"{dataset}: UniTVelo did not create layers['velocity']")
    if "velocity_graph" not in adata.uns:
        scv.tl.velocity_graph(adata)
    # Compute the saved velocity embedding explicitly rather than through the
    # stream plotting call shown in the tutorial.
    scv.tl.velocity_embedding(adata, basis="umap")

    missing_obsm = {"X_umap", "velocity_umap"} - set(adata.obsm)
    if missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete UniTVelo result; missing obsm={sorted(missing_obsm)}, "
            f"velocity_graph={'velocity_graph' in adata.uns}"
        )

    adata.uns["benchmark_run"] = {
        "operational_id": OUTPUT_DIRECTORY[fit_option], "FIT_OPTION": fit_option,
        "dataset": dataset, "unitvelo_version": pkg_resources.get_distribution("unitvelo").version,
    }
    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", dest="datasets", action="append", choices=sorted(DATASETS))
    parser.add_argument("--fit-option", dest="fit_options", action="append", choices=sorted(OUTPUT_DIRECTORY))
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    args = parser.parse_args()
    datasets = args.datasets or DATASETS_TO_RUN
    fit_options = args.fit_options or FIT_OPTIONS_TO_RUN
    if any(DATASETS[name].get("simulation") for name in datasets) and fit_options != ["1"]:
        raise ValueError("Simulation benchmark uses --fit-option 1")
    version = pkg_resources.get_distribution("unitvelo").version
    if version != "0.2.5.2":
        raise RuntimeError(f"UniTVelo 0.2.5.2 is required, found {version}")
    unknown_datasets = set(datasets) - set(DATASETS)
    unknown_options = set(fit_options) - set(OUTPUT_DIRECTORY)
    if unknown_datasets or unknown_options:
        raise ValueError(
            f"unknown datasets={sorted(unknown_datasets)}; "
            f"unknown FIT_OPTION values={sorted(unknown_options)}"
        )

    scv.settings.verbosity = 0
    print(f"UniTVelo version: {version}", flush=True)
    print(f"scVelo version: {scv.__version__}", flush=True)
    for fit_option in fit_options:
        for dataset in datasets:
            run_dataset(dataset, fit_option, args.data_root.expanduser().resolve(), args.output_root.expanduser().resolve())


if __name__ == "__main__":
    main()
