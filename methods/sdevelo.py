#!/usr/bin/env python3
"""Reproduce C17 (SDEvelo) with its direct SDENN workflow.

Raw inputs receive one scVelo ``filter_and_normalize -> moments`` pass.
The four fixed benchmark objects already contain the intended moments and are
passed to SDEvelo with ``args.process = False`` so they are never processed a
second time.  Select a physical GPU outside this script with
``CUDA_VISIBLE_DEVICES``; SDEvelo then uses device 0 within that visible set.
"""

from pathlib import Path
import sys

# Direct execution must import installed packages, not same-named launchers.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)

import scanpy as sc
import scvelo as scv
import sdevelo


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs" / "sdevelo"

# Edit this list for a focused rerun. SDEvelo's Config uses seed=0.
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

# These mappings are only the necessary data-format boundary for the 14
# benchmark inputs. Every input already provides spliced and unspliced RNA.
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
        "output_stem": "organoid_pre",
    },
    "retina": {
        "input": "retina/retina_full.h5ad",
        "label_key": "Annotation",
        "embedding_key": None,
        "preprocessed": False,
        "output_stem": "retina_full",
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


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / dataset
    output_stem = config.get("output_stem", dataset)
    output_path = output_dir / f"sdevelo_{output_stem}.h5ad"

    print(f"\n=== SDEvelo: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    if input_path.suffix.lower() == ".loom":
        adata = sc.read_loom(input_path, sparse=True)
    else:
        adata = sc.read_h5ad(input_path)

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
        print("using the fixed scVelo-preprocessed benchmark object", flush=True)
    else:
        # SDEvelo's supplied SDENN workflow starts from these scVelo moments.
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            raise KeyError(f"{dataset}: preprocessing did not provide X_pca for the benchmark embedding")
        adata.obsm["X_umap"] = adata.obsm["X_pca"].copy()

    # ``process=False`` is essential: preprocessing has already happened in
    # the visible block above.  The default Config keeps the upstream SDENN
    # seed (0), 300 epochs, mode=2, and its original SDE simulation.
    args = sdevelo.Config()
    args.process = False
    args.cuda_device = 0
    model = sdevelo.SDENN(args, adata)
    adata = model.train(args.nEpochs)
    adata.layers["velocity"] = adata.layers["sde_velocity"].copy()

    scv.tl.velocity_graph(adata, vkey="velocity")
    output_dir.mkdir(parents=True, exist_ok=True)
    scv.settings.figdir = str(output_dir)
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="clusters",
        vkey="velocity",
        save=f"_stream_{output_stem}.png",
        dpi=300,
        title="clusters",
    )

    required_layers = {"spliced", "unspliced", "Ms", "Mu", "sde_velocity", "velocity"}
    required_obsm = {"X_umap", "velocity_umap"}
    missing_layers = required_layers - set(adata.layers)
    missing_obsm = required_obsm - set(adata.obsm)
    if missing_layers or missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete SDEvelo result; missing layers={sorted(missing_layers)}, "
            f"obsm={sorted(missing_obsm)}, velocity_graph={'velocity_graph' in adata.uns}"
        )

    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    if sdevelo.__version__ != "0.2.12":
        raise RuntimeError(f"SDEvelo 0.2.12 is required, found {sdevelo.__version__}")
    unknown_datasets = set(DATASETS_TO_RUN) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets: {sorted(unknown_datasets)}")

    scv.settings.verbosity = 3
    scv.set_figure_params("scvelo", transparent=False)
    print(f"SDEvelo version: {sdevelo.__version__}", flush=True)
    print(f"scVelo version: {scv.__version__}", flush=True)
    print("global RNG seed: SDEvelo Config default (0)", flush=True)
    print("GPU selection: inherited from CUDA_VISIBLE_DEVICES", flush=True)

    for dataset in DATASETS_TO_RUN:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
