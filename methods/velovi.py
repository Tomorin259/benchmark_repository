#!/usr/bin/env python3
"""Reproduce C12 (veloVI) with the upstream tutorial workflow.

The April benchmark run used velovi 0.3.1, 2,000 HVGs for raw inputs,
1,000 maximum epochs, a batch size of 64, and 25 posterior samples.  These
explicit benchmark settings are retained below; the calculation order follows
the official tutorial: scVelo preprocessing, ``preprocess_data``, model setup,
training, posterior velocity, graph, and stream plot.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvelo as scv
import torch
import velovi
from velovi import VELOVI, preprocess_data


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs" / "velovi"

# Edit this list for a focused rerun. veloVI's tutorial does not prescribe an
# extra global seed, so none is set here.
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

# These mappings are the fixed data boundary of this benchmark, not a generic
# adapter. ``embedding_key=None`` means that preprocessing creates X_pca, which
# is then used as the stored embedding (the legacy Forebrain convention).
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
        "expected_shape": (1222, 2000),
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
        "expected_shape": (3365, 936),
    },
    "neuron": {
        "input": "neuron/neuron.h5ad",
        "label_key": "time_str",
        "embedding_key": None,
        "preprocessed": True,
        "expected_shape": (3066, 2000),
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
        "expected_shape": (3354, 2000),
    },
}

# Retain the published benchmark settings. The upstream tutorial's bare
# ``vae.train()`` call uses its package defaults instead.
MAX_EPOCHS = 1000
BATCH_SIZE = 64
N_POSTERIOR_SAMPLES = 25
N_JOBS = 8


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / dataset
    output_path = output_dir / f"veloVI_{dataset}_all.h5ad"

    print(f"\n=== veloVI: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    adata = sc.read(input_path, cache=False)
    adata.var_names_make_unique()

    missing_layers = {"spliced", "unspliced"} - set(adata.layers)
    if missing_layers:
        raise KeyError(f"{dataset}: missing count layers {sorted(missing_layers)}")
    if config["label_key"] not in adata.obs:
        raise KeyError(f"{dataset}: missing label key {config['label_key']!r}")
    adata.obs["clusters"] = adata.obs[config["label_key"]].copy()

    embedding_key = config["embedding_key"]
    if embedding_key is not None:
        if embedding_key not in adata.obsm:
            raise KeyError(f"{dataset}: missing embedding key {embedding_key!r}")
        adata.obsm["X_umap"] = adata.obsm[embedding_key].copy()

    if config["preprocessed"]:
        if adata.shape != config["expected_shape"]:
            raise ValueError(
                f"{dataset}: expected fixed benchmark shape "
                f"{config['expected_shape']}, found {adata.shape}"
            )
        missing_moments = {"Ms", "Mu"} - set(adata.layers)
        if missing_moments:
            raise KeyError(f"{dataset}: missing fixed moments {sorted(missing_moments)}")
        print("using the fixed scVelo-preprocessed benchmark object", flush=True)
    else:
        # The finalized benchmark uses the same 20-count / 2,000-HVG settings
        # as its historical veloVI run.
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            raise KeyError(f"{dataset}: preprocessing did not create X_pca")
        adata.obsm["X_umap"] = adata.obsm["X_pca"].copy()

    # Upstream veloVI tutorial: initialize moment-derived kinetic quantities,
    # register Ms/Mu, train the model, then draw posterior velocities.
    adata = preprocess_data(adata)
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train(max_epochs=MAX_EPOCHS, batch_size=BATCH_SIZE)

    adata.layers["velocity"] = np.asarray(
        vae.get_velocity(n_samples=N_POSTERIOR_SAMPLES, velo_statistic="mean")
    )
    adata.layers["velocity_u"] = np.asarray(
        vae.get_velocity(
            n_samples=N_POSTERIOR_SAMPLES,
            velo_statistic="mean",
            velo_mode="unspliced",
        )
    )

    scv.tl.velocity_graph(adata, vkey="velocity", n_jobs=N_JOBS)
    output_dir.mkdir(parents=True, exist_ok=True)
    # scVelo 0.3.1 expects a string here when it saves the tutorial STREAM plot.
    scv.settings.figdir = str(output_dir)
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="clusters",
        show=False,
        save=f"_{dataset}_stream_all.png",
    )
    plt.close("all")

    missing_layers = {"spliced", "unspliced", "Ms", "Mu", "velocity", "velocity_u"} - set(
        adata.layers
    )
    missing_obsm = {"X_umap", "velocity_umap"} - set(adata.obsm)
    if missing_layers or missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete veloVI result; missing layers={sorted(missing_layers)}, "
            f"obsm={sorted(missing_obsm)}, velocity_graph={'velocity_graph' in adata.uns}"
        )

    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    if velovi.__version__ != "0.3.1":
        raise RuntimeError(f"velovi 0.3.1 is required, found {velovi.__version__}")
    unknown_datasets = set(DATASETS_TO_RUN) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets: {sorted(unknown_datasets)}")

    torch.set_float32_matmul_precision("medium")
    if torch.cuda.is_available():
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.set_figure_params("scvelo")
    print(f"velovi version: {velovi.__version__}", flush=True)
    print(f"device: {'cuda' if torch.cuda.is_available() else 'cpu'}", flush=True)

    for dataset in DATASETS_TO_RUN:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
