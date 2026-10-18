#!/usr/bin/env python3
"""Reproduce C19 (GraphVelo on a scVelo dynamical-velocity input).

GraphVelo is applied after the visible scVelo dynamical workflow.  Its tangent
space model corrects the gene-space velocity and projects that correction to
the benchmark embedding.  This is distinct from C20, which uses Dynamo's
``M_s`` and ``velocity_S`` input instead.
"""

from pathlib import Path
import sys

# Direct execution must import installed packages, not same-named launchers.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)

import scanpy as sc
import scvelo as scv
from graphvelo.graph_velocity import GraphVelo
from graphvelo.utils import adj_to_knn


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs" / "graphvelo_scv"

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

N_JOBS = 8


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / dataset
    output_path = output_dir / f"graphvelo{dataset}.h5ad"

    print(f"\n=== GraphVelo / scVelo: {dataset} ===", flush=True)
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

    embedding_key = config["embedding_key"]
    if embedding_key is not None:
        if embedding_key not in adata.obsm:
            raise KeyError(f"{dataset}: missing embedding key {embedding_key!r}")
        adata.obsm["X_umap"] = adata.obsm[embedding_key].copy()

    if config["preprocessed"]:
        if adata.shape != config["expected_shape"]:
            raise ValueError(
                f"{dataset}: expected fixed benchmark shape {config['expected_shape']}, "
                f"found {adata.shape}"
            )
        missing_moments = {"Ms", "Mu"} - set(adata.layers)
        if missing_moments:
            raise KeyError(f"{dataset}: missing fixed preprocessing layers {sorted(missing_moments)}")
        print("using the fixed scVelo-preprocessed benchmark object", flush=True)
    else:
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            raise KeyError(f"{dataset}: preprocessing did not provide X_pca for the benchmark embedding")
        adata.obsm["X_umap"] = adata.obsm["X_pca"].copy()

    # scVelo supplies the initial gene-space velocity that GraphVelo corrects.
    scv.tl.recover_dynamics(adata, n_jobs=N_JOBS)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.latent_time(adata)

    output_dir.mkdir(parents=True, exist_ok=True)
    scv.settings.figdir = str(output_dir)
    scv.pl.velocity_embedding_stream(
        adata,
        color="clusters",
        basis="umap",
        save=f"PRE_plot_velocity_{dataset}.png",
    )

    if "connectivities" not in adata.obsp or "neighbors" not in adata.uns:
        raise KeyError(f"{dataset}: scVelo did not provide a neighbor graph for GraphVelo")
    indices, _ = adj_to_knn(adata.obsp["connectivities"])
    adata.uns["neighbors"]["indices"] = indices

    graphvelo = GraphVelo(adata, xkey="Ms", vkey="velocity")
    graphvelo.train()
    adata.layers["velocity_gv"] = graphvelo.project_velocity(adata.layers["Ms"])
    adata.obsm["gv_umap"] = graphvelo.project_velocity(adata.obsm["X_umap"])

    scv.pl.velocity_embedding_stream(
        adata,
        color="clusters",
        vkey="velocity_gv",
        basis="umap",
        X=adata.obsm["X_umap"],
        V=adata.obsm["gv_umap"],
        title="GraphVelo",
        figsize=(4, 3),
        save=f"graphvelo_plot_{dataset}.png",
    )

    required_layers = {"spliced", "unspliced", "Ms", "Mu", "velocity", "velocity_gv"}
    required_obsm = {"X_umap", "velocity_umap", "gv_umap"}
    missing_layers = required_layers - set(adata.layers)
    missing_obsm = required_obsm - set(adata.obsm)
    if missing_layers or missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete GraphVelo/scVelo result; "
            f"missing layers={sorted(missing_layers)}, obsm={sorted(missing_obsm)}, "
            f"velocity_graph={'velocity_graph' in adata.uns}"
        )

    adata.write_h5ad(output_path.as_posix())
    print(f"wrote: {output_path}", flush=True)


def main():
    if scv.__version__ != "0.2.4":
        raise RuntimeError(f"scVelo 0.2.4 is required, found {scv.__version__}")
    unknown_datasets = set(DATASETS_TO_RUN) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets: {sorted(unknown_datasets)}")

    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.set_figure_params("scvelo")
    print(f"scVelo version: {scv.__version__}", flush=True)
    for dataset in DATASETS_TO_RUN:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
