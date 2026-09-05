#!/usr/bin/env python3
"""Reproduce C08: TIVelo on the 14 benchmark datasets.

TIVelo's official entry point requires a cluster label and a two-dimensional
embedding. The mappings below make those dataset-specific inputs explicit;
the remaining sequence is the original scVelo preprocessing followed by the
published TIVelo call.
"""

from importlib.metadata import version
from pathlib import Path

import scanpy as sc
import scvelo as scv
import tivelo


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

# These mappings are part of the benchmark data contract. TIVelo uses
# clusters only for path inference and visualization; it still estimates a
# cell-level velocity from spliced and unspliced counts. The four entries
# marked preprocessed already contain scVelo moments and their established
# neighbor graph, so they must enter TIVelo unchanged.
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
        "preprocessed_shape": (1222, 2000),
    },
    "DentateGyrus_lamanno": {
        "input": "DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_tsne",
    },
    "ForebrainGlut": {
        "input": "ForebrainGlut/ForebrainGlut.loom",
        "label_key": "Clusters",
        "embedding_key": "X_pca",
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
    },
    "mousebrain": {
        "input": "mousebrain/mousebrain.h5ad",
        "label_key": "celltype",
        "embedding_key": "X_umap",
        "preprocessed_shape": (3365, 936),
    },
    "neuron": {
        "input": "neuron/neuron.h5ad",
        "label_key": "time_str",
        "embedding_key": "X_umap",
        "preprocessed_shape": (3066, 2000),
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
    },
    "rpe1_gen": {
        "input": "rpe1_gen/rpe1_gen.h5ad",
        "label_key": "phase",
        "embedding_key": "X_umap",
        "preprocessed_shape": (3354, 2000),
    },
}


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / "tivelo"
    output_path = output_dir / f"tivelo_{dataset}.h5ad"
    tivelo_cache = output_dir / dataset / "tivelo.h5ad"

    print(f"\n=== TIVelo: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    adata = sc.read(input_path, cache=False)
    adata.var_names_make_unique()

    missing_layers = {"spliced", "unspliced"} - set(adata.layers)
    if missing_layers:
        raise KeyError(f"{dataset}: missing required count layers {sorted(missing_layers)}")
    label_key = config["label_key"]
    if label_key not in adata.obs:
        raise KeyError(f"{dataset}: missing label key {label_key!r}")
    embedding_key = config["embedding_key"]
    if embedding_key not in adata.obsm:
        raise KeyError(f"{dataset}: missing embedding key {embedding_key!r}")

    adata.obs["clusters"] = adata.obs[label_key].astype(str).astype("category")
    adata.obsm["X_umap"] = adata.obsm[embedding_key].copy()

    if "preprocessed_shape" in config:
        expected_shape = config["preprocessed_shape"]
        if adata.shape != expected_shape:
            raise ValueError(
                f"{dataset}: expected established preprocessed shape "
                f"{expected_shape}, found {adata.shape}"
            )
        missing_moments = {"Ms", "Mu"} - set(adata.layers)
        if missing_moments:
            raise KeyError(
                f"{dataset}: missing established moments {sorted(missing_moments)}"
            )
        if "neighbors" not in adata.uns:
            raise KeyError(f"{dataset}: missing established scVelo neighbors")
        print("using established scVelo moments and neighbors", flush=True)
    else:
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    output_dir.mkdir(parents=True, exist_ok=True)
    # TIVelo silently reuses this file if present. Removing only this
    # method-owned cache makes an explicit rerun execute the model again.
    if tivelo_cache.exists():
        if not tivelo_cache.is_file():
            raise RuntimeError(f"{dataset}: unexpected TIVelo cache path {tivelo_cache}")
        tivelo_cache.unlink()

    adata = tivelo.tivelo(
        adata,
        group_key="clusters",
        emb_key="X_umap",
        data_name=dataset,
        save_folder=str(output_dir),
        njobs=-1,
        show_fig=False,
        filter_genes=True,
        constrain=True,
        loss_fun="mse",
        only_s=False,
        alpha_1=1,
        alpha_2=0.1,
        batch_size=1024,
        n_epochs=100,
        velocity_key="velocity",
        measure_performance=False,
        cluster_edges=None,
    )

    if "velocity" not in adata.layers and "velo_s_norm" in adata.layers:
        adata.layers["velocity"] = adata.layers["velo_s_norm"].copy()
    if "velocity_umap" not in adata.obsm and "velo_s_norm_umap" in adata.obsm:
        adata.obsm["velocity_umap"] = adata.obsm["velo_s_norm_umap"].copy()

    missing_layers = {"velocity"} - set(adata.layers)
    missing_obsm = {"X_umap", "velocity_umap"} - set(adata.obsm)
    if missing_layers or missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete TIVelo result; "
            f"layers={sorted(missing_layers)}, obsm={sorted(missing_obsm)}, "
            f"velocity_graph={'velocity_graph' in adata.uns}"
        )

    adata.write_h5ad(output_path)
    # Keep TIVelo's figures but remove its duplicate cache H5AD.
    if tivelo_cache.exists():
        tivelo_cache.unlink()
    print(f"wrote: {output_path}", flush=True)


def main():
    expected_versions = {
        "tivelo": "0.1.4",
        "scvelo": "0.3.1",
        "scanpy": "1.10.3",
    }
    actual_versions = {name: version(name) for name in expected_versions}
    if actual_versions != expected_versions:
        raise RuntimeError(
            f"TIVelo environment mismatch: expected={expected_versions}, "
            f"found={actual_versions}"
        )

    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.set_figure_params("scvelo")

    unknown_datasets = set(DATASETS_TO_RUN) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets={sorted(unknown_datasets)}")
    for dataset in DATASETS_TO_RUN:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
