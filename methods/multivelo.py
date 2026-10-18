#!/usr/bin/env python3
"""Reproduce C04: RNA-only MultiVelo on the 14 benchmark datasets.

Set DATA_ROOT to the archived input artifacts before running. The four
already-preprocessed benchmark objects are intentionally not preprocessed a
second time; every other dataset follows the scVelo preprocessing used before
the RNA-only MultiVelo fit.
"""

from pathlib import Path
import sys

# Direct execution must import installed packages, not same-named launchers.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)

import multivelo as mv
import scanpy as sc
import scvelo as scv


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs"

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

# Fixed dataset inputs and annotation mappings. All 14 inputs
# already expose spliced and unspliced. label_key supplies the standard
# clusters column saved with the method result. Forebrain receives X_pca after
# moments because its source loom has no precomputed embedding.
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
MULTIVELO_N_JOBS = 12


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / "multivelo"
    output_path = output_dir / f"multivelo_{dataset}.h5ad"

    print(f"\n=== MultiVelo RNA-only: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    # scVelo 0.3.3 delegates file loading to Scanpy.
    adata = sc.read(input_path, cache=False)
    adata.var_names_make_unique()

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
                f"{dataset}: expected fixed benchmark shape "
                f"{config['expected_shape']}, found {adata.shape}"
            )
        missing_moments = {"Ms", "Mu"} - set(adata.layers)
        if missing_moments:
            raise KeyError(f"{dataset}: missing fixed preprocessing layers {sorted(missing_moments)}")
        print("using the fixed scVelo-preprocessed benchmark object", flush=True)
    else:
        # The scVelo preprocessing preceding the RNA-only MultiVelo fit.
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            raise KeyError(f"{dataset}: preprocessing did not provide X_pca for the benchmark embedding")
        adata.obsm["X_umap"] = adata.obsm["X_pca"].copy()

    adata = mv.recover_dynamics_chrom(
        adata,
        adata_atac=None,
        max_iter=5,
        init_mode="invert",
        parallel=True,
        n_jobs=MULTIVELO_N_JOBS,
        save_plot=False,
        rna_only=True,
        fit=True,
        n_anchors=500,
        extra_color_key=None,
        verbose=False,
    )
    scv.tl.velocity_graph(adata, n_jobs=N_JOBS)
    # Compute the saved velocity embedding explicitly rather than through a
    # plotting side effect.
    scv.tl.velocity_embedding(adata, basis="umap")

    missing_layers = {"spliced", "unspliced", "velocity"} - set(adata.layers)
    missing_obsm = {"X_umap", "velocity_umap"} - set(adata.obsm)
    if missing_layers or missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete MultiVelo result; missing layers={sorted(missing_layers)}, "
            f"obsm={sorted(missing_obsm)}, velocity_graph={'velocity_graph' in adata.uns}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    if mv.__version__ != "0.1.0":
        raise RuntimeError(f"MultiVelo 0.1.0 is required, found {mv.__version__}")
    if scv.__version__ != "0.3.3":
        raise RuntimeError(f"scVelo 0.3.3 is required, found {scv.__version__}")
    unknown_datasets = set(DATASETS_TO_RUN) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets={sorted(unknown_datasets)}")

    # These tutorial settings also affect the automatic scale stored by
    # velocity_embedding, so they are part of the numerical reproduction.
    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.set_figure_params("scvelo")
    print(f"MultiVelo version: {mv.__version__}", flush=True)
    print(f"scVelo version: {scv.__version__}", flush=True)
    for dataset in DATASETS_TO_RUN:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
