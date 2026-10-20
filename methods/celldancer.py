#!/usr/bin/env python3
"""Reproduce C15 (cellDancer) with the published dataframe workflow.

The visible stages are: scVelo preprocessing, conversion of Ms/Mu and a
two-dimensional embedding to cellDancer's dataframe, cellDancer fitting,
transfer of its kinetics back to the benchmark AnnData, and scVelo's velocity
graph/embedding projection.
"""

import argparse
from pathlib import Path
import sys

# Direct execution must import installed packages, not same-named launchers.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvelo as scv
from scipy import sparse

import celldancer as cd


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))
from methods.simulation_inputs import simulation_inputs
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs" / "celldancer"

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

# ``embed_key`` is the two-dimensional representation passed to cellDancer.
# Forebrain deliberately uses the first two PCs: its native input has no UMAP,
# and its published cellDancer result consequently carries velocity_pca.
DATASETS = {
    "bonemarrow": {
        "input": "bonemarrow/bonemarrow.h5ad",
        "label_key": "clusters",
        "source_embedding": "X_tsne",
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": False,
    },
    "capo_pre": {
        "input": "capo_pre/capo_pre.h5ad",
        "label_key": "phase",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": True,
        "expected_shape": (1222, 2000),
    },
    "DentateGyrus_lamanno": {
        "input": "DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad",
        "label_key": "clusters",
        "source_embedding": "X_tsne",
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": False,
    },
    "ForebrainGlut": {
        "input": "ForebrainGlut/ForebrainGlut.loom",
        "label_key": "Clusters",
        "source_embedding": None,
        "embed_key": "X_pca",
        "velocity_basis": "pca",
        "preprocessed": False,
    },
    "Gastrulation": {
        "input": "Gastrulation/Gastrulation.h5ad",
        "label_key": "celltype",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": False,
    },
    "Pancreas": {
        "input": "Pancreas/Pancreas.h5ad",
        "label_key": "clusters",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": False,
    },
    "dentategyrus": {
        "input": "dentategyrus/dentategyrus.h5ad",
        "label_key": "clusters",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": False,
    },
    "hematopoiesis_dynamo": {
        "input": "hematopoiesis_dynamo/hematopoiesis_dynamo.h5ad",
        "label_key": "cell_type",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": False,
    },
    "mousebrain": {
        "input": "mousebrain/mousebrain.h5ad",
        "label_key": "celltype",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": True,
        "expected_shape": (3365, 936),
    },
    "neuron": {
        "input": "neuron/neuron.h5ad",
        "label_key": "time_str",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": True,
        "expected_shape": (3066, 2000),
    },
    "organoid": {
        "input": "organoid/organoid_pre.h5ad",
        "label_key": "cell_type",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": False,
    },
    "retina": {
        "input": "retina/retina_full.h5ad",
        "label_key": "Annotation",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": False,
    },
    "rpe1": {
        "input": "rpe1/rpe1.h5ad",
        "label_key": "cell_cycle_phase",
        "source_embedding": "X_RFP_GFP",
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": False,
    },
    "rpe1_gen": {
        "input": "rpe1_gen/rpe1_gen.h5ad",
        "label_key": "phase",
        "source_embedding": None,
        "embed_key": "X_umap",
        "velocity_basis": "umap",
        "preprocessed": True,
        "expected_shape": (3354, 2000),
    },
}

DATASETS.update({key: {**value, "source_embedding": None, "embed_key": "X_umap",
                        "velocity_basis": "umap"}
                 for key, value in simulation_inputs().items()})

MIN_SHARED_COUNTS = 20
N_TOP_GENES = 2000
CELL_DANCER_N_JOBS = 12
GRAPH_N_JOBS = 8
PERMUTATION_RATIO = 0.125


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / dataset
    output_path = OUTPUT_ROOT / f"celldancer_{dataset}.h5ad"
    if output_path.exists():
        raise FileExistsError(f"output exists: {output_path}; choose a new output root")
    dataframe_path = output_dir / "cell_type_u_s.csv"
    fit_dir = output_dir / "fit"

    print(f"\n=== cellDancer: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    adata = sc.read(input_path, cache=False)
    adata.var_names_make_unique()

    missing_layers = {"spliced", "unspliced"} - set(adata.layers)
    if missing_layers:
        raise KeyError(f"{dataset}: missing count layers {sorted(missing_layers)}")
    if config["label_key"] not in adata.obs:
        raise KeyError(f"{dataset}: missing label key {config['label_key']!r}")
    adata.obs["clusters"] = adata.obs[config["label_key"]].copy()

    source_embedding = config["source_embedding"]
    if source_embedding is not None:
        if source_embedding not in adata.obsm:
            raise KeyError(f"{dataset}: missing embedding key {source_embedding!r}")
        adata.obsm["X_umap"] = adata.obsm[source_embedding].copy()

    if config["preprocessed"]:
        if adata.shape != config["expected_shape"]:
            raise ValueError(
                f"{dataset}: expected fixed benchmark shape {config['expected_shape']}, "
                f"found {adata.shape}"
            )
        missing_moments = {"Ms", "Mu"} - set(adata.layers)
        if missing_moments:
            raise KeyError(f"{dataset}: missing fixed moments {sorted(missing_moments)}")
        print("using the fixed scVelo-preprocessed benchmark object", flush=True)
    else:
        scv.pp.filter_and_normalize(
            adata,
            min_shared_counts=MIN_SHARED_COUNTS,
            n_top_genes=N_TOP_GENES,
        )
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    embed_key = config["embed_key"]
    if embed_key not in adata.obsm:
        raise KeyError(f"{dataset}: missing cellDancer embedding {embed_key!r}")
    if adata.obsm[embed_key].shape[1] < 2:
        raise ValueError(f"{dataset}: {embed_key!r} must have at least two dimensions")

    # Official cellDancer input boundary: Ms/Mu, a label column, and the first
    # two dimensions of the declared embedding are written as a dataframe.
    output_dir.mkdir(parents=True, exist_ok=True)
    fit_dir.mkdir(parents=True, exist_ok=True)
    cell_dancer_df = cd.adata_to_df_with_embed(
        adata,
        us_para=["Mu", "Ms"],
        cell_type_para="clusters",
        embed_para=embed_key,
        save_path=str(dataframe_path),
    )
    loss_df, cell_dancer_df = cd.velocity(
        cell_dancer_df,
        permutation_ratio=PERMUTATION_RATIO,
        n_jobs=CELL_DANCER_N_JOBS,
        save_path=str(fit_dir),
    )
    if cell_dancer_df is None or len(cell_dancer_df) == 0:
        raise RuntimeError(
            f"{dataset}: cellDancer rejected every gene; inspect the Ms/Mu "
            "input and its per-gene quality diagnostics before continuing"
        )
    print(f"cellDancer fit rows: {len(loss_df)}", flush=True)

    # ``to_dynamo`` materializes cellDancer's fitted alpha/beta/gamma and
    # unspliced velocity. Copy them back by names, leaving any rejected genes
    # as NaN while retaining the benchmark's complete gene panel.
    adata_from_dancer = cd.utilities.to_dynamo(cell_dancer_df)
    cell_positions = adata_from_dancer.obs_names.get_indexer(adata.obs_names)
    if (cell_positions < 0).any():
        raise KeyError(f"{dataset}: cellDancer result is missing benchmark cells")
    common_genes = adata.var_names[adata.var_names.isin(adata_from_dancer.var_names)]
    if common_genes.empty:
        raise ValueError(f"{dataset}: cellDancer did not return any benchmark genes")
    input_positions = adata.var_names.get_indexer(common_genes)
    dancer_positions = adata_from_dancer.var_names.get_indexer(common_genes)

    for layer in ("celldancer_alpha", "celldancer_beta", "celldancer_gamma", "velocity_u"):
        adata.layers[layer] = np.full(adata.shape, np.nan)
    for target_layer, dancer_layer in (
        ("celldancer_alpha", "alpha"),
        ("celldancer_beta", "beta"),
        ("celldancer_gamma", "gamma"),
        ("velocity_u", "velocity_U"),
    ):
        values = adata_from_dancer.layers[dancer_layer]
        if sparse.issparse(values):
            values = values.toarray()
        adata.layers[target_layer][:, input_positions] = np.asarray(values)[
            cell_positions[:, None], dancer_positions
        ]

    adata = cd.export_velocity_to_dynamo(cell_dancer_df, adata)
    velocity_s = adata.layers.pop("velocity_S")
    adata.layers["velocity"] = velocity_s.toarray() if sparse.issparse(velocity_s) else np.asarray(velocity_s)

    velocity_basis = config["velocity_basis"]
    scv.tl.velocity_graph(adata, vkey="velocity", n_jobs=GRAPH_N_JOBS)
    scv.tl.velocity_embedding(adata, basis=velocity_basis, vkey="velocity")
    scv.settings.figdir = str(output_dir)
    scv.pl.velocity_embedding_stream(
        adata,
        basis=velocity_basis,
        vkey="velocity",
        color="clusters",
        show=False,
        save=f"_{dataset}_stream.png",
    )
    plt.close("all")

    required_layers = {
        "spliced",
        "unspliced",
        "velocity",
        "velocity_u",
        "celldancer_alpha",
        "celldancer_beta",
        "celldancer_gamma",
    }
    velocity_embedding_key = f"velocity_{velocity_basis}"
    missing_layers = required_layers - set(adata.layers)
    missing_obsm = {embed_key, velocity_embedding_key} - set(adata.obsm)
    if missing_layers or missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete cellDancer result; "
            f"missing layers={sorted(missing_layers)}, obsm={sorted(missing_obsm)}"
        )

    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    global DATA_ROOT, OUTPUT_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--dataset", dest="datasets", action="append", choices=sorted(DATASETS))
    args = parser.parse_args()
    DATA_ROOT = args.data_root.expanduser().resolve()
    OUTPUT_ROOT = args.output_root.expanduser().resolve()
    datasets = args.datasets or DATASETS_TO_RUN
    if scv.__version__ != "0.2.4":
        raise RuntimeError(f"scVelo 0.2.4 is required, found {scv.__version__}")
    unknown_datasets = set(datasets) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets: {sorted(unknown_datasets)}")

    for dataset in datasets:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
