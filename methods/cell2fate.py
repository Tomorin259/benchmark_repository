#!/usr/bin/env python3
"""Reproduce C09: Cell2fate on the 14 benchmark datasets.

This follows the Cell2fate tutorial order: select training genes from raw
spliced/unspliced counts, set up the model, train, export posterior summaries,
and compute the total RNA velocity with Cell2fate's scVelo routine.
"""

from importlib.metadata import version
from pathlib import Path

import cell2fate as c2f
import numpy as np
import scanpy as sc
import scvelo as scv
from scipy import sparse


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
RESULTS_ROOT = REPOSITORY_ROOT.parent / "results"
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

# Inputs ending in _counts contain raw UMI counts. For the three full-gene
# inputs, template preserves the benchmark cells, labels and visualization
# embedding before Cell2fate selects its own 2,000 training genes.
DATASETS = {
    "bonemarrow": {
        "input": "bonemarrow/bonemarrow.h5ad",
        "label_key": "clusters",
        "embedding_key": "X_tsne",
    },
    "capo_pre": {
        "input": "capo_pre/capo_pre_counts.h5ad",
        "template": "capo_pre/capo_pre.h5ad",
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
        "embedding_source": "scvelo_dynamic/scv_dynamic_ForebrainGlut.h5ad",
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
        "input": "mousebrain/mousebrain_counts.h5ad",
        "label_key": "celltype",
        "embedding_key": "X_umap",
    },
    "neuron": {
        "input": "neuron/neuron_counts.h5ad",
        "template": "neuron/neuron.h5ad",
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
    },
    "rpe1_gen": {
        "input": "rpe1_gen/rpe1_gen_counts.h5ad",
        "template": "rpe1_gen/rpe1_gen.h5ad",
        "label_key": "phase",
        "embedding_key": "X_umap",
    },
}


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / "cell2fate" / dataset
    output_path = output_dir / f"cell2fate_{dataset}.h5ad"

    print(f"\n=== Cell2fate: {dataset} ===", flush=True)
    print(f"raw input: {input_path}", flush=True)
    if input_path.suffix == ".loom":
        adata = sc.read_loom(input_path, sparse=True)
    else:
        adata = sc.read_h5ad(input_path)
    adata.var_names_make_unique()

    if adata.obs_names.has_duplicates or adata.var_names.has_duplicates:
        raise ValueError(f"{dataset}: cell or gene names are not unique")

    if "template" in config:
        template = sc.read_h5ad(DATA_ROOT / config["template"])
        positions = adata.obs_names.get_indexer(template.obs_names)
        if (positions < 0).any():
            missing = template.obs_names[positions < 0][:10].tolist()
            raise KeyError(f"{dataset}: raw counts are missing template cells {missing}")
        adata = adata[positions].copy()
        adata.obs = template.obs.copy()
        if "X_umap" in template.obsm:
            adata.obsm["X_umap"] = template.obsm["X_umap"].copy()

    for layer in ("spliced", "unspliced"):
        if layer not in adata.layers:
            raise KeyError(f"{dataset}: missing layers[{layer!r}]")
        matrix = adata.layers[layer]
        values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
        finite = bool(np.isfinite(values).all())
        nonnegative = finite and bool((values >= 0).all())
        integer_valued = finite and bool(
            np.allclose(values, np.rint(values), atol=1e-8, rtol=0)
        )
        deviation = (
            float(np.max(np.abs(values.astype(np.float64) - np.rint(values))))
            if values.size
            else 0.0
        )
        print(
            f"[counts] {layer}: dtype={matrix.dtype} shape={matrix.shape} "
            f"finite={finite} nonnegative={nonnegative} "
            f"integer_valued={integer_valued} "
            f"max_abs_fractional_deviation={deviation:.3g}",
            flush=True,
        )
        if not (finite and nonnegative and integer_valued):
            raise ValueError(f"{dataset}.{layer} is not raw integer-valued counts")
        if not np.issubdtype(matrix.dtype, np.integer):
            maximum = float(values.max()) if values.size else 0.0
            dtype = np.int32 if maximum <= np.iinfo(np.int32).max else np.int64
            converted = matrix.astype(dtype)
            if sparse.issparse(matrix):
                difference = matrix.tocsr().astype(np.float64) - converted.tocsr().astype(
                    np.float64
                )
                difference.eliminate_zeros()
                preserved = difference.nnz == 0
            else:
                preserved = np.array_equal(np.asarray(matrix), np.asarray(converted))
            if not preserved:
                raise AssertionError(f"{dataset}.{layer} changed during dtype conversion")
            adata.layers[layer] = converted
            print(f"[counts] {layer}: lossless {np.dtype(dtype)} conversion", flush=True)

    label_key = config["label_key"]
    if label_key not in adata.obs:
        raise KeyError(f"{dataset}: missing label key {label_key!r}")
    adata.obs["clusters"] = adata.obs[label_key].astype("category")

    if "embedding_key" in config:
        embedding_key = config["embedding_key"]
        if embedding_key not in adata.obsm:
            raise KeyError(f"{dataset}: missing embedding key {embedding_key!r}")
        adata.obsm["X_umap"] = adata.obsm[embedding_key].copy()
    else:
        embedding_reference = sc.read_h5ad(RESULTS_ROOT / config["embedding_source"])
        positions = embedding_reference.obs_names.get_indexer(adata.obs_names)
        if (positions < 0).any():
            raise KeyError(f"{dataset}: external embedding is missing cells")
        adata.obsm["X_umap"] = embedding_reference.obsm["X_pca"][positions].copy()

    print(f"training input before gene selection: {adata.shape}", flush=True)
    adata = c2f.utils.get_training_data(
        adata,
        cells_per_cluster=10**5,
        cluster_column="clusters",
        remove_clusters=[],
        min_shared_counts=20,
        n_var_genes=2000,
    )
    print(f"training input after gene selection: {adata.shape}", flush=True)

    c2f.Cell2fate_DynamicalModel.setup_anndata(
        adata,
        spliced_label="spliced",
        unspliced_label="unspliced",
    )
    n_modules = c2f.utils.get_max_modules(adata)
    model = c2f.Cell2fate_DynamicalModel(adata, n_modules=n_modules)
    model.train()
    adata = model.export_posterior(adata)

    output_dir.mkdir(parents=True, exist_ok=True)
    model.compute_and_plot_total_velocity_scvelo(
        adata,
        save=str(output_dir / f"{dataset}_velocity.png"),
        delete=False,
    )

    required_layers = {"velocity"}
    required_obsm = {"X_umap", "velocity_umap"}
    required_uns = {"velocity_graph", "velocity_graph_neg"}
    missing_layers = required_layers - set(adata.layers)
    missing_obsm = required_obsm - set(adata.obsm)
    missing_uns = required_uns - set(adata.uns)
    if missing_layers or missing_obsm or missing_uns:
        raise RuntimeError(
            f"{dataset}: incomplete Cell2fate result; "
            f"layers={sorted(missing_layers)}, obsm={sorted(missing_obsm)}, "
            f"uns={sorted(missing_uns)}"
        )

    for key in list(adata.layers):
        value = adata.layers[key]
        if hasattr(value, "detach"):
            adata.layers[key] = value.detach().cpu().numpy()
        elif hasattr(value, "numpy"):
            adata.layers[key] = value.numpy()

    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    expected_versions = {
        "cell2fate": "0.1a0",
        "scanpy": "1.9.1",
        "scvelo": "0.2.4",
    }
    actual_versions = {name: version(name) for name in expected_versions}
    if actual_versions != expected_versions:
        raise RuntimeError(
            f"Cell2fate environment mismatch: expected={expected_versions}, "
            f"found={actual_versions}"
        )

    scv.settings.verbosity = 3
    scv.set_figure_params("scvelo")

    unknown_datasets = set(DATASETS_TO_RUN) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets={sorted(unknown_datasets)}")
    for dataset in DATASETS_TO_RUN:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
