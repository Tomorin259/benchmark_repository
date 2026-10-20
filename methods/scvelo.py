#!/usr/bin/env python3
"""Reproduce C01--C03 with the scVelo tutorial workflow.

Set DATA_ROOT to the archived input artifacts before running. The four
already-preprocessed benchmark objects are intentionally not preprocessed a
second time; every other dataset follows scVelo's
filter_and_normalize -> moments sequence.
"""

import argparse
from pathlib import Path
import sys

# Direct execution must import installed packages, not same-named launchers.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)

import scanpy as sc
import scvelo as scv


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))
from methods.simulation_inputs import simulation_inputs
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs"

# Edit these two lists for a focused rerun.
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
MODES_TO_RUN = ["deterministic", "stochastic", "dynamical"]

# Fixed dataset inputs and annotation mappings. All 14 inputs
# already provide spliced and unspliced. label_key only supplies the standard
# clusters field used by the benchmark's method-native result files.
# embedding_key=None means that moments will create X_pca and that X_pca is
# then retained as the benchmark embedding, as in the legacy Forebrain run.
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

OUTPUT_PREFIX = {
    "deterministic": "scv_deterministic",
    "stochastic": "scv_stochastic",
    "dynamical": "scv_dynamic",
}
N_JOBS = 8


def run_dataset(dataset, mode, data_root, output_root):
    config = DATASETS[dataset]
    input_path = data_root / config["input"]
    operational_id = "scvelo_dynamic" if mode == "dynamical" else f"scvelo_{mode}"
    output_dir = output_root / operational_id
    output_path = output_dir / f"{OUTPUT_PREFIX[mode]}_{dataset}.h5ad"
    if output_path.exists():
        raise FileExistsError(f"output exists: {output_path}; choose a new output root")

    print(f"\n=== scVelo {mode}: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    # scVelo 0.3.3 delegates file loading to Scanpy.
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
        print("using the fixed scVelo-preprocessed benchmark object", flush=True)
    else:
        # scVelo tutorial preprocessing for raw benchmark inputs.
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            raise KeyError(f"{dataset}: preprocessing did not provide X_pca for the benchmark embedding")
        adata.obsm["X_umap"] = adata.obsm["X_pca"].copy()

    if mode == "dynamical":
        if config.get("simulation"):
            adata.var["velocity_genes"] = True
        scv.tl.recover_dynamics(adata, n_jobs=N_JOBS)
    scv.tl.velocity(adata, mode=mode)
    scv.tl.velocity_graph(adata, n_jobs=N_JOBS)
    # Compute the saved velocity embedding explicitly rather than through a
    # plotting side effect.
    scv.tl.velocity_embedding(adata, basis="umap")

    missing_layers = {"spliced", "unspliced", "velocity"} - set(adata.layers)
    missing_obsm = {"X_umap", "velocity_umap"} - set(adata.obsm)
    if missing_layers or missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete scVelo result; missing layers={sorted(missing_layers)}, "
            f"obsm={sorted(missing_obsm)}, velocity_graph={'velocity_graph' in adata.uns}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    adata.uns["benchmark_run"] = {
        "operational_id": operational_id, "mode": mode, "dataset": dataset,
        "scvelo_version": scv.__version__,
    }
    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", dest="datasets", action="append", choices=sorted(DATASETS))
    parser.add_argument("--mode", dest="modes", action="append", choices=sorted(OUTPUT_PREFIX))
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    args = parser.parse_args()
    datasets = args.datasets or DATASETS_TO_RUN
    modes = args.modes or MODES_TO_RUN
    if any(DATASETS[name].get("simulation") for name in datasets) and modes != ["dynamical"]:
        raise ValueError("Simulation benchmark uses --mode dynamical")
    if scv.__version__ != "0.3.3":
        raise RuntimeError(f"scVelo 0.3.3 is required, found {scv.__version__}")
    unknown_datasets = set(datasets) - set(DATASETS)
    unknown_modes = set(modes) - set(OUTPUT_PREFIX)
    if unknown_datasets or unknown_modes:
        raise ValueError(
            f"unknown datasets={sorted(unknown_datasets)}; unknown modes={sorted(unknown_modes)}"
        )

    # These tutorial settings also affect the automatic scale stored by
    # velocity_embedding, so they are part of the numerical reproduction.
    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.set_figure_params("scvelo")
    print(f"scVelo version: {scv.__version__}", flush=True)
    for mode in modes:
        for dataset in datasets:
            run_dataset(dataset, mode, args.data_root.expanduser().resolve(), args.output_root.expanduser().resolve())


if __name__ == "__main__":
    main()
