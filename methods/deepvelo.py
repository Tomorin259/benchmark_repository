#!/usr/bin/env python3
"""Reproduce C16 (DeepVelo GCN) with the benchmark's direct workflow.

Inputs and preprocessing follow the dataset configuration below.
Select the device with ``CUDA_VISIBLE_DEVICES``.
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# When this file is executed directly, Python puts ``methods/`` first on
# sys.path and shadows the installed DeepVelo package with this script.
# Remove only that script directory; the package must come from the audited
# environment.
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) in sys.path:
    sys.path.remove(str(SCRIPT_DIR))

import numpy as np
import scanpy as sc
import scvelo as scv
import torch

from deepvelo import Constants, train
from deepvelo.utils import update_dict, velocity
from deepvelo.utils.preprocess import autoset_coeff_s


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))
from methods.simulation_inputs import simulation_inputs
DATA_ROOT = REPOSITORY_ROOT.parent / "data"

# Edit this list for a focused rerun; the seed is configured below.
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

# Dataset inputs, preprocessing settings and annotation mappings.
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

SEED = 123
N_JOBS = 8


def seed_everything():
    """Reset all RNGs used by the historical per-dataset launcher."""
    torch.manual_seed(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    np.random.seed(SEED)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument(
        "--evaluation-role",
        choices=("cell", "gene", "kinetics"),
        required=True,
        help="paper analysis role; this controls the DeepVelo output head",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="empty directory receiving H5AD files, plots, and training logs",
    )
    parser.add_argument(
        "--dataset",
        dest="datasets",
        action="append",
        choices=sorted(DATASETS),
        help="run one dataset; repeat for a focused run (default: all 14)",
    )
    parser.add_argument(
        "--reset-seed-per-dataset",
        action="store_true",
        help="reset seed=123 at the start of every dataset, matching the historical launcher",
    )
    return parser.parse_args()


def run_dataset(dataset, output_root, pred_unspliced, evaluation_role, reset_seed_per_dataset):
    if reset_seed_per_dataset or DATASETS[dataset].get("simulation"):
        seed_everything()
    config = DATASETS[dataset]
    if config.get("simulation") and evaluation_role != "gene":
        raise ValueError("Simulation benchmark uses --evaluation-role gene")
    input_path = DATA_ROOT / config["input"]
    output_dir = output_root / dataset
    output_path = output_dir / f"deepvelo_{dataset}.h5ad"

    print(f"\n=== DeepVelo: {dataset} ===", flush=True)
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
        print("input moments: Ms, Mu", flush=True)
    else:
        # DeepVelo's official tutorial starts from this scVelo sequence.
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            raise KeyError(f"{dataset}: preprocessing did not provide X_pca for the benchmark embedding")
        adata.obsm["X_umap"] = adata.obsm["X_pca"].copy()

    # C16 selects the output head by evaluation role and stores logs and
    # checkpoints in the dataset-specific output directory.
    configs = {
        "name": "DeepVelo",
        "loss": {"args": {"coeff_s": autoset_coeff_s(adata)}},
        "arch": {
            "type": "VeloGCN",
            "args": {
                "layers": [64, 64],
                "dropout": 0.2,
                "fc_layer": False,
                "pred_unspliced": pred_unspliced,
            },
        },
    }
    configs = update_dict(Constants.default_configs, configs)
    configs["trainer"]["save_dir"] = str(output_dir / "training")

    # Official DeepVelo training starts from its initial velocity estimate.
    velocity(adata, mask_zero=False)
    train(adata, configs)
    if config.get("simulation"):
        adata.layers["velocity_u"] = adata.layers["velocity_unspliced"].copy()

    scv.tl.velocity_graph(adata, n_jobs=N_JOBS)
    output_dir.mkdir(parents=True, exist_ok=True)
    scv.settings.figdir = str(output_dir)
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="clusters",
        legend_fontsize=9,
        save=f"_stream_{dataset}.png",
        dpi=150,
    )

    required_layers = {
        "spliced",
        "unspliced",
        "Ms",
        "Mu",
        "velocity",
        "cell_specific_beta",
        "cell_specific_gamma",
        "variance_velocity",
    }
    if pred_unspliced:
        required_layers.update({"velocity_unspliced", "cell_specific_alpha"})
    required_obsm = {"X_umap", "velocity_umap"}
    missing_layers = required_layers - set(adata.layers)
    missing_obsm = required_obsm - set(adata.obsm)
    if missing_layers or missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete DeepVelo result; missing layers={sorted(missing_layers)}, "
            f"obsm={sorted(missing_obsm)}, velocity_graph={'velocity_graph' in adata.uns}"
        )
    if evaluation_role == "cell" and "velocity_unspliced" in adata.layers:
        raise RuntimeError(
            f"{dataset}: cell-level run unexpectedly produced velocity_unspliced"
        )

    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    global DATA_ROOT
    args = parse_args()
    DATA_ROOT = args.data_root.expanduser().resolve()
    if scv.__version__ != "0.2.5":
        raise RuntimeError(f"scVelo 0.2.5 is required, found {scv.__version__}")
    datasets = args.datasets or DATASETS_TO_RUN
    unknown_datasets = set(datasets) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets: {sorted(unknown_datasets)}")

    output_root = args.output_root.resolve()
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(
            f"output root is not empty: {output_root}; choose a new run directory"
        )
    output_root.mkdir(parents=True, exist_ok=True)
    pred_unspliced = args.evaluation_role in {"gene", "kinetics"}

    if not args.reset_seed_per_dataset:
        seed_everything()
    scv.settings.verbosity = 3
    scv.set_figure_params("scvelo", transparent=False)
    print(f"DeepVelo source: {train.__module__}", flush=True)
    print(f"global RNG seed: {SEED} (C16 historical configuration)", flush=True)
    print(
        f"evaluation role: {args.evaluation_role}; "
        f"pred_unspliced: {pred_unspliced}",
        flush=True,
    )
    print("GPU selection: inherited from CUDA_VISIBLE_DEVICES", flush=True)

    run_metadata = {
        "config_id": "C16",
        "operational_id": "deepvelo",
        "evaluation_role": args.evaluation_role,
        "pred_unspliced": pred_unspliced,
        "reset_seed_per_dataset": args.reset_seed_per_dataset,
        "simulation_seed_scope": "per_dataset",
        "seed": SEED,
        "n_jobs": N_JOBS,
        "datasets": datasets,
        "output_root": str(output_root),
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    (output_root / "run_config.json").write_text(
        json.dumps(run_metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    for dataset in datasets:
        run_dataset(
            dataset,
            output_root,
            pred_unspliced,
            args.evaluation_role,
            args.reset_seed_per_dataset,
        )

    run_metadata["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    (output_root / "run_config.json").write_text(
        json.dumps(run_metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
