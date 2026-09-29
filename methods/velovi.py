#!/usr/bin/env python3
"""Reproduce C12 (veloVI) with the upstream tutorial workflow.

Select the historical Cell or Gene workflow explicitly. Cell uses the package
training defaults and per-gene time scaling; Gene uses the separate historical
training settings and exports unscaled spliced and unspliced velocities.
Shared preprocessing and the absence of an extra global seed are retained.
"""

import argparse
from inspect import signature
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
import torch
import velovi
from velovi import VELOVI, preprocess_data


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"

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

# Gene's historical training settings. Cell retains the bare vae.train() call
# (velovi 0.3.1 defaults: max_epochs=500, batch_size=256).
GENE_MAX_EPOCHS = 1000
GENE_GPU_BATCH_SIZE = 64
GENE_CPU_BATCH_SIZE = 32
N_POSTERIOR_SAMPLES = 25
N_JOBS = 8


def run_dataset(dataset, evaluation_role, data_root, output_root):
    if evaluation_role not in {"cell", "gene"}:
        raise ValueError(f"unknown veloVI evaluation role: {evaluation_role}")
    config = DATASETS[dataset]
    input_path = data_root / config["input"]
    output_dir = output_root / dataset
    output_path = output_dir / f"veloVI_{dataset}_all.h5ad"

    print(f"\n=== veloVI {evaluation_role}: {dataset} ===", flush=True)
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
    if evaluation_role == "cell":
        train_kwargs = {}
    else:
        train_kwargs = {
            "max_epochs": GENE_MAX_EPOCHS,
            "batch_size": GENE_GPU_BATCH_SIZE if torch.cuda.is_available() else GENE_CPU_BATCH_SIZE,
        }
    vae.train(**train_kwargs)

    if evaluation_role == "cell":
        # Preserve the historical posterior-draw order and tutorial scaling.
        latent_time = np.asarray(vae.get_latent_time(n_samples=N_POSTERIOR_SAMPLES))
        if latent_time.shape != adata.shape:
            raise ValueError(f"{dataset}: latent time must have shape {adata.shape}")
        max_time = latent_time.max(axis=0)
        if not np.all(np.isfinite(max_time)) or np.any(max_time <= 0):
            raise ValueError(f"{dataset}: per-gene maximum latent time must be finite and positive")
        scaling = 20 / max_time
        adata.layers["velocity"] = np.asarray(
            vae.get_velocity(n_samples=N_POSTERIOR_SAMPLES, velo_statistic="mean")
        ) / scaling
        adata.layers["latent_time_velovi"] = latent_time
        rates = vae.get_rates()
        for rate in ("alpha", "beta", "gamma"):
            adata.var[f"fit_{rate}"] = np.asarray(rates[rate]) / scaling
        adata.var["fit_t_"] = (
            torch.nn.functional.softplus(vae.module.switch_time_unconstr)
            .detach().cpu().numpy()
        ) * scaling
        adata.layers["fit_t"] = latent_time * scaling
        adata.var["fit_scaling"] = 1.0
    else:
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

    defaults = signature(VELOVI.train).parameters
    adata.uns["benchmark_run"] = {
        "config_id": "C12",
        "evaluation_role": evaluation_role,
        "dataset": dataset,
        "velovi_version": velovi.__version__,
        "max_epochs": train_kwargs.get("max_epochs", defaults["max_epochs"].default),
        "batch_size": train_kwargs.get("batch_size", defaults["batch_size"].default),
        "n_posterior_samples": N_POSTERIOR_SAMPLES,
        "velocity_scaling": "divide_by_20_over_gene_max_latent_time" if evaluation_role == "cell" else "none",
    }

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

    required_layers = {"spliced", "unspliced", "Ms", "Mu", "velocity"}
    if evaluation_role == "cell":
        required_layers.update({"latent_time_velovi", "fit_t"})
        if "velocity_u" in adata.layers:
            raise RuntimeError(f"{dataset}: Cell-role output unexpectedly contains velocity_u")
    else:
        required_layers.add("velocity_u")
    missing_layers = required_layers - set(adata.layers)
    missing_obsm = {"X_umap", "velocity_umap"} - set(adata.obsm)
    if missing_layers or missing_obsm or "velocity_graph" not in adata.uns:
        raise RuntimeError(
            f"{dataset}: incomplete veloVI result; missing layers={sorted(missing_layers)}, "
            f"obsm={sorted(missing_obsm)}, velocity_graph={'velocity_graph' in adata.uns}"
        )

    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-role", choices=("cell", "gene"), required=True)
    parser.add_argument("--output-root", type=Path, required=True, help="empty directory for this role's run")
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--dataset", dest="datasets", action="append", choices=sorted(DATASETS))
    args = parser.parse_args()
    datasets = args.datasets or DATASETS_TO_RUN
    output_root = args.output_root.expanduser().resolve()
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(f"output root is not empty: {output_root}; choose a new run directory")

    if velovi.__version__ != "0.3.1":
        raise RuntimeError(f"velovi 0.3.1 is required, found {velovi.__version__}")
    unknown_datasets = set(datasets) - set(DATASETS)
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

    for dataset in datasets:
        run_dataset(dataset, args.evaluation_role, args.data_root.expanduser().resolve(), output_root)


if __name__ == "__main__":
    main()
