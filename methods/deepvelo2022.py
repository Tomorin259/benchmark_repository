#!/usr/bin/env python3
"""Reproduce C18, DeepVelo's scVelo-dynamics reconstruction configuration.

This is intentionally separate from C16's graph-convolutional DeepVelo.
It follows the archived DeepVelo notebook-style workflow: fit scVelo
dynamical velocities, train the published TensorFlow reconstruction network,
then evaluate its reconstructed velocity field on scVelo's selected velocity
genes.
"""

import argparse
from pathlib import Path
import sys
import time

# Direct execution must import installed packages, not same-named launchers.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)

import numpy as np
import scanpy as sc
import scvelo as scv
from scipy import sparse
from sklearn.model_selection import train_test_split


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))
from methods.simulation_inputs import simulation_inputs
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs" / "deepvelo2022"

# Edit this list for a focused rerun. The train/test split uses random_state=42.
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

# Dataset and name mappings are kept visible because this C18 result lineage
# historically used the source-file stems ``organoid_pre`` and ``retina_full``.
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
        "expected_shape": (3354, 2000),
    },
}

DATASETS.update(simulation_inputs())

N_JOBS = 8


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    output_dir = OUTPUT_ROOT / dataset
    output_stem = config.get("output_stem", dataset)
    output_path = output_dir / f"deepvelo_{output_stem}_recon.h5ad"
    if output_path.exists():
        raise FileExistsError(f"output exists: {output_path}; choose a new output root")
    model_path = output_dir / f"deepvelo_{output_stem}.h5"

    print(f"\n=== DeepVelo reconstruction: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    if input_path.suffix.lower() == ".loom":
        adata = sc.read_loom(input_path.as_posix(), sparse=True)
    else:
        adata = sc.read_h5ad(input_path.as_posix())
    adata_raw = adata.copy()

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
        # The original DeepVelo notebook carries a non-logged companion object
        # through gene selection to construct the reconstruction-network input.
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.filter_and_normalize(
            adata_raw,
            log=False,
            min_shared_counts=20,
            n_top_genes=2000,
        )
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            raise KeyError(f"{dataset}: preprocessing did not provide X_pca for the benchmark embedding")
        adata.obsm["X_umap"] = adata.obsm["X_pca"].copy()

    # The upstream reconstruction target is a scVelo dynamical velocity.
    if config.get("simulation"):
        adata.var["velocity_genes"] = True
    scv.tl.recover_dynamics(adata, n_jobs=N_JOBS)
    if config.get("simulation"):
        adata.var["velocity_genes"] = True
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata, n_jobs=N_JOBS)

    output_dir.mkdir(parents=True, exist_ok=True)
    scv.settings.figdir = str(output_dir)
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        save=f"_stream_{output_stem}_scvelo.png",
    )

    if "velocity_genes" not in adata.var:
        raise KeyError(f"{dataset}: scVelo did not produce var['velocity_genes']")
    if config.get("simulation"):
        adata.var["velocity_genes"] = True
    velocity_genes = adata.var["velocity_genes"].to_numpy(dtype=bool)
    if not velocity_genes.any():
        raise ValueError(f"{dataset}: scVelo selected no velocity genes")

    X_raw = adata_raw.X.toarray() if sparse.issparse(adata_raw.X) else np.asarray(adata_raw.X)
    X = np.tile(X_raw[:, velocity_genes], (5, 1))
    Y = np.tile(adata.layers["velocity"][:, velocity_genes], (5, 1))
    noise_sigma = (X_raw.std() / 70) ** 2
    X[adata_raw.n_obs :, :] += np.random.normal(0, noise_sigma, X[adata_raw.n_obs :, :].shape)
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=42)

    # TensorFlow sees the card selected by CUDA_VISIBLE_DEVICES.  Restricting
    # it to its first visible GPU and enabling memory growth is the upstream
    # notebook's device setup, expressed without a hard-coded physical index.
    import tensorflow as tf
    from tensorflow import keras

    gpus = tf.config.list_physical_devices("GPU")
    if gpus:
        tf.config.experimental.set_visible_devices(gpus[0], "GPU")
        tf.config.experimental.set_memory_growth(gpus[0], True)

    input_dim = keras.Input(shape=(X.shape[1],))
    encoded = keras.layers.Dense(64, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(input_dim)
    encoded = keras.layers.Dense(64, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(encoded)
    encoded = keras.layers.Dense(64, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(encoded)
    encoded = keras.layers.Dense(16, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(encoded)
    decoded = keras.layers.Dense(16, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(encoded)
    decoded = keras.layers.Dense(64, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(decoded)
    decoded = keras.layers.Dense(64, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(decoded)
    decoded = keras.layers.Dense(64, activity_regularizer=keras.regularizers.l1(1e-6))(decoded)
    decoded = keras.layers.Dense(Y.shape[1])(decoded)

    autoencoder = keras.Model(input_dim, decoded)
    autoencoder.compile(optimizer=keras.optimizers.Adam(learning_rate=0.00005), loss="mse")
    start = time.time()
    autoencoder.fit(
        X_train,
        y_train,
        epochs=50,
        batch_size=4,
        shuffle=True,
        validation_data=(X_test, y_test),
        callbacks=[keras.callbacks.EarlyStopping(monitor="val_loss", patience=3)],
    )
    autoencoder.save(model_path.as_posix())

    adata_recon = adata[:, velocity_genes].copy()
    adata_recon.layers["velocity_reconstructed"] = autoencoder.predict(X_raw[:, velocity_genes])
    scv.tl.velocity_graph(adata_recon, vkey="velocity_reconstructed", n_jobs=N_JOBS)
    scv.pl.velocity_embedding_stream(
        adata_recon,
        basis="umap",
        color="clusters",
        vkey="velocity_reconstructed",
        save=f"figure1_{output_stem}.png",
    )

    required_layers = {"spliced", "unspliced", "Ms", "Mu", "velocity", "velocity_reconstructed"}
    required_obsm = {"X_umap", "velocity_umap", "velocity_reconstructed_umap"}
    missing_layers = required_layers - set(adata_recon.layers)
    missing_obsm = required_obsm - set(adata_recon.obsm)
    if missing_layers or missing_obsm or "velocity_reconstructed_graph" not in adata_recon.uns:
        raise RuntimeError(
            f"{dataset}: incomplete DeepVelo reconstruction result; "
            f"missing layers={sorted(missing_layers)}, obsm={sorted(missing_obsm)}, "
            f"velocity_reconstructed_graph={'velocity_reconstructed_graph' in adata_recon.uns}"
        )

    adata_recon.write_h5ad(output_path.as_posix())
    print(f"wrote: {output_path}", flush=True)
    print(f"completed {dataset} in {time.time() - start:.1f}s", flush=True)


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

    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.set_figure_params("scvelo")
    print(f"scVelo version: {scv.__version__}", flush=True)
    print("GPU selection: inherited from CUDA_VISIBLE_DEVICES", flush=True)
    for dataset in datasets:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
