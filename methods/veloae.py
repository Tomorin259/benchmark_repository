#!/usr/bin/env python3
"""Reproduce C11: VeloAE latent-space velocity on the 14 datasets."""

from importlib.metadata import version
from pathlib import Path

import numpy as np
import scanpy as sc
import scvelo as scv
import torch
from scipy import sparse
from veloproj import estimate_ld_velocity, fit_model, get_parser, init_model, new_adata


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs"
GPU_INDEX = 3

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

# The four listed entries already contain the scVelo moments required by
# VeloAE. The remaining ten inputs are preprocessed with the parameters used
# by the original benchmark script.
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
    output_dir = OUTPUT_ROOT / "veloae"
    checkpoint_dir = output_dir / "checkpoints"
    output_path = output_dir / f"veloAE_{dataset}.h5ad"

    print(f"\n=== VeloAE: {dataset} ===", flush=True)
    print(f"input: {input_path}", flush=True)
    adata = sc.read(input_path, cache=False)
    adata.var_names_make_unique()

    missing_layers = {"spliced", "unspliced"} - set(adata.layers)
    if missing_layers:
        raise KeyError(f"{dataset}: missing required count layers {sorted(missing_layers)}")
    label_key = config["label_key"]
    embedding_key = config["embedding_key"]
    if label_key not in adata.obs:
        raise KeyError(f"{dataset}: missing label key {label_key!r}")
    if embedding_key not in adata.obsm:
        raise KeyError(f"{dataset}: missing embedding key {embedding_key!r}")
    adata.obs["clusters"] = adata.obs[label_key].copy()
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
            raise KeyError(f"{dataset}: missing moments {sorted(missing_moments)}")
        print("using established scVelo moments", flush=True)
    else:
        scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    args = get_parser().parse_args([f"--device=cuda:{GPU_INDEX}"])
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    args.output = str(output_dir)
    args.model_name = str(checkpoint_dir / f"veloAE_{dataset}.cpt")

    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(args.seed)
    np.random.seed(args.seed)
    torch.backends.cudnn.deterministic = True
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    print(f"VeloAE seed={args.seed}; device={device}", flush=True)

    # This supplies the raw-space stochastic velocity used as the auxiliary
    # supervision signal while fitting the VeloAE encoder.
    scv.tl.velocity(adata, vkey="stc_velocity", mode="stochastic")
    scv.tl.velocity_graph(adata, vkey="stc_velocity", n_jobs=args.scv_n_jobs)
    scv.tl.velocity_confidence(adata, vkey="stc_velocity")

    def tensor(matrix):
        values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
        return torch.FloatTensor(values).to(device)

    tensor_s = tensor(adata.layers["Ms"])
    tensor_u = tensor(adata.layers["Mu"])
    tensor_x = tensor(adata.X)
    tensor_v = tensor(adata.layers["stc_velocity"])

    model = init_model(adata, args, device)
    model = fit_model(
        args,
        adata,
        model,
        [tensor_s, tensor_u],
        tensor_v,
        [0, 1],
        device,
    )
    if not Path(args.model_name).is_file():
        raise RuntimeError(f"{dataset}: VeloAE did not save checkpoint {args.model_name}")

    model.load_state_dict(torch.load(args.model_name, map_location=device))
    model = model.to(device).eval()
    with torch.no_grad():
        x = model.encoder(tensor_x).cpu().numpy()
        s = model.encoder(tensor_s).cpu().numpy()
        u = model.encoder(tensor_u).cpu().numpy()
        v = estimate_ld_velocity(
            torch.from_numpy(s).to(device),
            torch.from_numpy(u).to(device),
            device=device,
            perc=[5, 95],
        ).cpu().numpy()

    adata = new_adata(
        adata,
        x,
        s,
        u,
        v,
        g_basis=args.ld_nb_g_src,
        n_nb_newadata=args.n_nb_newadata,
    )
    scv.tl.velocity_graph(adata, vkey="new_velocity", n_jobs=args.scv_n_jobs)
    scv.tl.velocity_confidence(adata, vkey="new_velocity")

    missing_layers = {"new_velocity"} - set(adata.layers)
    missing_obsm = {"X_umap", "new_velocity_umap"} - set(adata.obsm)
    missing_uns = {"new_velocity_graph", "new_velocity_graph_neg"} - set(adata.uns)
    if missing_layers or missing_obsm or missing_uns:
        raise RuntimeError(
            f"{dataset}: incomplete VeloAE result; "
            f"layers={sorted(missing_layers)}, obsm={sorted(missing_obsm)}, "
            f"uns={sorted(missing_uns)}"
        )

    adata.write_h5ad(output_path)
    print(f"wrote: {output_path}", flush=True)


def main():
    expected_versions = {
        "veloAE": "0.2.0",
        "scanpy": "1.9.8",
        "scvelo": "0.3.3",
        "torch": "2.4.1",
    }
    actual_versions = {name: version(name) for name in expected_versions}
    if actual_versions != expected_versions:
        raise RuntimeError(
            f"VeloAE environment mismatch: expected={expected_versions}, "
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
