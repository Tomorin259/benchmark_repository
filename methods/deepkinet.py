#!/usr/bin/env python3
"""Reproduce C13 (DeepKINET) with the upstream raw-count workflow.

DeepKINET requires raw spliced/unspliced counts even after gene selection.
For the four count-sensitive benchmark datasets, the established working
cell/gene panel is retained and aligned raw counts are substituted without
rounding.  The other ten datasets follow the official tutorial's single
scVelo ``filter_and_normalize`` call before the same raw-count substitution.

The official tutorial does not prescribe a global seed.  The published
benchmark result is frozen with seed 123 and deterministic PyTorch/CUDA
settings so that dataset order cannot change the formal result.
"""

import json
import os
import random
import sys
from pathlib import Path

# These variables must be set before importing torch or creating a CUDA context.
BENCHMARK_SEED = 123
OFFICIAL_NUMPY_SPLIT_SEED = 42
os.environ.setdefault("PYTHONHASHSEED", str(BENCHMARK_SEED))
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import scanpy as sc
import scvelo as scv
import torch
from scipy import sparse

import deepkinet


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPOSITORY_ROOT.parent / "data"
OUTPUT_ROOT = REPOSITORY_ROOT / "outputs" / "deepkinet"
CHECKPOINT_ROOT = OUTPUT_ROOT / "checkpoints"
RNG_ROOT = OUTPUT_ROOT / "rng_states"
MANIFEST_ROOT = OUTPUT_ROOT / "manifests"

# Edit this list for a focused rerun. The formal result uses every listed
# dataset and resets the benchmark seed before each independent dataset.
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

# These are the benchmark's fixed input contracts. A separate ``counts`` path
# means that the working object's cells/genes are retained while both raw count
# layers are aligned from that native count artifact.
DATASETS = {
    "bonemarrow": {"input": "bonemarrow/bonemarrow.h5ad", "counts": None},
    "capo_pre": {
        "input": "capo_pre/capo_pre.h5ad",
        "counts": "capo_pre/capo_pre_counts.h5ad",
        "expected_shape": (1222, 2000),
    },
    "DentateGyrus_lamanno": {
        "input": "DentateGyrus_lamanno/DentateGyrus_lamanno.h5ad",
        "counts": None,
    },
    "ForebrainGlut": {"input": "ForebrainGlut/ForebrainGlut.loom", "counts": None},
    "Gastrulation": {"input": "Gastrulation/Gastrulation.h5ad", "counts": None},
    "Pancreas": {"input": "Pancreas/Pancreas.h5ad", "counts": None},
    "dentategyrus": {"input": "dentategyrus/dentategyrus.h5ad", "counts": None},
    "hematopoiesis_dynamo": {
        "input": "hematopoiesis_dynamo/hematopoiesis_dynamo.h5ad",
        "counts": None,
    },
    "mousebrain": {
        "input": "mousebrain/mousebrain.h5ad",
        "counts": "mousebrain/mousebrain_counts.h5ad",
        "expected_shape": (3365, 936),
    },
    "neuron": {
        "input": "neuron/neuron.h5ad",
        "counts": "neuron/neuron_counts.h5ad",
        "expected_shape": (3066, 2000),
    },
    "organoid": {"input": "organoid/organoid_pre.h5ad", "counts": None},
    "retina": {"input": "retina/retina_full.h5ad", "counts": None},
    "rpe1": {"input": "rpe1/rpe1.h5ad", "counts": None},
    "rpe1_gen": {
        "input": "rpe1_gen/rpe1_gen.h5ad",
        "counts": "rpe1_gen/rpe1_gen_counts.h5ad",
        "expected_shape": (3354, 2000),
    },
}

MIN_SHARED_COUNTS = 20
N_TOP_GENES = 2000
EPOCHS = 2000
LOSS_MODE = "poisson"


def seed_everything(seed):
    """Reset all user-controlled RNGs before one dataset run."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    if hasattr(torch.backends.cuda, "matmul"):
        torch.backends.cuda.matmul.allow_tf32 = False
    if hasattr(torch.backends.cudnn, "allow_tf32"):
        torch.backends.cudnn.allow_tf32 = False
    torch.use_deterministic_algorithms(True)


def capture_rng_state():
    return {
        "python": random.getstate(),
        "numpy": np.random.get_state(),
        "torch_cpu": torch.get_rng_state(),
        "torch_cuda": torch.cuda.get_rng_state_all() if torch.cuda.is_available() else [],
    }


def read_data(path):
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".loom":
        return sc.read_loom(path, sparse=True)
    return sc.read_h5ad(path)


def lossless_int64(matrix, label):
    """Print and enforce DeepKINET's raw integer-count input invariant."""
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).reshape(-1)
    values = np.asarray(values)
    finite = bool(np.isfinite(values).all())
    nonnegative = finite and bool((values >= 0).all())
    integer_valued = finite and bool(np.allclose(values, np.rint(values), atol=1e-8, rtol=0))
    maximum_fraction = (
        float(np.abs(values.astype(np.float64) - np.rint(values)).max()) if values.size else 0.0
    )
    print(
        f"[counts check] {label}: dtype={matrix.dtype} shape={matrix.shape} "
        f"finite={finite} nonnegative={nonnegative} integer_valued={integer_valued} "
        f"max_abs_fractional_deviation={maximum_fraction:.3g}",
        flush=True,
    )
    if not (nonnegative and integer_valued):
        raise ValueError(f"{label}: refusing to round or cast non-integer raw counts")

    converted = matrix.astype(np.int64) if sparse.issparse(matrix) else np.asarray(matrix, dtype=np.int64)
    if sparse.issparse(matrix):
        difference = matrix.tocsr().astype(np.float64) - converted.tocsr().astype(np.float64)
        difference.eliminate_zeros()
        preserved = difference.nnz == 0
    else:
        preserved = np.array_equal(np.asarray(matrix), converted)
    if not preserved:
        raise AssertionError(f"{label}: values changed during int64 dtype conversion")
    print(f"[counts check] {label}: lossless int64 dtype conversion verified", flush=True)
    return converted


def restore_raw_counts(adata, raw_adata, dataset):
    """Align native raw counts to the current benchmark cells and gene panel."""
    if adata.obs_names.has_duplicates or raw_adata.obs_names.has_duplicates:
        raise ValueError(f"{dataset}: obs_names must be unique before count alignment")
    if adata.var_names.has_duplicates or raw_adata.var_names.has_duplicates:
        raise ValueError(f"{dataset}: var_names must be unique before count alignment")

    cells = raw_adata.obs_names.get_indexer(adata.obs_names)
    genes = raw_adata.var_names.get_indexer(adata.var_names)
    if (cells < 0).any():
        raise KeyError(f"{dataset}: raw counts missing working cells {adata.obs_names[cells < 0][:10].tolist()}")
    if (genes < 0).any():
        raise KeyError(f"{dataset}: raw counts missing selected genes {adata.var_names[genes < 0][:10].tolist()}")

    raw_subset = raw_adata[cells, genes]
    layer_aliases = {
        "spliced": ("spliced", "Spliced", "X_spliced"),
        "unspliced": ("unspliced", "Unspliced", "X_unspliced"),
    }
    for target_layer, source_layers in layer_aliases.items():
        source_layer = next((key for key in source_layers if key in raw_subset.layers), None)
        if source_layer is None:
            raise KeyError(f"{dataset}: raw count file has no {target_layer!r} layer")
        if source_layer != target_layer:
            print(f"layer mapping: {source_layer!r} -> {target_layer!r}", flush=True)
        adata.layers[target_layer] = lossless_int64(
            raw_subset.layers[source_layer],
            f"{dataset}.layers['{target_layer}']",
        )


def run_dataset(dataset):
    config = DATASETS[dataset]
    input_path = DATA_ROOT / config["input"]
    counts_path = DATA_ROOT / config["counts"] if config["counts"] else input_path
    output_path = OUTPUT_ROOT / f"DeepKINET_{dataset}.h5ad"
    checkpoint_path = CHECKPOINT_ROOT / f"{dataset}.deepkinet_opt.pt"
    manifest_path = MANIFEST_ROOT / f"{dataset}.json"
    initial_rng_path = RNG_ROOT / f"{dataset}.initial.pt"
    final_rng_path = RNG_ROOT / f"{dataset}.final.pt"

    seed_everything(BENCHMARK_SEED)
    CHECKPOINT_ROOT.mkdir(parents=True, exist_ok=True)
    RNG_ROOT.mkdir(parents=True, exist_ok=True)
    MANIFEST_ROOT.mkdir(parents=True, exist_ok=True)
    manifest = {
        "seed": BENCHMARK_SEED,
        "numpy_split_seed_in_official_package": OFFICIAL_NUMPY_SPLIT_SEED,
        "torch_deterministic_algorithms": True,
        "cudnn_deterministic": True,
        "cudnn_benchmark": False,
        "data_loader_num_workers": 2,
        "data_loader_shuffle": True,
        "epoch": EPOCHS,
        "loss_mode": LOSS_MODE,
        "batch_key": None,
        "dataset": dataset,
        "working_input": str(input_path.resolve()),
        "raw_counts_input": str(counts_path.resolve()),
        "checkpoint": str(checkpoint_path.resolve()),
        "output": str(output_path.resolve()),
        "initial_rng_state": str(initial_rng_path.resolve()),
        "final_rng_state": str(final_rng_path.resolve()),
        "status": "started",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    torch.save(capture_rng_state(), initial_rng_path)

    print(f"\n=== DeepKINET: {dataset} (seed={BENCHMARK_SEED}) ===", flush=True)
    print(f"working input: {input_path}", flush=True)
    print(f"raw counts: {counts_path}", flush=True)
    adata = read_data(input_path)
    raw_adata = read_data(counts_path)
    adata.var_names_make_unique()
    raw_adata.var_names_make_unique()

    if config["counts"] is None:
        # Official tutorial: select the benchmark gene panel first, then put
        # raw count layers back before calling the model.
        scv.pp.filter_and_normalize(
            adata,
            min_shared_counts=MIN_SHARED_COUNTS,
            n_top_genes=N_TOP_GENES,
        )
    else:
        if adata.shape != config["expected_shape"]:
            raise ValueError(
                f"{dataset}: expected fixed benchmark shape {config['expected_shape']}, "
                f"found {adata.shape}"
            )
        print("scVelo preprocessing: skipped for fixed benchmark working object", flush=True)

    restore_raw_counts(adata, raw_adata, dataset)
    print(f"model input: {adata.shape}", flush=True)

    # Official DeepKINET call. The checkpoint location is explicit so a run
    # cannot overwrite another dataset's training state or the repository root.
    adata, _ = deepkinet.estimate_kinetics(
        adata,
        epoch=EPOCHS,
        loss_mode=LOSS_MODE,
        checkpoint=str(checkpoint_path),
        batch_key=None,
    )

    required_layers = {
        "spliced",
        "unspliced",
        "DeepKINET_velocity",
        "splicing_rate",
        "degradation_rate",
    }
    required_obsm = {"latent_variable", "latent_velocity"}
    missing_layers = required_layers - set(adata.layers)
    missing_obsm = required_obsm - set(adata.obsm)
    if missing_layers or missing_obsm:
        raise RuntimeError(
            f"{dataset}: incomplete DeepKINET result; "
            f"missing layers={sorted(missing_layers)}, obsm={sorted(missing_obsm)}"
        )

    adata.uns["DeepKINET_reproducibility"] = {
        "seed": BENCHMARK_SEED,
        "numpy_split_seed_in_official_package": OFFICIAL_NUMPY_SPLIT_SEED,
        "torch_deterministic_algorithms": True,
        "cudnn_deterministic": True,
        "cudnn_benchmark": False,
        "data_loader_num_workers": 2,
        "data_loader_shuffle": True,
        "runner": "benchmark_repository/methods/deepkinet.py",
    }
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(".tmp.h5ad")
    adata.write_h5ad(temporary_path)
    temporary_path.replace(output_path)
    torch.save(capture_rng_state(), final_rng_path)
    manifest["status"] = "completed"
    manifest["shape"] = list(adata.shape)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"wrote: {output_path}", flush=True)


def main():
    try:
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass
    if scv.__version__ != "0.2.5":
        raise RuntimeError(f"scVelo 0.2.5 is required, found {scv.__version__}")
    unknown_datasets = set(DATASETS_TO_RUN) - set(DATASETS)
    if unknown_datasets:
        raise ValueError(f"unknown datasets: {sorted(unknown_datasets)}")
    print(f"DeepKINET source: {deepkinet.__file__}", flush=True)
    print(f"benchmark RNG seed: {BENCHMARK_SEED}", flush=True)

    for dataset in DATASETS_TO_RUN:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
