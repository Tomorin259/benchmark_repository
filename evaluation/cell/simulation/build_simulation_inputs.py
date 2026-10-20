#!/usr/bin/env python3
"""Build the shared Cell reference state for the three historical simulations.

Run with the Cell preparation dependencies in docs/environments.md, including
scVelo 0.2.5. The large reference moments
are written outside Git; compact neighbor graphs and the input manifest live in
this module. No method-result or historical simulation H5AD is modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np
import scanpy as sc
import scvelo as scv
from scipy import sparse


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
METHODS = {
    "scVelo_dyn": "scvelo_dynamic",
    "MultiVelo": "multivelo",
    "UniTVelo": "unitvelo",
    "TIvelo": "TIvelo",
    "veloVI": "velovi",
    "cellDancer": "celldancer",
    "DeepVelo": "deepvelo",
    "DeepVelo_scv": "deepvelo2022",
    "Graphvelo_scv": "graphvelo_scv",
}
METHOD_ORDER = tuple(METHODS)
SCENARIOS = ("murk", "branching", "all")
EXPECTED_CELLS = {"murk": 1000, "branching": 2000, "all": 2000}
INPUT_COLUMNS = (
    "config_id", "operational_id", "dataset_id", "evaluation_dataset",
    "result_path", "velocity_layer_key", "source_embedding_key",
    "velocity_embedding_key", "alignment_mode", "n_common",
    "evaluation_role", "pred_unspliced",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve(artifact_root: Path, locator: str) -> Path:
    prefix = "velocity:"
    if not locator.startswith(prefix):
        raise ValueError(f"invalid velocity locator: {locator}")
    relative = Path(locator[len(prefix):])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"unsafe velocity locator: {locator}")
    path = artifact_root / "velocity" / relative
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def h5_names(group: h5py.Group) -> np.ndarray:
    key = group.attrs.get("_index", "_index")
    return np.asarray(group[key][:]).astype(str)


def compact_graph(graph: sparse.spmatrix, name: str) -> dict:
    csr = graph.tocsr()
    return {
        f"{name}_data": csr.data,
        f"{name}_indices": csr.indices,
        f"{name}_indptr": csr.indptr,
        f"{name}_shape": np.asarray(csr.shape, dtype=np.int64),
    }


def read_gene_sources() -> dict:
    path = REPO / "evaluation/gene/results/gene_cbdir_reference30_all.csv"
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    result = {}
    for row in rows:
        dataset, method = row["dataset"], row["method"]
        if dataset not in SCENARIOS or method not in METHODS:
            continue
        if row["status"] != "ok":
            raise ValueError(f"gene source not successful: {dataset}/{method}")
        key = (dataset, method)
        source = (row["source_h5ad"], row["velocity_s_key"])
        if key in result and result[key] != source:
            raise ValueError(f"inconsistent source for {dataset}/{method}")
        result[key] = source
    expected = {(dataset, method) for dataset in SCENARIOS for method in METHODS}
    if set(result) != expected:
        raise ValueError(f"missing gene sources: {sorted(expected - set(result))}")
    return result


def build_one(dataset: str, source: Path, source_locator: str,
              output_root: Path, moments_root: Path, artifact_root: Path,
              gene_sources: dict) -> dict:
    adata = sc.read_h5ad(source)
    if adata.n_obs != EXPECTED_CELLS[dataset]:
        raise ValueError(f"{dataset}: unexpected cell count {adata.n_obs}")
    if adata.X is None:
        adata.X = adata.layers["spliced"].copy()
        expression_mode = "spliced_if_X_missing"
    else:
        expression_mode = "native_X"
    scv.pp.filter_and_normalize(
        adata, min_shared_counts=20, n_top_genes=min(2000, adata.n_vars), log=True
    )
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.pp.neighbors(
        adata, n_neighbors=30, n_pcs=30, use_rep="X_pca",
        metric="euclidean", random_state=0,
    )
    names = np.asarray(adata.obs_names).astype(str)
    vars_ = np.asarray(adata.var_names).astype(str)
    indices = np.asarray(adata.uns["neighbors"]["indices"], dtype=np.int64)
    with np.load(
        REPO / f"evaluation/gene/results/reference_neighbors/{dataset}.npz",
        allow_pickle=False,
    ) as gene_ref:
        if not np.array_equal(names, gene_ref["obs_names"]):
            raise ValueError(f"{dataset}: Gene/Cell reference cell names differ")
        if not np.array_equal(indices, gene_ref["indices"]):
            raise ValueError(f"{dataset}: Gene/Cell reference 30-NN indices differ")
    for method in METHOD_ORDER:
        artifact, _ = gene_sources[dataset, method]
        method_source = resolve(
            artifact_root,
            artifact.replace("artifact://velocity/", "velocity:", 1),
        )
        with h5py.File(method_source) as handle:
            if not np.array_equal(names, h5_names(handle["obs"])):
                raise ValueError(f"{dataset}/{method}: cell order differs")
            if not np.array_equal(vars_, h5_names(handle["var"])):
                raise ValueError(f"{dataset}/{method}: gene order differs")
            if not np.allclose(
                adata.obsm["X_umap"], handle["obsm"]["X_umap"][:],
                rtol=0, atol=1e-5,
            ):
                raise ValueError(f"{dataset}/{method}: UMAP differs")
    source_hash = sha256(source)
    neighbor_metadata = {
        "dataset": dataset,
        "source": source_locator,
        "source_sha256": source_hash,
        "scvelo_version": scv.__version__,
        "n_obs": int(adata.n_obs),
        "n_vars_after_filter": int(adata.n_vars),
        "parameters": {
            "min_shared_counts": 20, "n_top_genes_requested": min(2000, adata.n_vars),
            "log": True, "n_neighbors": 30, "n_pcs": 30,
            "use_rep": "X_pca", "metric": "euclidean", "random_state": 0,
        },
        "gene_reference_indices_exact": True,
    }
    neighbor_path = output_root / "reference_neighbors" / f"{dataset}.npz"
    neighbor_path.parent.mkdir(parents=True, exist_ok=True)
    neighbor_temp = neighbor_path.with_suffix(".tmp.npz")
    np.savez_compressed(
        neighbor_temp,
        obs_names=names,
        indices=indices,
        reference_pca=np.asarray(adata.obsm["X_pca"], dtype=np.float32),
        reference_var_names=vars_,
        metadata=np.asarray(json.dumps(neighbor_metadata, sort_keys=True)),
        **compact_graph(adata.obsp["distances"], "distances"),
        **compact_graph(adata.obsp["connectivities"], "connectivities"),
    )
    os.replace(neighbor_temp, neighbor_path)

    moments_root.mkdir(parents=True, exist_ok=True)
    moments_path = moments_root / f"{dataset}.npz"
    moments_temp = moments_path.with_suffix(".tmp.npz")
    moments_metadata = {
        "dataset": dataset,
        "source": source_locator,
        "source_sha256": source_hash,
        "expression_mode": expression_mode,
        "preprocessing": "scvelo_0.2.5/filter_and_normalize(min_shared_counts=20,n_top_genes=min(2000,n_vars),log=True)",
        "smoothing_graph": "simulation/reference_neighbors/fixed_30nn",
    }
    np.savez_compressed(
        moments_temp,
        obs_names=names, var_names=vars_,
        Ms=np.asarray(adata.layers["Ms"], dtype=np.float32),
        Mu=np.asarray(adata.layers["Mu"], dtype=np.float32),
        metadata=np.asarray(json.dumps(moments_metadata, sort_keys=True)),
    )
    os.replace(moments_temp, moments_path)
    return {
        "source": source_locator,
        "source_sha256": source_hash,
        "reference_neighbors_sha256": sha256(neighbor_path),
        "reference_moments_sha256": sha256(moments_path),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--reference-moments-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=HERE)
    args = parser.parse_args()
    config = json.loads((HERE / "config.json").read_text(encoding="utf-8"))
    core = json.loads((REPO / "evaluation/cell/config.json").read_text(encoding="utf-8"))
    if config["reference"] != core["reference"]:
        raise ValueError("simulation and canonical Cell reference contracts differ")
    gene_sources = read_gene_sources()
    registry = {
        row["operational_id"]: row["config_id"]
        for row in csv.DictReader(
            (REPO / "methods/configurations.tsv").open(encoding="utf-8"),
            delimiter="\t",
        )
    }
    rows = []
    manifest = {"scvelo_version": scv.__version__, "scenarios": {}}
    for dataset in SCENARIOS:
        dataset_config = config["datasets"][dataset]
        source_locator = dataset_config["reference_umap_path"]
        source = resolve(args.artifact_root, source_locator)
        manifest["scenarios"][dataset] = build_one(
            dataset, source, source_locator, args.output_root,
            args.reference_moments_root, args.artifact_root, gene_sources,
        )
        for method in METHOD_ORDER:
            external_locator, layer_key = gene_sources[dataset, method]
            operational = METHODS[method]
            rows.append({
                "config_id": registry[operational],
                "operational_id": operational,
                "dataset_id": dataset_config["dataset_id"],
                "evaluation_dataset": dataset,
                "result_path": external_locator.replace(
                    "artifact://velocity/", "velocity:", 1
                ),
                "velocity_layer_key": layer_key,
                "source_embedding_key": "X_umap",
                "velocity_embedding_key": "velocity_umap",
                "alignment_mode": (
                    "reference_graph_projection" if method == "DeepVelo" else "direct"
                ),
                "n_common": str(EXPECTED_CELLS[dataset]),
                "evaluation_role": "simulation_gene_matched_cell",
                "pred_unspliced": "true" if method == "DeepVelo" else "not_applicable",
            })
    args.output_root.mkdir(parents=True, exist_ok=True)
    with (args.output_root / "inputs.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.output_root / "reference_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Built {len(rows)} inputs for {len(SCENARIOS)} scenarios")


if __name__ == "__main__":
    main()
