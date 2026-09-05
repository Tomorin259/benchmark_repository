#!/usr/bin/env python3
"""Evaluate cell-level trajectory metrics from explicit result inputs.

The evaluator consumes a relative result-input table and a dataset contract.
It keeps method-native velocity layers at the input boundary, computes any
missing embedding required by the contract, and writes one compact row per
method/dataset pair.  It does not train methods or silently select a nearby
result file.
"""

from __future__ import annotations

import argparse
import io
import json
import re
from contextlib import redirect_stdout
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MODULE_DIR / "config.json"
DEFAULT_INPUTS = MODULE_DIR / "result_inputs.tsv"


def load_config(path: Path = DEFAULT_CONFIG) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def resolve_artifact(artifact_root: Path, relative_path: str) -> Path:
    path = Path(relative_path)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Result paths must be relative to artifact root: {relative_path}")
    candidate = artifact_root / path
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    # External artifact roots are often assembled from symlinks to a mounted
    # archive. Lexical traversal is rejected above; a controlled symlink is a
    # valid artifact-root boundary and should not be rejected after resolving.
    return candidate.resolve()


def read_result_inputs(path: Path = DEFAULT_INPUTS) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {
        "config_id",
        "operational_id",
        "dataset_id",
        "dataset",
        "result_path",
        "velocity_layer_key",
        "source_embedding_key",
        "velocity_embedding_key",
        "evaluation_embedding",
        "alignment_mode",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if frame.duplicated(["operational_id", "dataset"]).any():
        raise ValueError(f"{path}: duplicate method/dataset input")
    return frame


def dense(values):
    if sparse.issparse(values):
        return values.toarray()
    return np.asarray(values)


def prepare_labels(adata: ad.AnnData, dataset: str, config: dict) -> str:
    dataset_config = config["datasets"][dataset]
    target = dataset_config["label"]
    label_config = config["labels"]
    source = label_config.get("dataset_specific_sources", {}).get(dataset)
    if source and source in adata.obs:
        adata.obs[target] = adata.obs[source].copy()
    elif target not in adata.obs:
        for candidate in label_config.get("fallback_sources", []):
            if candidate in adata.obs:
                adata.obs[target] = adata.obs[candidate].copy()
                break
    if target not in adata.obs:
        raise KeyError(f"{dataset}: label column {target!r} is missing")
    return target


def prepare_velocity(adata: ad.AnnData, row: pd.Series) -> None:
    layer = row["velocity_layer_key"]
    if layer not in adata.layers:
        raise KeyError(f"{row['operational_id']}: missing layers/{layer}")
    values = dense(adata.layers[layer]).astype(np.float64, copy=False)
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    adata.layers["velocity"] = np.clip(values, -1e10, 1e10)


def _derive_neighbor_indices(adata: ad.AnnData) -> None:
    if "indices" in adata.uns.get("neighbors", {}):
        return
    graph = adata.obsp.get("distances")
    if graph is None:
        graph = adata.obsp.get("connectivities")
    if graph is None or not sparse.issparse(graph):
        raise KeyError("neighbors['indices'] is missing and no sparse neighbor graph exists")
    graph = graph.tocsr()
    indices = np.empty(adata.n_obs, dtype=object)
    for index in range(adata.n_obs):
        start, end = graph.indptr[index], graph.indptr[index + 1]
        indices[index] = graph.indices[start:end].astype(np.int64, copy=False)
    adata.uns.setdefault("neighbors", {})["indices"] = indices


def ensure_neighbors(adata: ad.AnnData, scv, n_neighbors: int) -> None:
    if "neighbors" not in adata.uns:
        if adata.X is None and "spliced" in adata.layers:
            adata.X = adata.layers["spliced"].copy()
        scv.pp.moments(
            adata,
            n_pcs=min(30, adata.n_vars),
            n_neighbors=n_neighbors,
        )
        scv.pp.neighbors(adata)
    elif "indices" not in adata.uns["neighbors"]:
        adata.uns.pop("neighbors", None)
        adata.obsp.pop("connectivities", None)
        adata.obsp.pop("distances", None)
        scv.pp.neighbors(adata)
    _derive_neighbor_indices(adata)


def ensure_velocity_graph(adata: ad.AnnData, scv, n_jobs: int = 8) -> None:
    if "velocity_graph" in adata.uns:
        return
    for key in ("Velocity_graph", "velo_s_norm_graph", "sde_velocity_graph"):
        if key in adata.uns:
            adata.uns["velocity_graph"] = adata.uns[key].copy()
            return
    if "pearson_transition_matrix" in adata.obsp:
        adata.uns["velocity_graph"] = adata.obsp["pearson_transition_matrix"].copy()
        return
    scv.tl.velocity_graph(adata, vkey="velocity", n_jobs=n_jobs)


def prepare_embedding(
    adata: ad.AnnData,
    row: pd.Series,
    scv,
    n_neighbors: int,
    n_jobs: int,
) -> str:
    embedding = row["evaluation_embedding"]
    source_key = row["source_embedding_key"]
    target_key = f"X_{embedding}"
    if target_key not in adata.obsm:
        if source_key not in adata.obsm:
            raise KeyError(f"{row['operational_id']}: missing obsm/{source_key}")
        adata.obsm[target_key] = np.asarray(adata.obsm[source_key]).copy()
    velocity_embedding_key = row["velocity_embedding_key"]
    runtime_computed = velocity_embedding_key.startswith("computed:")
    if not runtime_computed and velocity_embedding_key in adata.obsm:
        adata.obsm[f"velocity_{embedding}"] = np.asarray(
            adata.obsm[velocity_embedding_key]
        ).copy()
        return embedding
    ensure_neighbors(adata, scv, n_neighbors)
    if row["operational_id"] == "DeepKINET":
        # The historical evaluator explicitly reran neighbors for this result
        # family before constructing the missing velocity embedding.
        scv.pp.neighbors(adata)
    ensure_velocity_graph(adata, scv, n_jobs=n_jobs)
    scv.tl.velocity_embedding(
        adata, basis=embedding, vkey="velocity", all_comps=True
    )
    return embedding


def consistency_scores(
    velocity: np.ndarray, neighbor_indices, n_neighbors: int
) -> np.ndarray:
    """Return the mean cosine similarity to each cell's nearest neighbours."""
    velocity = np.asarray(velocity, dtype=np.float64)
    scores = np.zeros(velocity.shape[0], dtype=float)
    for index in range(velocity.shape[0]):
        neighbors = np.asarray(neighbor_indices[index])[1 : n_neighbors + 1]
        if neighbors.size == 0:
            continue
        current = velocity[index]
        neighbor_values = velocity[neighbors]
        valid = np.isfinite(current) & np.isfinite(neighbor_values).all(axis=0)
        if not valid.any():
            continue
        current = current[valid]
        neighbor_values = neighbor_values[:, valid]
        current_norm = np.linalg.norm(current)
        neighbor_norms = np.linalg.norm(neighbor_values, axis=1)
        denominators = current_norm * neighbor_norms
        valid_norms = denominators > 0
        if valid_norms.any():
            scores[index] = np.mean(
                neighbor_values[valid_norms].dot(current) / denominators[valid_norms]
            )
    return scores


def parse_unitvelo_means(output: str) -> tuple[float, float]:
    cbdir = re.search(
        r"# Cross-Boundary Direction Correctness \(A->B\).*?Total Mean:\s*([^\n\r]+)",
        output,
        flags=re.DOTALL,
    )
    iccoh = re.search(
        r"# In-cluster Coherence.*?Total Mean:\s*([^\n\r]+)",
        output,
        flags=re.DOTALL,
    )
    if cbdir is None or iccoh is None:
        raise ValueError("unitvelo.evaluate output did not contain both Total Mean values")
    return float(cbdir.group(1).strip()), float(iccoh.group(1).strip())


def evaluate_one(
    artifact_root: Path,
    row: pd.Series,
    config: dict,
    dependencies: dict,
    progress: bool = False,
) -> dict:
    scv = dependencies["scv"]
    unitvelo = dependencies["unitvelo"]
    compute_batc = dependencies["compute_batc"]
    dataset = row["dataset"]
    dataset_config = config["datasets"][dataset]
    path = resolve_artifact(artifact_root, row["result_path"])
    adata = ad.read_h5ad(path)
    prepare_velocity(adata, row)
    label = prepare_labels(adata, dataset, config)
    ensure_neighbors(adata, scv, config["metrics"]["n_neighbors"])
    embedding = prepare_embedding(
        adata,
        row,
        scv,
        config["metrics"]["n_neighbors"],
        config["metrics"]["velocity_graph_n_jobs"],
    )
    neighbors = adata.uns["neighbors"]["indices"]
    consistency = consistency_scores(
        dense(adata.layers["velocity"]),
        neighbors,
        config["metrics"]["n_neighbors"],
    )

    evaluation_adata = adata
    if row["operational_id"].lower() == "dynamo" and "velocity_genes" in adata.var:
        evaluation_adata = adata[:, adata.var["velocity_genes"].astype(bool)].copy()
    output = io.StringIO()
    with redirect_stdout(output):
        unitvelo_result = unitvelo.evaluate(
            evaluation_adata,
            dataset_config["cluster_edges"],
            label,
            k_velocity="velocity",
            x_emb=f"X_{embedding}",
        )
    cbdir, iccoh = parse_unitvelo_means(output.getvalue())

    batc_config = dataset_config.get("batc", {})
    batc_embedding = batc_config.get("embedding", embedding)
    batc_velocity_key = "velocity"
    if batc_config.get("recompute_embedding"):
        scv.pp.neighbors(adata, n_neighbors=config["metrics"]["n_neighbors"], n_pcs=min(30, adata.n_vars))
        scv.tl.umap(adata)
        scv.tl.velocity_graph(adata, vkey="velocity", n_jobs=config["metrics"]["velocity_graph_n_jobs"])
        scv.tl.velocity_embedding(adata, basis=batc_embedding, vkey="velocity", all_comps=True)
        batc_velocity_key = f"velocity_{batc_embedding}"
    batc_score = compute_batc(
        adata,
        dataset=dataset,
        cluster_edges=dataset_config["cluster_edges"],
        cluster_key=label,
        n_bins=dataset_config["batc_n_bins"],
        embedding_key=batc_embedding,
        velocity_key=batc_velocity_key,
        progress=progress,
    )
    if isinstance(batc_score, (list, tuple)):
        batc_score = float(np.mean(batc_score))
    return {
        "dataset": dataset,
        "method": row["operational_id"],
        "consistency": float(np.mean(consistency)),
        "cbdir": float(cbdir),
        "iccoh": float(iccoh),
        "batc": float(batc_score),
        "result_path": row["result_path"],
        "velocity_layer_key": row["velocity_layer_key"],
        "unitvelo_raw": unitvelo_result,
    }


def load_dependencies():
    import scvelo as scv
    import unitvelo
    try:
        from batc import compute_batc
    except ModuleNotFoundError as error:
        if error.name != "batc":
            raise
        from .batc import compute_batc

    return {"scv": scv, "unitvelo": unitvelo, "compute_batc": compute_batc}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--inputs", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--methods", nargs="*")
    parser.add_argument("--datasets", nargs="*")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.check_only == (args.output is not None):
        parser.error("use exactly one of --output or --check-only")

    config = load_config(args.config)
    rows = read_result_inputs(args.inputs)
    if args.methods:
        rows = rows[rows["operational_id"].isin(args.methods)]
    if args.datasets:
        rows = rows[rows["dataset"].isin(args.datasets)]
    if rows.empty:
        raise ValueError("No result inputs selected")
    artifact_root = args.artifact_root.resolve()
    if args.check_only:
        for _, row in rows.iterrows():
            resolve_artifact(artifact_root, row["result_path"])
        print(f"Validated {len(rows)} result inputs")
        return

    dependencies = load_dependencies()
    results = []
    for _, row in rows.iterrows():
        print(f"[{row['dataset']}] {row['operational_id']}", flush=True)
        results.append(
            evaluate_one(
                artifact_root,
                row,
                config,
                dependencies,
                progress=args.progress,
            )
        )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = ["dataset", "method", "consistency", "cbdir", "iccoh", "batc"]
    pd.DataFrame(results).loc[:, columns].to_csv(output, index=False)
    print(f"Wrote {len(results)} cell-level rows to {output}")


if __name__ == "__main__":
    main()
