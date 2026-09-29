#!/usr/bin/env python3
"""Prepare one audited cell-metric input with fixed 30-NN evaluation state.

This program runs only in the indicators environment (scvelo 0.2.5).  It
never changes a method result: it writes a separate temporary AnnData object
that contains the configured velocity layer, canonical neighbours, a canonical
UMAP, and a freshly computed scvelo velocity graph.  It preserves each
method's stored 2-D velocity projection except for explicitly registered
coordinate-incompatible method/dataset pairs.
The paired TiVelo program uses the temporary object only for its two official
metrics.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
from contextlib import redirect_stdout
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics.pairwise import cosine_similarity


HERE = Path(__file__).resolve().parent
CELL_EVALUATION_ROOT = HERE
sys.path.insert(0, str(CELL_EVALUATION_ROOT))

from batc import compute_batc  # noqa: E402
from metrics import (  # noqa: E402
    consistency_scores,
    dense,
    load_config,
    parse_unitvelo_means,
    prepare_labels,
)
from reference_neighbors import (  # noqa: E402
    DEFAULT_REFERENCE_NEIGHBORS_ROOT,
    REFERENCE_PARAMS,
    apply_reference_neighbors,
    load_reference_neighbors,
)


def read_inputs(path: Path) -> pd.DataFrame:
    inputs = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {
        "operational_id",
        "evaluation_dataset",
        "result_path",
        "velocity_layer_key",
    }
    missing = required.difference(inputs.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    return inputs


def resolve_source(artifact_root: Path, source_path: str) -> Path:
    prefix, separator, relative = source_path.partition(":")
    if prefix != "velocity" or not separator:
        raise ValueError(f"unsupported source result path: {source_path}")
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise ValueError(f"unsafe source result path: {source_path}")
    candidate = artifact_root / "velocity" / relative_path
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    return candidate


def select_row(inputs: pd.DataFrame, method: str, dataset: str) -> pd.Series:
    matched = inputs[
        (inputs["operational_id"] == method)
        & (inputs["evaluation_dataset"] == dataset)
    ]
    if len(matched) != 1:
        raise ValueError(
            f"expected exactly one input row for {method}/{dataset}, found {len(matched)}"
        )
    return matched.iloc[0]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_npz_atomic(path: Path, **payload) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp.npz")
    np.savez_compressed(temporary, **payload)
    os.replace(temporary, path)


def prepare_velocity(adata: ad.AnnData, layer_key: str) -> int:
    if layer_key not in adata.layers:
        raise KeyError(f"missing layers/{layer_key!r}")
    values = dense(adata.layers[layer_key]).astype(np.float64, copy=False)
    invalid = int(np.size(values) - np.count_nonzero(np.isfinite(values)))
    # This preserves the historical evaluator's explicit invalid-value policy.
    adata.layers["velocity"] = np.nan_to_num(
        values, nan=0.0, posinf=0.0, neginf=0.0
    )
    return invalid


def canonical_umap(
    adata: ad.AnnData,
    dataset: str,
    inputs: pd.DataFrame,
    artifact_root: Path,
    cache_root: Path,
    reference_neighbors_root: Path,
    config: dict,
) -> dict:
    """Load one registered dataset-level UMAP, not a method-specific embedding."""
    cache_root.mkdir(parents=True, exist_ok=True)
    cache_path = cache_root / f"{dataset}.npz"
    reference_locator = config["datasets"][dataset].get("reference_umap_path")
    if reference_locator is None:
        reference = select_row(inputs, "scvelo_deterministic", dataset)
        reference_locator = reference["result_path"]
    reference_path = resolve_source(artifact_root, reference_locator)
    reference_neighbor_path = reference_neighbors_root / f"{dataset}.npz"
    reference_neighbor_sha = sha256(reference_neighbor_path)

    cache_matches = False
    if cache_path.is_file():
        with np.load(cache_path, allow_pickle=False) as cached:
            cached_source = (
                str(cached["source_path"]) if "source_path" in cached.files else ""
            )
            cached_source_sha = (
                str(cached["source_sha256"])
                if "source_sha256" in cached.files
                else ""
            )
            cached_neighbor_sha = (
                str(cached["reference_neighbor_sha256"])
                if "reference_neighbor_sha256" in cached.files
                else ""
            )
            construction = (
                str(cached["construction"])
                if "construction" in cached.files
                else ""
            )
            cache_matches = (
                cached_source == reference_locator
                and cached_source_sha == sha256(reference_path)
            ) or (
                construction == "scanpy_umap_random_state_0_on_fixed_reference_graph"
                and cached_neighbor_sha == reference_neighbor_sha
            )

    if not cache_matches:
        reference_adata = ad.read_h5ad(reference_path, backed="r")
        try:
            if "X_umap" not in reference_adata.obsm:
                raise KeyError(f"{reference_path}: missing obsm/'X_umap'")
            reference_umap = np.asarray(reference_adata.obsm["X_umap"], dtype=np.float64)
            if reference_umap.ndim == 2 and reference_umap.shape[1] == 2:
                save_npz_atomic(
                    cache_path,
                    obs_names=np.asarray(reference_adata.obs_names).astype(str),
                    X_umap=reference_umap,
                    source_path=np.asarray(reference_locator),
                    source_sha256=np.asarray(sha256(reference_path)),
                    reference_neighbor_sha256=np.asarray(""),
                    construction=np.asarray("source_X_umap"),
                )
        finally:
            reference_adata.file.close()

    # Forebrain's historical X_umap is actually 30-D PCA.  Generate one
    # deterministic dataset-level 2-D representation from the audited
    # reference PCA and graph, never from a method-specific output.
    if not cache_path.is_file():
        import scanpy as sc
        from reference_neighbors import load_reference_neighbors
        reference_graph = load_reference_neighbors(
            dataset, output_root=reference_neighbors_root
        )
        with np.load(
            reference_neighbors_root / f"{dataset}.npz", allow_pickle=False
        ) as reference_asset:
            reference_pca = np.asarray(reference_asset["reference_pca"], dtype=np.float32)
        umap_adata = ad.AnnData(
            X=np.zeros((len(reference_graph["obs_names"]), 1), dtype=np.float32)
        )
        umap_adata.obs_names = reference_graph["obs_names"]
        umap_adata.obsm["X_pca"] = reference_pca
        umap_adata.obsp["distances"] = reference_graph["distances"]
        umap_adata.obsp["connectivities"] = reference_graph["connectivities"]
        umap_adata.uns["neighbors"] = {
            "indices": reference_graph["indices"],
            "distances_key": "distances",
            "connectivities_key": "connectivities",
            "params": {
                "n_neighbors": REFERENCE_PARAMS["n_neighbors"],
                "n_pcs": REFERENCE_PARAMS["n_pcs"],
                "use_rep": REFERENCE_PARAMS["use_rep"],
                "method": "umap",
            },
        }
        sc.tl.umap(umap_adata, random_state=0)
        coordinates = np.asarray(umap_adata.obsm["X_umap"], dtype=np.float64)
        if coordinates.ndim != 2 or coordinates.shape[1] != 2:
            raise ValueError(f"{dataset}: generated embedding is not 2-D")
        save_npz_atomic(
            cache_path,
            obs_names=np.asarray(reference_graph["obs_names"]).astype(str),
            X_umap=coordinates,
            source_path=np.asarray("generated_from_fixed_reference_graph"),
            source_sha256=np.asarray(""),
            reference_neighbor_sha256=np.asarray(reference_neighbor_sha),
            construction=np.asarray("scanpy_umap_random_state_0_on_fixed_reference_graph"),
        )

    with np.load(cache_path, allow_pickle=False) as cached:
        names = np.asarray(cached["obs_names"]).astype(str)
        coordinates = np.asarray(cached["X_umap"], dtype=np.float64)
        if coordinates.ndim != 2 or coordinates.shape[1] != 2:
            raise ValueError(
                f"{dataset}: cached evaluation embedding must be 2-D, got {coordinates.shape}"
            )
        provenance = {
            "reference_umap_source": str(cached["source_path"]),
            "reference_umap_sha256": str(cached["source_sha256"]),
            "reference_neighbor_sha256": (
                str(cached["reference_neighbor_sha256"])
                if "reference_neighbor_sha256" in cached.files
                else ""
            ),
            "reference_umap_construction": (
                str(cached["construction"])
                if "construction" in cached.files
                else "legacy_cached_source_X_umap"
            ),
        }
    current_names = np.asarray(adata.obs_names).astype(str)
    if len(np.unique(current_names)) != len(current_names):
        raise ValueError("method output has duplicate cell identifiers")
    if not set(current_names).issubset(set(names)):
        raise ValueError("method cells are not a subset of the canonical UMAP cells")
    lookup = {name: index for index, name in enumerate(names)}
    adata.obsm["X_umap"] = coordinates[
        np.asarray([lookup[name] for name in current_names], dtype=np.int64)
    ]
    return provenance


def load_reference_moments(dataset: str, root: Path) -> dict:
    path = root / f"{dataset}.npz"
    if not path.is_file():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=False) as data:
        return {
            "path": path,
            "obs_names": np.asarray(data["obs_names"]).astype(str),
            "var_names": np.asarray(data["var_names"]).astype(str),
            "Ms": np.asarray(data["Ms"], dtype=np.float32),
            "Mu": np.asarray(data["Mu"], dtype=np.float32),
            "metadata": json.loads(str(data["metadata"])),
        }


def unique_casefold_lookup(names: np.ndarray, label: str) -> dict:
    lookup = {}
    duplicates = set()
    for position, name in enumerate(names):
        token = str(name).casefold()
        if token in lookup:
            duplicates.add(token)
        else:
            lookup[token] = position
    if duplicates:
        preview = ", ".join(sorted(duplicates)[:5])
        raise ValueError(f"{label}: duplicate case-insensitive gene names ({preview})")
    return lookup


def align_cells_to_reference(adata: ad.AnnData, reference_names: np.ndarray, row: pd.Series):
    method_names = np.asarray(adata.obs_names).astype(str)
    if len(np.unique(method_names)) != len(method_names):
        raise ValueError("method output has duplicate cell identifiers")
    method_lookup = {name: index for index, name in enumerate(method_names)}
    common_names = [name for name in reference_names if name in method_lookup]
    expected_common = int(row["n_common"])
    if len(common_names) != expected_common:
        raise ValueError(
            f"manifest common-cell mismatch: expected {expected_common}, found {len(common_names)}"
        )
    if len(common_names) < REFERENCE_PARAMS["n_neighbors"]:
        raise ValueError("fewer cells than the required 30 reference neighbors")
    method_positions = np.asarray([method_lookup[name] for name in common_names], dtype=np.int64)
    reference_lookup = {name: index for index, name in enumerate(reference_names)}
    reference_positions = np.asarray(
        [reference_lookup[name] for name in common_names], dtype=np.int64
    )
    aligned = adata[method_positions].copy()
    return aligned, reference_positions, {
        "reference_cell_count": int(len(reference_names)),
        "method_cell_count_before_alignment": int(len(method_names)),
        "scored_cell_count": int(len(common_names)),
        "reference_cells_excluded": int(len(reference_names) - len(common_names)),
        "nonreference_method_cells_excluded": int(len(method_names) - len(common_names)),
    }


GENE_ID_FIELDS = (
    "var_names",
    "query",
    "scopes",
    "gene_name_mapping",
    "symbol",
    "_id",
    "gene_ids",
    "gene_id",
    "Gene",
    "ensembl_id",
    "gene_symbol",
)


def align_genes_to_reference(
    adata: ad.AnnData,
    reference_names: np.ndarray,
    extra_fields: dict = None,
):
    reference_lookup = unique_casefold_lookup(reference_names, "reference moments")
    mappings = {}
    used_method_positions = set()
    field_counts = {}
    fields = list(GENE_ID_FIELDS)
    if extra_fields:
        fields.extend(extra_fields)
    for field in fields:
        if field == "var_names":
            values = np.asarray(adata.var_names).astype(str)
        elif extra_fields and field in extra_fields:
            values = np.asarray(extra_fields[field]).astype(str)
        elif field in adata.var:
            values = np.asarray(adata.var[field]).astype(str)
        else:
            continue
        if len(values) != adata.n_vars:
            continue
        matched = 0
        for method_position, value in enumerate(values):
            reference_position = reference_lookup.get(str(value).casefold())
            if (
                reference_position is not None
                and reference_position not in mappings
                and method_position not in used_method_positions
            ):
                mappings[reference_position] = method_position
                used_method_positions.add(method_position)
                matched += 1
        field_counts[field] = matched
    reference_positions = np.asarray(sorted(mappings), dtype=np.int64)
    method_positions = np.asarray(
        [mappings[position] for position in reference_positions], dtype=np.int64
    )
    if len(reference_positions) == 0:
        raise ValueError("no method genes can be mapped to the reference state")
    aligned = adata[:, method_positions].copy()
    return aligned, reference_positions, {
        "reference_gene_count": int(len(reference_names)),
        "scored_gene_count": int(len(reference_positions)),
        "gene_mapping": "one_to_one_casefold/reference_order",
        "mapping_field_counts": field_counts,
        "coordinate_space": "reference_gene_space",
    }


def graphvelo_dynamo_parent_fields(
    adata: ad.AnnData,
    dataset: str,
    inputs: pd.DataFrame,
    artifact_root: Path,
) -> dict:
    parent_row = select_row(inputs, "dynamo", dataset)
    parent_path = resolve_source(artifact_root, parent_row["result_path"])
    parent = ad.read_h5ad(parent_path, backed="r")
    try:
        # The Dynamo parent may contain symbols that differ only by case.  Do
        # not assign one arbitrarily: use an exact symbol when available, then
        # a case-folded symbol only when it is unique.
        parent_names = np.asarray(parent.var_names).astype(str)
        exact_positions = {}
        casefold_positions = {}
        for position, name in enumerate(parent_names):
            exact_positions.setdefault(name, []).append(position)
            casefold_positions.setdefault(name.casefold(), []).append(position)

        def parent_position(name: str) -> int:
            exact = exact_positions.get(name, [])
            if len(exact) == 1:
                return exact[0]
            folded = casefold_positions.get(name.casefold(), [])
            return folded[0] if len(folded) == 1 else -1

        parent_positions = np.asarray(
            [
                parent_position(str(name))
                for name in np.asarray(adata.var_names).astype(str)
            ],
            dtype=np.int64,
        )
        fields = {}
        for field in GENE_ID_FIELDS[1:]:
            if field not in parent.var:
                continue
            inherited = np.full(adata.n_vars, "", dtype=object)
            valid = parent_positions >= 0
            inherited[valid] = np.asarray(parent.var[field]).astype(str)[
                parent_positions[valid]
            ]
            fields["parent_dynamo_" + field] = inherited
        return fields
    finally:
        parent.file.close()


def inject_reference_state(
    adata: ad.AnnData,
    reference: dict,
    cell_positions: np.ndarray,
    gene_positions: np.ndarray,
) -> None:
    adata.layers["Ms"] = reference["Ms"][cell_positions][:, gene_positions]
    adata.layers["Mu"] = reference["Mu"][cell_positions][:, gene_positions]


def apply_reference_graph(
    adata: ad.AnnData,
    dataset: str,
    cell_positions: np.ndarray,
    reference_root: Path,
) -> dict:
    stored = load_reference_neighbors(dataset, output_root=reference_root)
    if len(cell_positions) == len(stored["obs_names"]):
        applied, detail = apply_reference_neighbors(
            adata, dataset, output_root=reference_root
        )
        if not applied:
            raise ValueError(detail)
        return {"mode": "stored_exact_30nn", "detail": detail}

    with np.load(reference_root / f"{dataset}.npz", allow_pickle=False) as data:
        reference_pca = np.asarray(data["reference_pca"], dtype=np.float32)
    adata.obsm["X_pca"] = reference_pca[cell_positions]
    import scvelo as scv
    scv.pp.neighbors(
        adata,
        n_neighbors=REFERENCE_PARAMS["n_neighbors"],
        n_pcs=REFERENCE_PARAMS["n_pcs"],
        use_rep=REFERENCE_PARAMS["use_rep"],
        metric=REFERENCE_PARAMS["metric"],
        random_state=REFERENCE_PARAMS["random_state"],
    )
    return {
        "mode": "recomputed_subset_30nn_from_reference_pca",
        "detail": f"recomputed 30NN for {adata.n_obs} retained reference cells",
    }


def clear_derived_graphs(adata: ad.AnnData) -> None:
    """Drop stale transition graphs without touching native 2-D velocities."""
    for key in list(adata.uns):
        if key in {
            "velocity_graph",
            "velocity_graph_neg",
            "Velocity_graph",
            "Velocity_graph_neg",
            "velo_s_norm_graph",
            "sde_velocity_graph",
            "new_velocity_graph",
            "new_velocity_graph_neg",
        }:
            del adata.uns[key]
    adata.obsp.pop("pearson_transition_matrix", None)


def install_native_2d_velocity(
    adata: ad.AnnData,
    row: pd.Series,
    dataset: str,
    native_coordinates: np.ndarray,
) -> dict:
    """Install a checked native 2-D velocity field in the canonical UMAP frame."""
    if dataset == "forebrain" or row["alignment_mode"] == "reference_graph_projection":
        return {
            "mode": "reference_graph_projection",
            "source_velocity_embedding_key": row["velocity_embedding_key"],
        }

    if row["alignment_mode"] != "direct":
        raise ValueError(
            f"{row['operational_id']}/{dataset}: non-direct 2-D geometry requires "
            "an explicit evaluation contract"
        )
    velocity_key = row["velocity_embedding_key"]
    if velocity_key not in adata.obsm:
        raise KeyError(
            f"{row['operational_id']}/{dataset}: missing obsm/{velocity_key!r}"
        )
    velocity = np.asarray(adata.obsm[velocity_key], dtype=np.float64)
    canonical = np.asarray(adata.obsm["X_umap"], dtype=np.float64)
    if native_coordinates.shape != canonical.shape or native_coordinates.shape[1] != 2:
        raise ValueError(
            f"{row['operational_id']}/{dataset}: native coordinate shape "
            f"{native_coordinates.shape} is not the canonical 2-D shape {canonical.shape}"
        )
    if velocity.shape != canonical.shape:
        raise ValueError(
            f"{row['operational_id']}/{dataset}: velocity shape {velocity.shape} "
            f"does not match canonical UMAP {canonical.shape}"
        )
    max_error = float(np.max(np.abs(native_coordinates - canonical)))
    # H5AD round trips store some persisted reference UMAPs as float32.
    if not np.allclose(native_coordinates, canonical, rtol=0.0, atol=1e-5):
        raise ValueError(
            f"{row['operational_id']}/{dataset}: native and canonical UMAP differ "
            f"(max_abs_error={max_error})"
        )
    adata.obsm["velocity_umap"] = velocity.astype(np.float32, copy=False)
    return {
        "mode": "native_stored_2d_velocity",
        "source_velocity_embedding_key": velocity_key,
        "native_coordinate_key": row["source_embedding_key"],
        "native_coordinate_max_abs_error": max_error,
    }


def tivelo_frozen_scores(
    adata: ad.AnnData,
    cluster_key: str,
    cluster_edges: list,
) -> dict:
    """Exact TiVelo 0.1.4 formulas, retained only for cross-environment QA."""
    labels = adata.obs[cluster_key]
    neighbor_indices = adata.uns["neighbors"]["indices"]
    highdim_raw = {}
    trans_raw = {}
    for source, target in cluster_edges:
        selected = np.asarray(labels == source)
        selected_neighbors = neighbor_indices[selected]
        boundary_nodes = [
            nodes[np.asarray(labels.iloc[nodes] == target)]
            for nodes in selected_neighbors
        ]
        source_velocities = adata.layers["velocity"][selected]
        source_expression = adata.layers["Ms"][selected]
        highdim_values = []
        trans_values = []
        graph_rows = adata.uns["velocity_graph"][selected]
        for x_pos, velocity, graph_row, nodes in zip(
            source_expression, source_velocities, graph_rows, boundary_nodes
        ):
            if len(nodes) == 0:
                continue
            displacement = adata.layers["Ms"][nodes] - x_pos
            highdim_values.append(
                float(cosine_similarity(displacement, velocity.reshape(1, -1)).mean())
            )
            trans_values.append(float(graph_row.toarray()[:, nodes].mean()))
        edge = f"{source}->{target}"
        highdim_raw[edge] = highdim_values
        trans_raw[edge] = trans_values

    def aggregate(raw: dict) -> float:
        values = [np.mean(item) for item in raw.values() if len(item)]
        values = [value for value in values if np.isfinite(value)]
        return float(np.mean(values)) if values else float("nan")

    return {
        "highdim_cbdir": aggregate(highdim_raw),
        "transcosine": aggregate(trans_raw),
        "highdim_cbdir_raw": highdim_raw,
        "transcosine_raw": trans_raw,
    }


def evaluate_four_native_metrics(adata: ad.AnnData, dataset: str, label: str, config: dict) -> dict:
    import unitvelo

    neighbors = adata.uns["neighbors"]["indices"]
    consistency = consistency_scores(
        dense(adata.layers["velocity"]),
        neighbors,
        REFERENCE_PARAMS["n_neighbors"],
    )

    evaluation_adata = adata
    if adata.uns["evaluation_metadata"]["method"].lower() == "dynamo" and "velocity_genes" in adata.var:
        evaluation_adata = adata[:, adata.var["velocity_genes"].astype(bool)].copy()
    output = io.StringIO()
    with redirect_stdout(output):
        unitvelo_result = unitvelo.evaluate(
            evaluation_adata,
            config["datasets"][dataset]["cluster_edges"],
            label,
            k_velocity="velocity",
            x_emb="X_umap",
        )
    cbdir, iccoh = parse_unitvelo_means(output.getvalue())
    cbdir_raw = unitvelo_result["Cross-Boundary Direction Correctness (A->B)"]
    cbdir_edge_raw = {
        f"{source}->{target}": [float(value) for value in values]
        for (source, target), values in cbdir_raw.items()
    }
    cbdir_edge_means = [np.mean(values) for values in cbdir_raw.values() if len(values)]
    # Match UnitVelo summary_scores: exact-zero edge means are omitted.
    cbdir_included = [value for value in cbdir_edge_means if value]
    cbdir_from_edges = float(np.mean(cbdir_included)) if cbdir_included else float("nan")
    if not np.isclose(cbdir, cbdir_from_edges, rtol=0.0, atol=1e-12, equal_nan=True):
        raise ValueError(
            "unitvelo CBDir aggregate differs from its native edge dict: "
            f"stdout={cbdir!r}, edge_mean={cbdir_from_edges!r}"
        )
    batc_score, batc_details = compute_batc(
        adata,
        dataset=dataset,
        cluster_edges=config["datasets"][dataset]["cluster_edges"],
        cluster_key=label,
        n_bins=config["datasets"][dataset]["batc_n_bins"],
        embedding_key="umap",
        velocity_key="velocity",
        store_in_adata=False,
        progress=False,
        return_details=True,
    )
    if not np.isclose(
        float(batc_score), float(batc_details["overall"]), rtol=0.0, atol=1e-12
    ):
        raise ValueError("BATC aggregate changed while collecting edge diagnostics")
    return {
        "consistency": float(np.mean(consistency)),
        "cbdir": float(cbdir),
        "iccoh": float(iccoh),
        "batc": float(np.mean(batc_score)),
        "cbdir_edge_raw": cbdir_edge_raw,
        "batc_edge_details": batc_details["edges"],
    }


def json_ready(value):
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--reference-neighbors-root",
        type=Path,
        default=DEFAULT_REFERENCE_NEIGHBORS_ROOT,
        help="fixed reference-30 assets; defaults to the repository copy",
    )
    parser.add_argument("--reference-moments-root", type=Path, required=True)
    parser.add_argument("--reference-umap-cache", type=Path, required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--prepared-output", type=Path, required=True)
    parser.add_argument("--indicators-output", type=Path, required=True)
    parser.add_argument("--n-jobs", type=int, default=8)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scvelo as scv

    config = load_config(args.config)
    inputs = read_inputs(args.inputs)
    row = select_row(inputs, args.method, args.dataset)
    source_path = resolve_source(args.artifact_root, row["result_path"])
    adata = ad.read_h5ad(source_path)
    invalid_velocity_values = prepare_velocity(adata, row["velocity_layer_key"])
    reference = load_reference_moments(args.dataset, args.reference_moments_root)
    adata, reference_cell_positions, cell_alignment = align_cells_to_reference(
        adata, reference["obs_names"], row
    )
    source_embedding_key = row["source_embedding_key"]
    if source_embedding_key not in adata.obsm:
        raise KeyError(
            f"{args.method}/{args.dataset}: missing obsm/{source_embedding_key!r}"
        )
    native_coordinates = np.asarray(
        adata.obsm[source_embedding_key], dtype=np.float64
    ).copy()
    if args.method == "veloae":
        if "Ms" not in adata.layers or adata.layers["Ms"].shape != adata.layers["velocity"].shape:
            raise ValueError("VeloAE native latent Ms does not match new_velocity")
        if "Mu" not in adata.layers or adata.layers["Mu"].shape != adata.layers["velocity"].shape:
            raise ValueError("VeloAE native latent Mu does not match new_velocity")
        gene_alignment = {
            "reference_gene_count": int(len(reference["var_names"])),
            "scored_gene_count": int(adata.n_vars),
            "gene_mapping": "not_applicable",
            "mapping_field_counts": {},
            "coordinate_space": "method_native_latent",
        }
        reference_state = {
            "path": None,
            "expression_mode": None,
            "smoothing_graph": "method_native_latent_state",
        }
    else:
        extra_gene_fields = (
            graphvelo_dynamo_parent_fields(
                adata, args.dataset, inputs, args.artifact_root
            )
            if args.method == "graphvelo_dynamo"
            else None
        )
        adata, reference_gene_positions, gene_alignment = align_genes_to_reference(
            adata, reference["var_names"], extra_fields=extra_gene_fields
        )
        inject_reference_state(
            adata, reference, reference_cell_positions, reference_gene_positions
        )
        reference_state = {
            "path": f"external:cell-reference-moments/{args.dataset}.npz",
            "expression_mode": reference["metadata"]["expression_mode"],
            "smoothing_graph": (
                "reference_neighbors/scvelo_0.2.5/"
                f"{REFERENCE_PARAMS['n_neighbors']}nn"
            ),
        }
    label = prepare_labels(adata, args.dataset, config)
    neighbor_application = apply_reference_graph(
        adata,
        args.dataset,
        reference_cell_positions,
        args.reference_neighbors_root,
    )
    umap_provenance = canonical_umap(
        adata,
        args.dataset,
        inputs,
        args.artifact_root,
        args.reference_umap_cache,
        args.reference_neighbors_root,
        config,
    )
    projection = install_native_2d_velocity(
        adata, row, args.dataset, native_coordinates
    )
    clear_derived_graphs(adata)
    scv.tl.velocity_graph(
        adata,
        vkey="velocity",
        xkey="Ms",
        n_neighbors=REFERENCE_PARAMS["n_neighbors"],
        n_recurse_neighbors=REFERENCE_PARAMS["n_recurse_neighbors"],
        mode_neighbors=REFERENCE_PARAMS["mode_neighbors"],
        n_jobs=args.n_jobs,
    )
    if projection["mode"] == "reference_graph_projection":
        scv.tl.velocity_embedding(
            adata, basis="umap", vkey="velocity", all_comps=True
        )
    metadata = {
        "method": args.method,
        "dataset": args.dataset,
        "result_path": row["result_path"],
        "result_sha256": sha256(source_path),
        "velocity_layer_key": row["velocity_layer_key"],
        "invalid_velocity_values_replaced_with_zero": invalid_velocity_values,
        "neighbor_graph": (
            "reference_neighbors/scvelo_0.2.5/"
            f"{REFERENCE_PARAMS['n_neighbors']}nn"
        ),
        "velocity_graph": (
            "scvelo_0.2.5/"
            f"n_neighbors={REFERENCE_PARAMS['n_neighbors']}"
            f"/n_recurse_neighbors={REFERENCE_PARAMS['n_recurse_neighbors']}"
            f"/{REFERENCE_PARAMS['mode_neighbors']}"
        ),
        "scvelo_version": scv.__version__,
        "label_key": label,
        "neighbor_application": neighbor_application,
        "cell_alignment": cell_alignment,
        "gene_alignment": gene_alignment,
        "reference_moments": reference_state,
        "velocity_embedding": projection,
        **umap_provenance,
    }
    adata.uns["evaluation_metadata"] = metadata
    frozen = tivelo_frozen_scores(
        adata, label, config["datasets"][args.dataset]["cluster_edges"]
    )
    native = evaluate_four_native_metrics(adata, args.dataset, label, config)
    plt.close("all")
    args.prepared_output.parent.mkdir(parents=True, exist_ok=True)
    args.indicators_output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.prepared_output, compression="gzip")
    payload = {
        "metadata": metadata,
        "metrics": {**native, **frozen},
    }
    args.indicators_output.write_text(
        json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "method": args.method,
                "dataset": args.dataset,
                "metrics": {
                    key: payload["metrics"][key]
                    for key in (
                        "consistency", "cbdir", "iccoh", "batc",
                        "highdim_cbdir", "transcosine",
                    )
                },
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
