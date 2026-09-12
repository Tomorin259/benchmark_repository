#!/usr/bin/env python3
"""Core signed-transition implementation for the four-method ensemble."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scvelo as scv
from scipy import sparse

from .geometric_median import (
    row_geometric_median4,
    validate_row_geometric_median4,
)


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MODULE_DIR / "config.json"

SCVELO_VERSION = "0.3.4"
H5AD_SCHEMA_VERSION = 6
FIXED_REFERENCE_NEIGHBOR_SCHEMA_VERSION = 1
METHOD_NAMES = ("TIvelo", "veloVI", "DeepVelo", "UniTVelo")
VARIANT_MEAN4 = "signed_t_mean4"
VARIANT_ROW_GEOMETRIC_MEDIAN4 = "signed_t_row_geometric_median4"
VARIANT_NAMES = (VARIANT_MEAN4, VARIANT_ROW_GEOMETRIC_MEDIAN4)
EMBEDDING_BASIS = "umap"
TEMPORARY_VKEY = "_ensemble_signed_input"

GRAPH_KWARGS = {
    "vkey": "velocity",
    "xkey": "Ms",
}
TRANSITION_KWARGS = {
    "scale": 10,
    "self_transitions": True,
    "use_negative_cosines": True,
    "weight_diffusion": 0,
    "n_neighbors": None,
}
FIXED_REFERENCE_NEIGHBOR_BUILD = {
    "n_neighbors": 30,
    "n_pcs": 30,
    "use_rep": "X_pca",
    "metric": "euclidean",
    "random_state": 0,
}

VARIANT_SPECS = {
    VARIANT_MEAN4: {
        "label": "Signed-T Mean-4",
        "definition": "Arithmetic mean of the four signed off-diagonal transition entries.",
        "method_weight": "fixed equal method weights of 1/4; no learned weights",
        "directory": "mean4",
    },
    VARIANT_ROW_GEOMETRIC_MEDIAN4: {
        "label": "Signed-T RowGeometricMedian-4",
        "plot_label": "GeoMedian-4",
        "definition": (
            "For each source cell, take the Euclidean geometric median of the "
            "four complete signed off-diagonal transition rows; preserve its "
            "native L1 mass without post-hoc renormalization."
        ),
        "method_weight": (
            "no learned method-identity weights; permutation-invariant convex "
            "solver coefficients are determined once per source cell by the "
            "four-row geometric-median geometry"
        ),
        "directory": "row_geometric_median4",
    },
}


def validate_scvelo_runtime() -> None:
    if str(scv.__version__) != SCVELO_VERSION:
        raise RuntimeError(
            f"ensemble requires scVelo {SCVELO_VERSION}, found {scv.__version__}"
        )


def protocol_metadata() -> dict:
    return {
        "scvelo_version": SCVELO_VERSION,
        "velocity_graph": copy.deepcopy(GRAPH_KWARGS),
        "transition_matrix": {
            key: "None" if value is None else copy.deepcopy(value)
            for key, value in TRANSITION_KWARGS.items()
        },
        "velocity_embedding": {
            "basis": EMBEDDING_BASIS,
            "autoscale": True,
            "explicit_transition_operator": True,
        },
        "fixed_neighbors": {
            **FIXED_REFERENCE_NEIGHBOR_BUILD,
            "asset_role": "dataset-level fixed reference graph",
            "reference_pca_installed": False,
        },
    }


def dense_float64(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values, dtype=np.float64)


def assert_unique_names(adata: ad.AnnData, label: str) -> None:
    if not adata.obs_names.is_unique:
        raise ValueError(f"{label} has duplicate cell identifiers")
    if not adata.var_names.is_unique:
        raise ValueError(f"{label} has duplicate gene identifiers")


def load_config(path: Path = DEFAULT_CONFIG) -> dict:
    with Path(path).open(encoding="utf-8") as handle:
        config = json.load(handle)
    if config.get("schema_version") != 1:
        raise ValueError("unsupported ensemble config schema")
    if tuple(config.get("methods", ())) != METHOD_NAMES:
        raise ValueError("config method order must be the four frozen methods")
    if set(config.get("datasets", {})) != {
        "bonemarrow",
        "capofibro",
        "dentategyrus",
        "dentate_lamanno",
        "forebrain",
        "gastrulation",
        "hematopoiesis",
        "mousebrain",
        "neuron",
        "organoid",
        "pancreas",
        "retina",
        "rpe1",
        "rpe1_gen",
    }:
        raise ValueError("config must contain exactly the frozen 14 datasets")
    return config


def _json_scalar(value):
    if isinstance(value, np.generic):
        value = value.item()
    elif isinstance(value, np.ndarray) and value.ndim == 0:
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return value


def _canonical_dataset(value) -> str:
    normalized = str(value).strip().lower().replace("-", "_")
    return {"dentategyrus_lamanno": "dentate_lamanno"}.get(normalized, normalized)


def _validate_sparse_asset(data, name: str, n_obs: int) -> sparse.csr_matrix:
    try:
        shape = tuple(int(value) for value in np.asarray(data[f"{name}_shape"]).ravel())
        values = np.asarray(data[f"{name}_data"])
        indices = np.asarray(data[f"{name}_indices"])
        indptr = np.asarray(data[f"{name}_indptr"])
    except KeyError as error:
        raise ValueError(f"fixed-neighbor asset is missing CSR field {error.args[0]!r}") from error
    if shape != (n_obs, n_obs):
        raise ValueError(f"fixed {name} matrix has an invalid shape")
    if values.ndim != 1 or indices.ndim != 1 or indptr.ndim != 1:
        raise ValueError(f"fixed {name} CSR arrays must be one-dimensional")
    if not np.issubdtype(indices.dtype, np.integer) or not np.issubdtype(indptr.dtype, np.integer):
        raise ValueError(f"fixed {name} CSR indices must be integers")
    if indptr.size != n_obs + 1 or indptr[0] != 0 or np.any(indptr[1:] < indptr[:-1]):
        raise ValueError(f"fixed {name} CSR indptr is invalid")
    if indptr[-1] != values.size or values.size != indices.size:
        raise ValueError(f"fixed {name} CSR arrays have inconsistent lengths")
    if np.any(indices < 0) or np.any(indices >= n_obs) or not np.isfinite(values).all():
        raise ValueError(f"fixed {name} CSR data are invalid")
    matrix = sparse.csr_matrix((values, indices, indptr), shape=shape)
    matrix.sort_indices()
    return matrix


def _validate_neighbor_metadata(metadata: dict) -> str:
    if not isinstance(metadata, dict) or metadata.get("dataset") is None:
        raise ValueError("fixed-neighbor metadata must contain a dataset")
    parameters = metadata.get("parameters")
    if not isinstance(parameters, dict):
        raise ValueError("fixed-neighbor metadata must contain parameters")
    for key, expected in FIXED_REFERENCE_NEIGHBOR_BUILD.items():
        actual = _json_scalar(parameters.get(key))
        if type(actual) is not type(expected) or actual != expected:
            raise ValueError(f"fixed-neighbor metadata {key!r} must be {expected!r}")
    return _canonical_dataset(metadata["dataset"])


def load_fixed_neighbor_asset(path: Path) -> dict:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"missing fixed-neighbor asset: {path}")
    with np.load(path, allow_pickle=False) as data:
        try:
            obs_names = np.asarray(data["obs_names"]).astype(str)
            indices = np.asarray(data["indices"])
            metadata = json.loads(str(_json_scalar(data["metadata"])))
        except KeyError as error:
            raise ValueError(f"{path}: missing fixed-neighbor field {error.args[0]!r}") from error
        if len(np.unique(obs_names)) != len(obs_names):
            raise ValueError("fixed-neighbor asset has duplicate cell identifiers")
        n_obs = len(obs_names)
        expected_shape = (n_obs, FIXED_REFERENCE_NEIGHBOR_BUILD["n_neighbors"])
        if indices.shape != expected_shape or not np.issubdtype(indices.dtype, np.integer):
            raise ValueError(f"fixed-neighbor indices must have shape {expected_shape}")
        if np.any(indices < 0) or np.any(indices >= n_obs):
            raise ValueError("fixed-neighbor indices contain out-of-range cells")
        if not np.array_equal(indices[:, 0], np.arange(n_obs)):
            raise ValueError("fixed-neighbor indices must store the source cell first")
        dataset = _validate_neighbor_metadata(metadata)
        asset = {
            "path": path,
            "dataset": dataset,
            "obs_names": obs_names,
            "indices": indices.astype(np.int64, copy=False),
            "distances": _validate_sparse_asset(data, "distances", n_obs),
            "connectivities": _validate_sparse_asset(data, "connectivities", n_obs),
            "metadata": metadata,
            "sha256": sha256(path),
        }
    return asset


def install_fixed_neighbors(adata: ad.AnnData, asset: dict) -> dict:
    """Install one fixed graph and map it into the AnnData cell order."""
    assert_unique_names(adata, "reference")
    asset_names = np.asarray(asset["obs_names"]).astype(str)
    current_names = np.asarray(adata.obs_names).astype(str)
    if len(np.unique(asset_names)) != len(asset_names):
        raise ValueError("fixed-neighbor asset has duplicate cell identifiers")
    if set(asset_names) != set(current_names):
        raise ValueError("fixed-neighbor asset cell set differs from AnnData")
    indices = np.asarray(asset["indices"], dtype=np.int64)
    n_obs = len(asset_names)
    if indices.shape != (n_obs, FIXED_REFERENCE_NEIGHBOR_BUILD["n_neighbors"]):
        raise ValueError("fixed-neighbor indices have an invalid shape")
    if np.any(indices < 0) or np.any(indices >= n_obs):
        raise ValueError("fixed-neighbor indices contain out-of-range cells")
    if not np.array_equal(indices[:, 0], np.arange(n_obs)):
        raise ValueError("fixed-neighbor indices must store the source cell first")
    for name in ("distances", "connectivities"):
        matrix = sparse.csr_matrix(asset[name])
        if matrix.shape != (n_obs, n_obs):
            raise ValueError(f"fixed {name} matrix has an invalid shape")

    asset_position = {name: index for index, name in enumerate(asset_names)}
    asset_rows = np.asarray(
        [asset_position[name] for name in current_names], dtype=np.int64
    )
    current_position = np.empty(n_obs, dtype=np.int64)
    current_position[asset_rows] = np.arange(n_obs, dtype=np.int64)
    mapped_indices = current_position[indices[asset_rows]]

    adata.obsp["distances"] = (
        sparse.csr_matrix(asset["distances"])[asset_rows][:, asset_rows].tocsr()
    )
    adata.obsp["connectivities"] = (
        sparse.csr_matrix(asset["connectivities"])[asset_rows][:, asset_rows].tocsr()
    )
    adata.uns["neighbors"] = {
        "indices": mapped_indices,
        "distances_key": "distances",
        "connectivities_key": "connectivities",
        "params": {
            "method": "umap",
            **FIXED_REFERENCE_NEIGHBOR_BUILD,
            "reference_dataset": asset["dataset"],
        },
        "fixed_reference_neighbors": {
            "schema_version": FIXED_REFERENCE_NEIGHBOR_SCHEMA_VERSION,
            "artifact_id": f"neighbor:{asset['dataset']}",
            "dataset": asset["dataset"],
            "sha256": asset["sha256"],
            "build": dict(FIXED_REFERENCE_NEIGHBOR_BUILD),
            "reference_pca": {
                "installed": False,
                "semantics": (
                    "asset reference_pca is not installed; only graph and indices are installed"
                ),
            },
        },
    }
    return {
        "n_obs": n_obs,
        "indices_shape": list(mapped_indices.shape),
        "distances_nnz": int(adata.obsp["distances"].nnz),
        "connectivities_nnz": int(adata.obsp["connectivities"].nnz),
        "parameters": protocol_metadata()["fixed_neighbors"],
        "stable_metadata": copy.deepcopy(
            adata.uns["neighbors"]["fixed_reference_neighbors"]
        ),
    }


def _velocity_gene_mask(adata: ad.AnnData, method_name: str) -> np.ndarray:
    if "velocity_genes" not in adata.var:
        velocity = dense_float64(adata.layers["velocity"])
        return np.isfinite(velocity).all(axis=0) & np.any(velocity != 0, axis=0)
    raw = np.asarray(adata.var["velocity_genes"])
    if raw.ndim != 1:
        raise ValueError(f"{method_name}: velocity_genes must be one-dimensional")
    if raw.dtype.kind == "b":
        return raw.astype(bool, copy=True)
    if raw.dtype.kind in "iuf":
        values = np.asarray(raw, dtype=np.float64)
        if not np.isfinite(values).all() or not np.isin(values, [0.0, 1.0]).all():
            raise ValueError(
                f"{method_name}: velocity_genes must contain only 0/1 values"
            )
        return values.astype(bool)
    raise TypeError(f"{method_name}: velocity_genes must be boolean or numeric 0/1")


def align_method(
    reference: ad.AnnData,
    method_adata: ad.AnnData,
    method_name: str,
) -> ad.AnnData:
    """Align cells, neighbors, PCA and UMAP while preserving native gene space."""
    assert_unique_names(reference, "reference")
    assert_unique_names(method_adata, method_name)
    missing = reference.obs_names.difference(method_adata.obs_names)
    extra = method_adata.obs_names.difference(reference.obs_names)
    if len(missing) or len(extra):
        raise ValueError(
            f"{method_name} cell set differs from reference: "
            f"missing={len(missing)}, extra={len(extra)}"
        )
    for layer in ("Ms", "velocity"):
        if layer not in method_adata.layers:
            raise KeyError(f"{method_name} is missing layers['{layer}']")

    aligned = method_adata[reference.obs_names, :].copy()
    for key in list(aligned.obsp.keys()):
        del aligned.obsp[key]
    for key in ("distances", "connectivities"):
        if key not in reference.obsp:
            raise KeyError(f"reference is missing obsp['{key}']")
        aligned.obsp[key] = reference.obsp[key].copy()
    aligned.uns["neighbors"] = copy.deepcopy(reference.uns["neighbors"])
    for key in ("X_pca", "X_umap"):
        if key not in reference.obsm:
            raise KeyError(f"reference is missing obsm['{key}']")
        aligned.obsm[key] = np.asarray(reference.obsm[key]).copy()

    mask = _velocity_gene_mask(aligned, method_name)
    if mask.shape != (aligned.n_vars,) or not np.any(mask):
        raise ValueError(f"{method_name} has no valid velocity genes")
    aligned.var["velocity_genes"] = mask
    return aligned


def compute_velocity_graph(
    adata: ad.AnnData,
    *,
    n_jobs: int,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Run scVelo's graph routine using the frozen tutorial defaults."""
    validate_scvelo_runtime()
    n_jobs = int(n_jobs)
    if n_jobs < 1:
        raise ValueError("n_jobs must be positive")
    adata.uns.pop("velocity_graph", None)
    adata.uns.pop("velocity_graph_neg", None)
    scv.tl.velocity_graph(
        adata,
        **GRAPH_KWARGS,
        n_jobs=n_jobs,
        show_progress_bar=False,
    )
    positive = sparse.csr_matrix(adata.uns["velocity_graph"]).copy()
    negative = sparse.csr_matrix(adata.uns["velocity_graph_neg"]).copy()
    positive.eliminate_zeros()
    negative.eliminate_zeros()
    return positive, negative


def compute_signed_transition(
    adata: ad.AnnData,
    positive: sparse.spmatrix,
    negative: sparse.spmatrix,
) -> sparse.csr_matrix:
    """Convert signed cosine graphs through scVelo's signed transition path."""
    validate_scvelo_runtime()
    graph_key = f"{TEMPORARY_VKEY}_graph"
    negative_key = f"{TEMPORARY_VKEY}_graph_neg"
    if graph_key in adata.uns or negative_key in adata.uns:
        raise KeyError("temporary ensemble graph fields already exist")
    adata.uns[graph_key] = sparse.csr_matrix(positive).astype(np.float64)
    adata.uns[negative_key] = sparse.csr_matrix(negative).astype(np.float64)
    try:
        transition = scv.tl.transition_matrix(
            adata,
            vkey=TEMPORARY_VKEY,
            **TRANSITION_KWARGS,
        )
    finally:
        adata.uns.pop(graph_key, None)
        adata.uns.pop(negative_key, None)
    transition = sparse.csr_matrix(transition).astype(np.float64)
    transition.eliminate_zeros()
    if transition.shape[0] != transition.shape[1]:
        raise ValueError("signed transition matrix must be square")
    return transition


def offdiagonal(matrix) -> sparse.csr_matrix:
    value = sparse.csr_matrix(matrix).astype(np.float64).copy()
    value.setdiag(0.0)
    value.eliminate_zeros()
    return value


def row_l1_norm(matrix) -> np.ndarray:
    return np.asarray(np.abs(sparse.csr_matrix(matrix)).sum(axis=1)).ravel()


def align_method_values(
    matrices: list[sparse.spmatrix],
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Align four sparse operators on the deterministic union edge support."""
    if len(matrices) != len(METHOD_NAMES):
        raise ValueError("exactly four transition matrices are required")
    normalized = []
    shape = None
    for matrix in matrices:
        current = sparse.csr_matrix(matrix).astype(np.float64).copy()
        current.eliminate_zeros()
        current.sort_indices()
        if current.shape[0] != current.shape[1]:
            raise ValueError("transition matrices must be square")
        if shape is None:
            shape = current.shape
        if current.shape != shape:
            raise ValueError("transition matrices must share one shape")
        normalized.append(current)

    n_obs = shape[0]
    keyed = []
    for matrix in normalized:
        coo = matrix.tocoo(copy=False)
        keys = coo.row.astype(np.int64) * np.int64(n_obs) + coo.col.astype(np.int64)
        order = np.argsort(keys, kind="stable")
        keyed.append((keys[order], np.asarray(coo.data, dtype=np.float64)[order]))
    union = np.unique(np.concatenate([keys for keys, _ in keyed]))
    values = np.zeros((len(normalized), union.size), dtype=np.float64)
    for method_id, (keys, data) in enumerate(keyed):
        values[method_id, np.searchsorted(union, keys)] = data
    rows = union // np.int64(n_obs)
    cols = union % np.int64(n_obs)
    template = sparse.csr_matrix(
        (np.ones(union.size, dtype=np.float64), (rows, cols)),
        shape=shape,
    )
    template.sort_indices()
    return template, values


def aggregate_offdiagonal(
    template: sparse.spmatrix,
    values: np.ndarray,
    variant: str,
) -> tuple[np.ndarray, dict]:
    """Aggregate four complete signed off-diagonal rows."""
    template = sparse.csr_matrix(template)
    values = np.asarray(values, dtype=np.float64)
    if values.shape != (len(METHOD_NAMES), template.nnz):
        raise ValueError("aggregation input must have shape (4, template.nnz)")
    if not np.isfinite(values).all():
        raise ValueError("aggregation input contains NaN or infinity")
    if variant == VARIANT_MEAN4:
        aggregate = values.mean(axis=0)
        diagnostics = {"rule": "arithmetic mean", "n_methods": 4}
    elif variant == VARIANT_ROW_GEOMETRIC_MEDIAN4:
        aggregate, diagnostics, _ = row_geometric_median4(
            template, values, return_cell_diagnostics=True
        )
        diagnostics["invariants"] = validate_row_geometric_median4()
        diagnostics["n_methods"] = 4
    else:
        raise ValueError(f"unknown ensemble variant: {variant}")
    return aggregate, diagnostics


def complete_signed_operator(
    template: sparse.spmatrix,
    aggregate: np.ndarray,
    *,
    tolerance: float = 1e-10,
) -> tuple[sparse.csr_matrix, dict]:
    """Place the remaining signed-L1 mass on the diagonal."""
    template = sparse.csr_matrix(template)
    aggregate = np.asarray(aggregate, dtype=np.float64)
    if aggregate.ndim != 1 or aggregate.size != template.nnz:
        raise ValueError("aggregate does not match the common edge support")
    if not np.isfinite(aggregate).all():
        raise ValueError("aggregate contains NaN or infinity")

    off = sparse.csr_matrix(
        (aggregate.copy(), template.indices.copy(), template.indptr.copy()),
        shape=template.shape,
    )
    off.setdiag(0.0)
    off.eliminate_zeros()
    mass = row_l1_norm(off)
    if np.any(mass > 1.0 + tolerance):
        maximum = float(np.max(mass))
        raise ValueError(f"aggregated signed-L1 row mass exceeds one: {maximum}")
    diagonal = np.maximum(1.0 - mass, 0.0)
    result = (off + sparse.diags(diagonal, format="csr")).tocsr()
    result.eliminate_zeros()
    validate_signed_operator(result, "completed ensemble operator")
    return result, {
        "directional_mass_mean": float(mass.mean()) if mass.size else 0.0,
        "neutral_mass_mean": float(diagonal.mean()) if diagonal.size else 0.0,
        "n_rows_rescaled": 0,
        "n_edges": int(template.nnz),
    }


def validate_signed_operator(
    matrix,
    label: str,
    *,
    tolerance: float = 1e-8,
) -> None:
    matrix = sparse.csr_matrix(matrix)
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{label} must be square")
    if matrix.nnz and not np.isfinite(matrix.data).all():
        raise ValueError(f"{label} contains NaN or infinity")
    masses = row_l1_norm(matrix)
    if not np.allclose(masses, 1.0, atol=tolerance, rtol=0.0):
        raise ValueError(
            f"{label} is not signed-L1 normalized: "
            f"min={masses.min()}, max={masses.max()}"
        )


def _expression_displacement(operator, expression) -> np.ndarray:
    expression = dense_float64(expression)
    if not np.isfinite(expression).all():
        raise ValueError("expression layer contains NaN or infinity")
    displacement = sparse.csr_matrix(operator).dot(expression) - expression
    return dense_float64(displacement)


def project_with_scvelo(
    reference: ad.AnnData,
    operator: sparse.spmatrix,
    *,
    variant: str,
    logical_inputs: dict[str, dict],
) -> tuple[ad.AnnData, int]:
    """Materialize the aggregate and obtain its explicit-T UMAP projection."""
    validate_scvelo_runtime()
    if variant not in VARIANT_SPECS:
        raise ValueError(f"unknown ensemble variant: {variant}")
    operator = sparse.csr_matrix(operator).astype(np.float64)
    validate_signed_operator(operator, "projection operator")
    if "Ms" not in reference.layers:
        raise KeyError("reference is missing layers['Ms']")
    if "X_umap" not in reference.obsm:
        raise KeyError("reference is missing obsm['X_umap']")

    result = reference.copy()
    result.obsp["ensemble_transition_signed"] = operator.copy()
    result.layers["ensemble_velocity"] = _expression_displacement(
        operator, result.layers["Ms"]
    ).astype(np.float32)
    result.var["ensemble_velocity_genes"] = np.ones(result.n_vars, dtype=bool)
    scv.tl.velocity_embedding(
        result,
        basis=EMBEDDING_BASIS,
        vkey="ensemble_velocity",
        T=operator.copy(),
        autoscale=True,
        all_comps=True,
    )
    embedded_key = f"ensemble_velocity_{EMBEDDING_BASIS}"
    embedded = np.asarray(result.obsm[embedded_key], dtype=np.float64)
    invalid = ~np.isfinite(embedded).all(axis=1)
    if np.any(invalid):
        embedded = embedded.copy()
        embedded[invalid] = 0.0
        result.obsm[embedded_key] = embedded.astype(np.float32)

    result.uns["ensemble_signed_params"] = {
        "schema_version": H5AD_SCHEMA_VERSION,
        "workflow": "scvelo_034_signed_t_ensemble",
        "aggregation_variant": variant,
        "aggregation_label": VARIANT_SPECS[variant]["label"],
        "aggregation_definition": VARIANT_SPECS[variant]["definition"],
        "methods": list(METHOD_NAMES),
        "method_weight": VARIANT_SPECS[variant]["method_weight"],
        "scvelo_protocol": protocol_metadata(),
        "diagonal_completion": "T_ii = 1 - sum_j!=i(abs(T_ij))",
        "operator_semantics": (
            "signed L1-normalized scVelo projection operator; not a "
            "non-negative Markov probability matrix"
        ),
        "compatibility_velocity": {
            "field": "layers['ensemble_velocity']",
            "definition": "ensemble_transition_signed @ Ms - Ms",
            "expression_layer": "Ms",
            "expression_source": "canonical reference",
            "compatibility_only": True,
            "biological_interpretation": (
                "not an aggregation of the methods' gene-level velocity layers"
            ),
        },
        "fixed_reference_neighbors": copy.deepcopy(
            result.uns.get("neighbors", {}).get("fixed_reference_neighbors", {})
        ),
        "fixed_neighbor_semantics": (
            "only distances, connectivities, and neighbor indices are installed; "
            "the asset reference_pca is not installed into obsm['X_pca']"
        ),
        "logical_inputs": copy.deepcopy(logical_inputs),
    }
    return result, int(invalid.sum())


def write_h5ad(adata: ad.AnnData, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    adata.write_h5ad(temporary, compression="gzip")
    temporary.replace(path)


def save_stream_plot(
    adata: ad.AnnData,
    path: Path,
    *,
    label_key: str,
    title: str,
) -> None:
    import matplotlib.pyplot as plt

    path = Path(path)
    if label_key not in adata.obs:
        raise KeyError(f"stream color key is missing: {label_key}")
    axis = scv.pl.velocity_embedding_stream(
        adata,
        basis=EMBEDDING_BASIS,
        vkey="ensemble_velocity",
        color=label_key,
        title=title,
        show=False,
    )
    figure = axis.figure if hasattr(axis, "figure") else plt.gcf()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    figure.savefig(temporary, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError(f"empty stream plot: {path}")
    temporary.replace(path)
    plt.close(figure)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
