#!/usr/bin/env python3
"""Validate the public fixed-reference Cell evaluation bundle."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy import sparse


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MODULE_DIR / "config.json"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
WINDOWS_ABSOLUTE_RE = re.compile(r"^[A-Za-z]:[\\/]")
PRIVATE_PATH_MARKERS = ("/" + "mnt" + "/", "/" + "Users" + "/", chr(92) + "Users" + chr(92))

EXPECTED_DATASET_COUNT = 14
EXPECTED_METHOD_COUNT = 20
EXPECTED_DATASET_KEYS = {
    "bonemarrow", "capofibro", "DentateGyrus_lamanno", "dentategyrus",
    "forebrain", "gastrulation", "hematopoiesis", "mousebrain", "neuron",
    "organoid", "pancreas", "retina", "rpe1", "rpe1_gen",
}
EXPECTED_METHOD_KEYS = {
    "scvelo_deterministic", "scvelo_stochastic", "scvelo_dynamic", "multivelo",
    "unitvelo", "unitvelo_mode2", "dynamo", "TIvelo", "cell2fate", "TFvelo",
    "veloae", "veloVI", "DeepKINET", "TSvelo", "celldancer", "deepvelo",
    "SDEvelo", "deepvelo2022", "graphvelo_scv", "graphvelo_dynamo",
}
EXPECTED_PRIMARY_COLUMNS = [
    "dataset", "method", "consistency", "iccoh", "cbdir", "transcosine", "batc",
]
EXPECTED_SUPPLEMENTARY_COLUMNS = [
    "dataset", "method", "highdim_cbdir",
    "transcosine_evaluable_edges", "transcosine_evaluable_source_cells",
    "highdim_cbdir_evaluable_edges", "highdim_cbdir_evaluable_source_cells",
]
EXPECTED_AUDIT_COLUMNS = [
    "dataset", "method", "result_path", "result_sha256",
    "velocity_layer_key", "invalid_velocity_values_replaced_with_zero", "status",
]
REFERENCE_ASSET_FIELDS = {
    "obs_names",
    "indices",
    "reference_pca",
    "reference_var_names",
    "metadata",
    "distances_data",
    "distances_indices",
    "distances_indptr",
    "distances_shape",
    "connectivities_data",
    "connectivities_indices",
    "connectivities_indptr",
    "connectivities_shape",
}


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"missing canonical table: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        return list(reader.fieldnames), list(reader)


def _check_table(
    path: Path,
    expected_columns: list[str],
    expected_keys: set[tuple[str, str]],
) -> list[dict[str, str]]:
    columns, rows = _read_csv(path)
    if columns != expected_columns:
        raise ValueError(
            f"{path}: columns differ from config: expected {expected_columns}, got {columns}"
        )
    keys = {(row.get("dataset", ""), row.get("method", "")) for row in rows}
    if keys != expected_keys or len(rows) != len(expected_keys):
        raise ValueError(f"{path}: expected one row for each dataset/method pair")
    return rows


def _finite(row: dict[str, str], column: str, path: Path) -> None:
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: {column} is not numeric") from error
    if not math.isfinite(value):
        raise ValueError(f"{path}: {column} is not finite")


def _nonnegative_int(row: dict[str, str], column: str, path: Path) -> None:
    try:
        value = int(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: {column} is not an integer") from error
    if value < 0:
        raise ValueError(f"{path}: {column} is negative")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _check_logical_locator(value: str, label: str) -> None:
    if (
        not isinstance(value, str)
        or not value
        or Path(value).is_absolute()
        or WINDOWS_ABSOLUTE_RE.match(value)
        or any(marker in value for marker in PRIVATE_PATH_MARKERS)
    ):
        raise ValueError(f"{label}: expected a non-empty logical locator")
    if ".." in Path(value).parts:
        raise ValueError(f"{label}: parent traversal is not allowed")


def _positive_int(value, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label}: expected a positive integer")
    return value


def _text_vector(values, expected_length: int, label: str) -> np.ndarray:
    values = np.asarray(values)
    if values.ndim != 1 or values.shape[0] != expected_length:
        raise ValueError(
            f"{label}: expected a one-dimensional vector of length {expected_length}, "
            f"got {values.shape}"
        )
    if values.dtype.kind == "S":
        values = values.astype("U")
    values = values.astype(str)
    if any(not value for value in values):
        raise ValueError(f"{label}: contains an empty identifier")
    if len(set(values.tolist())) != expected_length:
        raise ValueError(f"{label}: contains duplicate identifiers")
    return values


def _validate_csr(data, name: str, n_obs: int, path: Path) -> None:
    fields = {
        f"{name}_data", f"{name}_indices", f"{name}_indptr", f"{name}_shape",
    }
    missing = sorted(fields.difference(data.files))
    if missing:
        raise ValueError(f"{path}: {name} CSR payload is missing {missing}")
    shape = np.asarray(data[f"{name}_shape"])
    if (
        shape.ndim != 1
        or shape.shape[0] != 2
        or not np.issubdtype(shape.dtype, np.integer)
        or tuple(int(value) for value in shape) != (n_obs, n_obs)
    ):
        raise ValueError(f"{path}: {name} has invalid matrix shape {shape.tolist()}")
    values = np.asarray(data[f"{name}_data"])
    indices = np.asarray(data[f"{name}_indices"])
    indptr = np.asarray(data[f"{name}_indptr"])
    if values.ndim != 1 or indices.ndim != 1 or indptr.ndim != 1:
        raise ValueError(f"{path}: {name} CSR components must be one-dimensional")
    if not np.issubdtype(values.dtype, np.number):
        raise ValueError(f"{path}: {name} data must be numeric")
    if not np.issubdtype(indices.dtype, np.integer):
        raise ValueError(f"{path}: {name} indices must be integer")
    if not np.issubdtype(indptr.dtype, np.integer):
        raise ValueError(f"{path}: {name} indptr must be integer")
    if not np.isfinite(values).all():
        raise ValueError(f"{path}: {name} data contains non-finite values")
    if indices.shape[0] != values.shape[0] or indptr.shape[0] != n_obs + 1:
        raise ValueError(f"{path}: {name} CSR component lengths are inconsistent")
    if (
        indptr[0] != 0
        or indptr[-1] != values.shape[0]
        or np.any(np.diff(indptr) < 0)
    ):
        raise ValueError(f"{path}: {name} indptr is malformed")
    if indices.size and (indices.min() < 0 or indices.max() >= n_obs):
        raise ValueError(f"{path}: {name} indices are out of range")
    matrix = sparse.csr_matrix((values, indices, indptr), shape=(n_obs, n_obs))
    matrix.check_format(full_check=True)


def _validate_reference_asset(
    path: Path,
    dataset: str,
    entry: dict,
    manifest_contract: dict,
) -> None:
    n_obs = entry["n_obs"]
    try:
        with np.load(path, allow_pickle=False) as data:
            fields = set(data.files)
            if fields != REFERENCE_ASSET_FIELDS:
                missing = sorted(REFERENCE_ASSET_FIELDS.difference(fields))
                extra = sorted(fields.difference(REFERENCE_ASSET_FIELDS))
                raise ValueError(
                    f"{path}: NPZ fields differ from the frozen structure; "
                    f"missing={missing}, extra={extra}"
                )

            obs_names = _text_vector(data["obs_names"], n_obs, f"{dataset} obs_names")
            indices = np.asarray(data["indices"])
            if not np.issubdtype(indices.dtype, np.integer):
                raise ValueError(f"{path}: neighbor indices must be integer")
            if indices.shape != (n_obs, manifest_contract["n_neighbors"]):
                raise ValueError(f"{path}: invalid neighbor-index shape {indices.shape}")
            if indices.size and (indices.min() < 0 or indices.max() >= n_obs):
                raise ValueError(f"{path}: neighbor indices are out of range")
            if not np.array_equal(indices[:, 0], np.arange(n_obs)):
                raise ValueError(f"{path}: neighbor indices must store the source cell first")

            reference_pca = np.asarray(data["reference_pca"])
            if (
                reference_pca.ndim != 2
                or reference_pca.shape[0] != n_obs
                or reference_pca.shape[1] < manifest_contract["n_pcs"]
                or not np.issubdtype(reference_pca.dtype, np.number)
                or not np.isfinite(reference_pca).all()
            ):
                raise ValueError(
                    f"{path}: reference_pca has invalid shape or values "
                    f"{reference_pca.shape}"
                )
            reference_var_names = _text_vector(
                data["reference_var_names"],
                entry["n_vars_after_filter"],
                f"{dataset} reference_var_names",
            )

            metadata_raw = np.asarray(data["metadata"])
            if metadata_raw.size != 1:
                raise ValueError(f"{path}: metadata must be a scalar JSON string")
            metadata_value = metadata_raw.reshape(-1)[0]
            if isinstance(metadata_value, bytes):
                metadata_value = metadata_value.decode("utf-8")
            metadata = json.loads(str(metadata_value))
            if not isinstance(metadata, dict):
                raise ValueError(f"{path}: metadata must decode to an object")
            if metadata.get("dataset") != dataset or metadata.get("n_obs") != n_obs:
                raise ValueError(f"{path}: metadata dataset/n_obs do not match manifest")
            if metadata.get("n_top_genes_used") != entry["n_top_genes_used"]:
                raise ValueError(f"{path}: metadata n_top_genes_used does not match manifest")
            if metadata.get("n_vars_after_filter") != entry["n_vars_after_filter"]:
                raise ValueError(f"{path}: metadata n_vars_after_filter does not match manifest")
            parameters = metadata.get("parameters")
            if not isinstance(parameters, dict) or any(
                parameters.get(key) != value
                for key, value in manifest_contract.items()
            ):
                raise ValueError(f"{path}: metadata parameters do not match manifest contract")
            source = metadata.get("source")
            _check_logical_locator(source, f"{path}: metadata source")
            if not source.startswith("velocity:"):
                raise ValueError(f"{path}: metadata source must use velocity: locator")

            _validate_csr(data, "distances", n_obs, path)
            _validate_csr(data, "connectivities", n_obs, path)
    except Exception as error:
        if isinstance(error, ValueError) and str(error).startswith(f"{path}:"):
            raise
        raise ValueError(f"{dataset}: unable to load reference-neighbor NPZ {path}: {error}") from error


def _validate_reference_neighbors(root: Path, datasets: list[str], reference: dict) -> None:
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError(f"missing reference-neighbor manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise ValueError("reference-neighbor manifest has an unsupported schema version")
    n_neighbors = _positive_int(reference.get("n_neighbors"), "config reference n_neighbors")
    n_pcs = _positive_int(reference.get("n_pcs"), "config reference n_pcs")
    expected_contract = {
        "metric": reference.get("metric"),
        "min_shared_counts": 20,
        "n_neighbors": n_neighbors,
        "n_pcs": n_pcs,
        "n_top_genes_requested": 2000,
        "random_state": reference.get("random_state"),
        "use_rep": reference.get("use_rep"),
    }
    if manifest.get("contract") != expected_contract:
        raise ValueError("reference-neighbor manifest contract differs from config")
    manifest_datasets = manifest.get("datasets")
    if not isinstance(manifest_datasets, dict) or set(manifest_datasets) != set(datasets):
        raise ValueError("reference-neighbor manifest does not cover the configured datasets")
    for dataset in datasets:
        entry = manifest_datasets[dataset]
        if not isinstance(entry, dict):
            raise ValueError(f"{dataset}: reference-neighbor manifest entry is not an object")
        filename = entry.get("file", "")
        relative = Path(filename)
        if relative.is_absolute() or relative.name != filename or relative.suffix != ".npz":
            raise ValueError(f"{dataset}: invalid reference-neighbor filename")
        asset = root / relative
        if not asset.is_file():
            raise ValueError(f"{dataset}: missing reference-neighbor asset {asset}")
        n_obs = _positive_int(entry.get("n_obs"), f"{dataset}: n_obs")
        n_vars_after_filter = _positive_int(
            entry.get("n_vars_after_filter"), f"{dataset}: n_vars_after_filter"
        )
        n_top_genes_used = _positive_int(
            entry.get("n_top_genes_used"), f"{dataset}: n_top_genes_used"
        )
        expected_sha = entry.get("sha256", "")
        if not SHA256_RE.fullmatch(expected_sha):
            raise ValueError(f"{dataset}: invalid reference-neighbor SHA256")
        if _sha256(asset) != expected_sha:
            raise ValueError(f"{dataset}: reference-neighbor SHA256 mismatch")
        shape = entry.get("indices_shape")
        if shape != [n_obs, n_neighbors]:
            raise ValueError(f"{dataset}: reference-neighbor shape is not reference-{n_neighbors}")
        _validate_reference_asset(
            asset,
            dataset,
            {
                **entry,
                "n_obs": n_obs,
                "n_vars_after_filter": n_vars_after_filter,
                "n_top_genes_used": n_top_genes_used,
            },
            expected_contract,
        )


def _validate_audit_sources(audit_rows: list[dict[str, str]], expected_keys: set[tuple[str, str]]) -> None:
    """Match audit paths, hashes, and layer keys to the input manifest."""
    inputs_path = MODULE_DIR / "inputs.tsv"
    if not inputs_path.is_file():
        raise ValueError("Cell input manifest is required for audit validation")
    with inputs_path.open(encoding="utf-8", newline="") as handle:
        input_reader = csv.DictReader(handle, delimiter="\t")
        input_rows = list(input_reader)
    required = {"operational_id", "evaluation_dataset", "result_path", "sha256", "velocity_layer_key"}
    if not input_reader.fieldnames or not required.issubset(input_reader.fieldnames):
        raise ValueError(f"{inputs_path}: missing audit source columns")
    by_key = {}
    for row in input_rows:
        key = (row.get("evaluation_dataset", ""), row.get("operational_id", ""))
        if key in by_key:
            raise ValueError(f"{inputs_path}: duplicate audit source key {key}")
        by_key[key] = row
    if set(by_key) != expected_keys:
        raise ValueError(f"{inputs_path}: audit source key space differs from the result tables")
    for row in audit_rows:
        key = (row["dataset"], row["method"])
        expected = (by_key[key]["result_path"], by_key[key]["sha256"], by_key[key]["velocity_layer_key"])
        actual = (row["result_path"], row["result_sha256"], row["velocity_layer_key"])
        if actual != expected:
            raise ValueError(f"cell_metrics_audit.csv: source contract mismatch for {key}")


def _validate_streams(root: Path, methods: list[str]) -> None:
    manifest_path = root / "render_manifest.csv"
    columns, rows = _read_csv(manifest_path)
    required = {"dataset", "method", "result_path"}
    if not required.issubset(columns):
        raise ValueError(f"{manifest_path}: missing stream manifest columns")
    if len(rows) != 2 * len(methods) or len(set(row["dataset"] for row in rows)) != 2:
        raise ValueError(f"{manifest_path}: expected two complete stream datasets")
    expected_methods = set(methods)
    by_dataset = {}
    for row in rows:
        dataset = row["dataset"]
        method = row["method"]
        by_dataset.setdefault(dataset, set()).add(method)
        if method not in expected_methods:
            raise ValueError(f"{manifest_path}: unknown method {method}")
        for value in row.values():
            _check_logical_locator(value, f"{manifest_path}:{dataset}/{method}") if "/" in value else None
        image = root / dataset / f"{method}.png"
        if not image.is_file():
            raise ValueError(f"missing stream image: {image}")
    if any(methods_set != expected_methods for methods_set in by_dataset.values()):
        raise ValueError(f"{manifest_path}: incomplete method coverage")


def _validate_forebrain_cbdir(
    config: dict, result_root: Path, reference_neighbors_root: Path,
    reference_embeddings_root: Path, primary_rows: list,
) -> None:
    try:
        from .reference_embeddings import load_forebrain_frames
    except ImportError:
        from reference_embeddings import load_forebrain_frames
    frames = load_forebrain_frames(reference_neighbors_root, reference_embeddings_root)
    cbdir = config["datasets"]["forebrain"]["cbdir"]
    if (cbdir.get("embedding"), cbdir.get("dimensions"), cbdir.get("coordinate_sha256")) != (
        "pca", 30, frames["metadata"]["pca_coordinate_sha256"],
    ):
        raise ValueError("Forebrain CBDir config must bind the frozen PCA30 frame")
    manifest_path = result_root / "forebrain_cbdir_pca_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    edge_path = result_root / "forebrain_cbdir_pca_edge_details.csv"
    if _sha256(edge_path) != manifest["edge_details_sha256"]:
        raise ValueError("Forebrain PCA CBDir edge table SHA256 mismatch")
    if _sha256(result_root / "cell_metrics.csv") != manifest["primary_table_sha256"]:
        raise ValueError("Forebrain PCA CBDir primary table SHA256 mismatch")
    for name, digest in manifest["guard_tables"].items():
        if _sha256(result_root / name) != digest:
            raise ValueError(f"Forebrain promotion changed protected table {name}")
    _, rows = _read_csv(edge_path)
    expected = {
        (method, str(source), str(target)) for method in config["method_order"]
        for source, target in config["datasets"]["forebrain"]["cluster_edges"]
    }
    if len(rows) != 120 or {(r["method"], r["source"], r["target"]) for r in rows} != expected:
        raise ValueError("Forebrain PCA CBDir must contain exactly 20x6 edge rows")
    primary = {r["method"]: float(r["cbdir"]) for r in primary_rows if r["dataset"] == "forebrain"}
    for method in config["method_order"]:
        edges = [row for row in rows if row["method"] == method]
        values = []
        for row in edges:
            _nonnegative_int(row, "evaluable_source_cells", edge_path)
            if row["dataset"] != "forebrain":
                raise ValueError("Forebrain PCA edge table contains a different dataset")
            included = row["included_in_mean"] == "True"
            if row["score"]:
                _finite(row, "score", edge_path)
                if included != (float(row["score"]) != 0.0):
                    raise ValueError("Forebrain CBDir zero-edge exclusion differs from UnitVelo")
            elif included:
                raise ValueError("Forebrain CBDir cannot include an empty edge")
            if included:
                values.append(float(row["score"]))
        if not values or not np.isclose(np.mean(values), primary[method], rtol=0.0, atol=1e-12):
            raise ValueError(f"Forebrain PCA edge mean differs from CBDir for {method}")


def validate_results(
    config_path: Path = DEFAULT_CONFIG,
    result_root: Path = MODULE_DIR / "results",
    reference_neighbors_root: Path = MODULE_DIR / "reference_neighbors",
    streams_root: Path = MODULE_DIR / "streams",
    reference_embeddings_root: Path = MODULE_DIR / "reference_embeddings",
) -> None:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != 1:
        raise ValueError("Cell config has an unsupported schema version")
    expected_runtime = {"input_manifest": "inputs.tsv"}
    if config.get("runtime") != expected_runtime:
        raise ValueError("Cell config runtime paths differ from the module-local contract")
    datasets_config = config.get("datasets")
    methods = config.get("method_order")
    if not isinstance(datasets_config, dict) or not isinstance(methods, list):
        raise ValueError("Cell config is missing datasets or method_order")
    datasets = list(datasets_config)
    if len(datasets) != EXPECTED_DATASET_COUNT or set(datasets) != EXPECTED_DATASET_KEYS:
        raise ValueError("Cell config must contain exactly the frozen 14-dataset key set")
    if len(methods) != EXPECTED_METHOD_COUNT or set(methods) != EXPECTED_METHOD_KEYS:
        raise ValueError("Cell config must contain exactly the frozen 20-method key set")
    result_tables = config.get("result_tables")
    if not isinstance(result_tables, dict):
        raise ValueError("Cell config is missing result_tables")
    expected_table_columns = {
        "primary": EXPECTED_PRIMARY_COLUMNS,
        "supplementary": EXPECTED_SUPPLEMENTARY_COLUMNS,
        "audit": EXPECTED_AUDIT_COLUMNS,
    }
    for table_name, columns in expected_table_columns.items():
        if result_tables.get(table_name) != columns:
            raise ValueError(f"Cell config has an unexpected {table_name} table schema")
    reference = config.get("reference")
    if not isinstance(reference, dict):
        raise ValueError("Cell config is missing reference")
    expected_reference = {
        "n_neighbors": 30,
        "n_pcs": 30,
        "use_rep": "X_pca",
        "metric": "euclidean",
        "random_state": 0,
        "n_recurse_neighbors": 1,
        "mode_neighbors": "distances",
    }
    if any(reference.get(key) != value for key, value in expected_reference.items()):
        raise ValueError("Cell config differs from the fixed reference-neighbor contract")
    expected_keys = {(dataset, method) for dataset in datasets for method in methods}
    if len(expected_keys) != EXPECTED_DATASET_COUNT * EXPECTED_METHOD_COUNT:
        raise ValueError("Cell config does not define the frozen 14x20 key space")

    primary_path = result_root / "cell_metrics.csv"
    primary_rows = _check_table(
        primary_path, config["result_tables"]["primary"], expected_keys
    )
    for row in primary_rows:
        for column in config["result_tables"]["primary"][2:]:
            _finite(row, column, primary_path)

    supplementary_path = result_root / "cell_metrics_supplementary.csv"
    supplementary_rows = _check_table(
        supplementary_path,
        config["result_tables"]["supplementary"],
        expected_keys,
    )
    for row in supplementary_rows:
        _finite(row, "highdim_cbdir", supplementary_path)
        for column in config["result_tables"]["supplementary"]:
            if column.endswith("_edges") or column.endswith("_cells"):
                _nonnegative_int(row, column, supplementary_path)

    audit_path = result_root / "cell_metrics_audit.csv"
    audit_rows = _check_table(
        audit_path, config["result_tables"]["audit"], expected_keys
    )
    for row in audit_rows:
        _check_logical_locator(row["result_path"], f"{audit_path}:result_path")
        if not SHA256_RE.fullmatch(row["result_sha256"]):
            raise ValueError(f"{audit_path}: invalid result SHA256")
        if not row["result_path"].startswith("velocity:"):
            raise ValueError(f"{audit_path}: result_path must use velocity: locator")
        if not row["velocity_layer_key"] or row["status"] != "ok":
            raise ValueError(f"{audit_path}: invalid audit status/layer")
        _nonnegative_int(row, "invalid_velocity_values_replaced_with_zero", audit_path)
    _validate_audit_sources(audit_rows, expected_keys)

    edge_columns, edge_rows = _read_csv(result_root / "cell_metrics_edge_details.csv")
    required_edge_columns = {
        "dataset", "method", "metric", "source", "target", "score",
        "evaluable_source_cells",
    }
    if not required_edge_columns.issubset(edge_columns) or not edge_rows:
        raise ValueError("edge-detail table is missing or empty")
    parity_columns, parity_rows = _read_csv(result_root / "tivelo_parity.csv")
    if not {"dataset", "method", "metric", "raw_values_match", "aggregate_abs_difference"}.issubset(parity_columns):
        raise ValueError("TiVelo parity table is missing required columns")
    if len(parity_rows) != 2 * len(expected_keys):
        raise ValueError("TiVelo parity table has an unexpected row count")
    for row in parity_rows:
        if row["raw_values_match"] != "True" or float(row["aggregate_abs_difference"]) > 1e-7:
            raise ValueError("TiVelo parity check failed")

    _validate_reference_neighbors(
        reference_neighbors_root, datasets, reference
    )
    _validate_streams(streams_root, methods)
    _validate_forebrain_cbdir(
        config, result_root, reference_neighbors_root, reference_embeddings_root, primary_rows
    )

    print(
        f"Cell canonical bundle valid: {len(datasets)} datasets x "
        f"{len(methods)} methods = {len(expected_keys)} rows; "
        f"{len(datasets)} reference graphs; {2 * len(methods)} streams per dataset set"
    )


if __name__ == "__main__":
    validate_results()
