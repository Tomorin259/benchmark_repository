#!/usr/bin/env python3
"""Validate the Phase 3A cell-level configuration and result contract.

The validator deliberately uses only the Python standard library.  It checks
metadata and tabular contracts; it does not open, parse, or otherwise inspect
H5AD/AnnData files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "evaluation" / "cell" / "dataset_config.json"
DATASET_MANIFEST_PATH = ROOT / "manifests" / "datasets.tsv"
METHOD_MANIFEST_PATH = ROOT / "manifests" / "method_configurations.tsv"
INPUT_MANIFEST_PATH = ROOT / "manifests" / "cell_result_inputs.tsv"
RESULTS_PATH = ROOT / "results_summary" / "cell_metrics" / "all_metrics.csv"
PATH_ALIASES_PATH = ROOT / "manifests" / "path_aliases.tsv"

EXPECTED_ROW_COUNT = 280
EXPECTED_CONFIG_SCHEMA_VERSION = 1
EXPECTED_CONFIG_SOURCE_SHA256 = "c79929214882afc5d8dcb5e219b8a8f16785eb5a7938143f6e1e5c43d1885c50"
EXPECTED_METRIC_BASELINE_SHA256 = "f08e6f3f32c179541b80dda24c7cb2c7a0434701c929fa92f02c8073cb97faff"
EXPECTED_DATASET_SOURCE_KEYS = (
    "bonemarrow",
    "capofibro",
    "DentateGyrus_lamanno",
    "dentategyrus",
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
)
EXPECTED_DATASET_MANIFEST_IDS = {
    "bonemarrow": "bonemarrow",
    "capofibro": "capo_pre",
    "DentateGyrus_lamanno": "dentategyrus_lamanno",
    "dentategyrus": "dentategyrus",
    "forebrain": "forebrain",
    "gastrulation": "gastrulation",
    "hematopoiesis": "hematopoiesis",
    "mousebrain": "mousebrain",
    "neuron": "scnt_neuron",
    "organoid": "sceu_organoid",
    "pancreas": "pancreas",
    "retina": "retina",
    "rpe1": "sceu_rpe1",
    "rpe1_gen": "rpe1_gen",
}
EXPECTED_METHOD_IDS = (
    "scvelo_deterministic",
    "scvelo_stochastic",
    "scvelo_dynamic",
    "multivelo",
    "unitvelo",
    "unitvelo_mode2",
    "dynamo",
    "TIvelo",
    "cell2fate",
    "TFvelo",
    "veloae",
    "veloVI",
    "DeepKINET",
    "TSvelo",
    "celldancer",
    "deepvelo",
    "SDEvelo",
    "deepvelo2022",
    "graphvelo_scv",
    "graphvelo_dynamo",
)
EXPECTED_BATC_N_BINS = {
    "bonemarrow": 30,
    "capofibro": 5,
    "DentateGyrus_lamanno": 20,
    "dentategyrus": 30,
    "forebrain": 20,
    "gastrulation": 20,
    "hematopoiesis": 20,
    "mousebrain": 20,
    "neuron": 20,
    "organoid": 40,
    "pancreas": 30,
    "retina": 40,
    "rpe1": 10,
    "rpe1_gen": 20,
}
EXPECTED_LABELS = {
    **{dataset: "clusters" for dataset in EXPECTED_DATASET_SOURCE_KEYS},
    "organoid": "cell_type",
}
EXPECTED_SPECIFIC_LABEL_SOURCES = {
    "rpe1": "cell_cycle_phase",
    "retina": "Annotation",
    "neuron": "time_str",
    "hematopoiesis": "cell_type",
    "forebrain": "Clusters",
    "gastrulation": "celltype",
}
EXPECTED_METHOD_EMBEDDINGS = {
    method.lower(): ("umap_ori" if method.lower() == "dynamo" else "umap")
    for method in EXPECTED_METHOD_IDS
}
GENERIC_LABEL_FALLBACKS = ("celltype", "Clusters", "phase", "Annotation", "cell_type")
INPUT_COLUMNS = (
    "config_id",
    "operational_id",
    "dataset_id",
    "evaluation_dataset",
    "source_result_path",
    "result_path",
    "path_resolution",
    "velocity_layer_key",
    "source_embedding_key",
    "velocity_embedding_key",
    "evaluation_embedding",
    "alignment_mode",
    "n_reference",
    "n_method",
    "n_common",
    "source_embedding_max_error",
    "affine_relative_rmse",
    "source_nonfinite_velocity_values",
    "aligned_nonfinite_velocity_values",
    "file_size_bytes",
    "sha256",
)
RESULT_COLUMNS = ("dataset", "method", "consistency", "cbdir", "iccoh", "batc")
PATH_RESOLUTIONS = {"direct", "owner_confirmed_move"}
ALIGNMENT_MODES = {"direct", "affine", "graph_reproject"}
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
WINDOWS_ABSOLUTE_RE = re.compile(r"^[A-Za-z]:[\\/]")
PRIVATE_PATH_MARKERS = ("/mnt/second19T/", "/Users/", "\\Users\\")


class ValidationError(RuntimeError):
    """Raised for a contract violation."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def _read_tsv_columns(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    required = tuple(required)
    try:
        handle = path.open(newline="", encoding="utf-8-sig")
    except OSError as exc:
        raise ValidationError(f"{path}: cannot read manifest: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        missing = sorted(set(required).difference(columns))
        if missing:
            raise ValidationError(f"{path}: missing required columns: {missing}")
        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValidationError(f"{path}:{line_number}: non-rectangular row")
            rows.append({column: row[column] for column in required})
    return rows


def _read_table(path: Path, delimiter: str, expected_columns: tuple[str, ...]) -> list[dict[str, str]]:
    try:
        handle = path.open(newline="", encoding="utf-8-sig")
    except OSError as exc:
        raise ValidationError(f"{path}: cannot read table: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        columns = reader.fieldnames or []
        if columns != list(expected_columns):
            raise ValidationError(
                f"{path}: schema mismatch; expected {list(expected_columns)}, found {columns}"
            )
        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValidationError(f"{path}:{line_number}: non-rectangular row")
            empty = [column for column in expected_columns if row[column] == ""]
            if empty:
                raise ValidationError(f"{path}:{line_number}: empty fields: {empty}")
            rows.append(row)
    return rows


def _load_json(path: Path) -> dict:
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(f"{path}: invalid JSON: {exc}") from exc
    _require(isinstance(value, dict), f"{path}: top-level value must be an object")
    return value


def _is_non_bool_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _validate_edge_structure(dataset: str, edges: object) -> None:
    _require(isinstance(edges, list) and edges, f"dataset {dataset}: cluster_edges must be a non-empty list")
    seen: set[tuple[tuple[str, object], tuple[str, object]]] = set()
    for edge_index, edge in enumerate(edges, start=1):
        _require(
            isinstance(edge, list) and len(edge) == 2,
            f"dataset {dataset}: cluster_edges[{edge_index}] must be a two-item array",
        )
        normalized: list[tuple[str, object]] = []
        for endpoint_index, endpoint in enumerate(edge, start=1):
            valid = isinstance(endpoint, str) or _is_non_bool_int(endpoint)
            _require(
                valid,
                f"dataset {dataset}: cluster_edges[{edge_index}][{endpoint_index}] must be a string or integer",
            )
            if isinstance(endpoint, str):
                _require(
                    bool(endpoint.strip()),
                    f"dataset {dataset}: cluster_edges[{edge_index}] contains an empty endpoint",
                )
                normalized.append(("str", endpoint))
            else:
                normalized.append(("int", endpoint))
        edge_key = (normalized[0], normalized[1])
        _require(edge_key not in seen, f"dataset {dataset}: duplicate cluster edge {edge}")
        seen.add(edge_key)


def _validate_config(config: dict, dataset_manifest: list[dict[str, str]], method_manifest: list[dict[str, str]]) -> None:
    _require(
        config.get("schema_version") == EXPECTED_CONFIG_SCHEMA_VERSION,
        f"dataset_config.json: schema_version must be {EXPECTED_CONFIG_SCHEMA_VERSION}",
    )
    source = config.get("source")
    _require(isinstance(source, dict), "dataset_config.json: source must be an object")
    source_path = source.get("path")
    _require(source_path == "velocity:indicators/cell/all_metrics_main.py", "dataset_config.json: unexpected source path")
    _reject_private_path(source_path, "dataset_config.json source.path")
    _require(source.get("sha256") == EXPECTED_CONFIG_SOURCE_SHA256, "dataset_config.json: unexpected source sha256")
    _require(isinstance(source.get("symbols"), list) and source["symbols"], "dataset_config.json: source.symbols must be non-empty")

    _require(config.get("dataset_manifest") == "manifests/datasets.tsv", "dataset_config.json: unexpected dataset_manifest")
    _require(config.get("method_manifest") == "manifests/method_configurations.tsv", "dataset_config.json: unexpected method_manifest")
    _require(_is_non_bool_int(config.get("batc_default_n_bins")), "dataset_config.json: batc_default_n_bins must be an integer")
    _require(config["batc_default_n_bins"] > 0, "dataset_config.json: batc_default_n_bins must be positive")

    datasets = config.get("datasets")
    _require(isinstance(datasets, dict), "dataset_config.json: datasets must be an object")
    _require(set(datasets) == set(EXPECTED_DATASET_SOURCE_KEYS), "dataset_config.json: expected exactly the 14 source dataset keys")
    _require(len(datasets) == 14, "dataset_config.json: expected 14 datasets")

    dataset_rows_by_id = {row["dataset_id"]: row for row in dataset_manifest}
    _require(len(dataset_rows_by_id) == len(dataset_manifest), "manifests/datasets.tsv: duplicate dataset_id")
    for source_dataset in EXPECTED_DATASET_SOURCE_KEYS:
        entry = datasets[source_dataset]
        _require(isinstance(entry, dict), f"dataset {source_dataset}: configuration must be an object")
        expected_manifest_id = EXPECTED_DATASET_MANIFEST_IDS[source_dataset]
        _require(entry.get("dataset_id") == expected_manifest_id, f"dataset {source_dataset}: dataset_id must be {expected_manifest_id!r}")
        _require(entry["dataset_id"] in dataset_rows_by_id, f"dataset {source_dataset}: dataset_id is not a foreign key in datasets.tsv")
        _require(
            dataset_rows_by_id[entry["dataset_id"]]["role"] in {"real_cell", "real_cell_and_kinetic"},
            f"dataset {source_dataset}: dataset_id does not refer to a cell-level dataset",
        )

        _require(isinstance(entry.get("emb"), str) and bool(entry["emb"].strip()), f"dataset {source_dataset}: emb must be a non-empty string")
        _require(entry["emb"] == "umap", f"dataset {source_dataset}: emb must be the extracted default 'umap'")
        _require(entry.get("label") == EXPECTED_LABELS[source_dataset], f"dataset {source_dataset}: unexpected label")
        _validate_edge_structure(source_dataset, entry.get("cluster_edges"))
        _require(entry.get("batc_n_bins") == EXPECTED_BATC_N_BINS[source_dataset], f"dataset {source_dataset}: unexpected BATC n_bins")
        expected_n_bins_source = "DEFAULT_BATC_N_BINS" if source_dataset == "rpe1_gen" else "BATC_N_BINS_BY_DATASET"
        _require(entry.get("batc_n_bins_source") == expected_n_bins_source, f"dataset {source_dataset}: unexpected BATC n_bins source")

        label_source = entry.get("label_source")
        _require(isinstance(label_source, dict), f"dataset {source_dataset}: label_source must be an object")
        _require(label_source.get("target") == entry["label"], f"dataset {source_dataset}: label_source.target must equal label")
        fallback_sources = label_source.get("fallback_sources")
        _require(isinstance(fallback_sources, list), f"dataset {source_dataset}: label_source.fallback_sources must be a list")
        if source_dataset == "organoid":
            _require(label_source.get("strategy") == "preexisting_target_only", "dataset organoid: label source must be preexisting_target_only")
            _require(label_source.get("dataset_specific_source") is None, "dataset organoid: no source column is applied")
            _require(fallback_sources == [], "dataset organoid: generic fallback must not be claimed")
            mapper = label_source.get("legacy_mapper")
            _require(isinstance(mapper, dict), "dataset organoid: legacy_mapper must be recorded")
            _require(mapper.get("name") == "CELL_MAPPER_ORGANOID", "dataset organoid: unexpected legacy mapper name")
            _require(mapper.get("source_column") == "som_cluster_id", "dataset organoid: unexpected legacy mapper source column")
            _require(mapper.get("target_column") == "cell_type", "dataset organoid: unexpected legacy mapper target column")
            _require(mapper.get("status") == "legacy_defined_but_not_applied", "dataset organoid: mapper status must be legacy_defined_but_not_applied")
            _require(mapper.get("legacy_defined_but_not_applied") is True, "dataset organoid: legacy mapper flag must be true")
            _require(mapper.get("applies_at_runtime") is False, "dataset organoid: mapper must not be claimed as runtime-applied")
            mapping = mapper.get("mapping")
            _require(isinstance(mapping, dict) and set(mapping) == {str(i) for i in range(1, 12)}, "dataset organoid: legacy mapper mapping must contain keys 1-11")
        else:
            expected_specific = EXPECTED_SPECIFIC_LABEL_SOURCES.get(source_dataset)
            _require(fallback_sources == list(GENERIC_LABEL_FALLBACKS), f"dataset {source_dataset}: generic fallback order differs from prepare_clusters")
            _require(label_source.get("dataset_specific_source") == expected_specific, f"dataset {source_dataset}: unexpected dataset-specific label source")
            expected_strategy = "dataset_specific_then_generic_fallback" if expected_specific else "generic_fallback_if_target_missing"
            _require(label_source.get("strategy") == expected_strategy, f"dataset {source_dataset}: unexpected label-source strategy")
            _require("legacy_mapper" not in label_source, f"dataset {source_dataset}: legacy organoid mapper must not be attached")

    method_ids_from_manifest = [row["operational_id"] for row in method_manifest]
    _require(len(method_ids_from_manifest) == len(set(method_ids_from_manifest)), "manifests/method_configurations.tsv: duplicate operational_id")
    _require(set(method_ids_from_manifest) == set(EXPECTED_METHOD_IDS), "manifests/method_configurations.tsv: operational IDs differ from the 20-method contract")
    overrides = config.get("method_embedding_overrides")
    _require(isinstance(overrides, dict), "dataset_config.json: method_embedding_overrides must be an object")
    normalized_overrides = {str(method).lower(): value for method, value in overrides.items()}
    _require(len(normalized_overrides) == 20, "dataset_config.json: expected 20 method embedding overrides")
    _require(set(normalized_overrides) == set(EXPECTED_METHOD_EMBEDDINGS), "dataset_config.json: method override keys differ from the 20 operational methods")
    for method, expected_embedding in EXPECTED_METHOD_EMBEDDINGS.items():
        _require(normalized_overrides[method] == expected_embedding, f"method {method}: unexpected embedding override")


def _reject_private_path(value: object, context: str) -> None:
    _require(isinstance(value, str), f"{context}: path must be a string")
    _require(not value.startswith("/"), f"{context}: absolute paths are not allowed")
    _require(not WINDOWS_ABSOLUTE_RE.match(value), f"{context}: Windows absolute paths are not allowed")
    _require(not any(marker in value for marker in PRIVATE_PATH_MARKERS), f"{context}: private path is not allowed")


def _parse_non_negative_float(value: str, context: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValidationError(f"{context}: expected a numeric value") from exc
    _require(math.isfinite(parsed) and parsed >= 0, f"{context}: expected a finite non-negative value")
    return parsed


def _validate_manifest_row(
    row: dict[str, str],
    line_number: int,
    config_by_id: dict[str, str],
    dataset_by_evaluation: dict[str, str],
    embedding_by_method: dict[str, str],
    render_aliases: dict[str, str],
) -> None:
    context = f"cell_result_inputs.tsv:{line_number}"
    _require(row["config_id"] in config_by_id, f"{context}: unknown config_id {row['config_id']!r}")
    _require(
        config_by_id[row["config_id"]] == row["operational_id"],
        f"{context}: operational_id does not match config_id",
    )
    _require(
        dataset_by_evaluation.get(row["evaluation_dataset"]) == row["dataset_id"],
        f"{context}: dataset_id does not match evaluation_dataset",
    )
    _reject_private_path(row["source_result_path"], context + " source_result_path")
    _reject_private_path(row["result_path"], context + " result_path")
    _require(row["path_resolution"] in PATH_RESOLUTIONS, f"{context}: path_resolution must be one of {sorted(PATH_RESOLUTIONS)}")
    _require(row["source_result_path"].startswith("velocity:"), f"{context}: source_result_path must use velocity:")
    _require(row["result_path"].startswith("velocity:"), f"{context}: result_path must use velocity:")
    source = _normalise_locator(row["source_result_path"])
    result = _normalise_locator(row["result_path"])
    if row["path_resolution"] == "direct":
        _require(source == result, f"{context}: direct path must be unchanged")
    else:
        _require(source != result, f"{context}: moved path must differ from source")
        _require(render_aliases.get(source) == result, f"{context}: moved path is not supported by path_aliases.tsv")

    size_text = row["file_size_bytes"]
    _require(size_text.isdigit() and int(size_text) > 0, f"{context}: file_size_bytes must be a positive integer")
    _require(SHA256_RE.fullmatch(row["sha256"]) is not None, f"{context}: sha256 must be 64 hexadecimal characters")
    for field in ("velocity_layer_key", "source_embedding_key", "velocity_embedding_key"):
        _require(bool(row[field].strip()), f"{context}: {field} must be non-empty")
    expected_embedding = embedding_by_method[row["operational_id"].lower()]
    _require(row["evaluation_embedding"] == expected_embedding, f"{context}: unexpected evaluation_embedding")
    _require(row["alignment_mode"] in ALIGNMENT_MODES, f"{context}: unsupported alignment_mode")
    for field in ("n_reference", "n_method", "n_common"):
        _require(row[field].isdigit() and int(row[field]) > 0, f"{context}: {field} must be a positive integer")
    _require(int(row["n_common"]) <= int(row["n_method"]), f"{context}: n_common exceeds n_method")
    for field in (
        "source_embedding_max_error",
        "affine_relative_rmse",
        "source_nonfinite_velocity_values",
        "aligned_nonfinite_velocity_values",
    ):
        _parse_non_negative_float(row[field], f"{context}: {field}")


def _validate_cell_result_inputs(
    config: dict, method_manifest: list[dict[str, str]]
) -> tuple[list[dict[str, str]], set[tuple[str, str]]]:
    rows = _read_table(INPUT_MANIFEST_PATH, "\t", INPUT_COLUMNS)
    _require(len(rows) == EXPECTED_ROW_COUNT, f"{INPUT_MANIFEST_PATH}: expected {EXPECTED_ROW_COUNT} rows, found {len(rows)}")
    config_by_id = {row["config_id"]: row["operational_id"] for row in method_manifest}
    _require(len(config_by_id) == 20, "method configuration manifest must contain 20 unique config IDs")
    dataset_by_evaluation = {
        evaluation: entry["dataset_id"] for evaluation, entry in config["datasets"].items()
    }
    embedding_by_method = {
        method.lower(): embedding for method, embedding in config["method_embedding_overrides"].items()
    }
    render_aliases = _load_render_aliases()
    contract_keys: set[tuple[str, str]] = set()
    metric_keys: set[tuple[str, str]] = set()
    canonical_paths: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        _validate_manifest_row(
            row,
            line_number,
            config_by_id,
            dataset_by_evaluation,
            embedding_by_method,
            render_aliases,
        )
        contract_key = (row["config_id"], row["dataset_id"])
        _require(contract_key not in contract_keys, f"cell_result_inputs.tsv:{line_number}: duplicate contract key {contract_key}")
        contract_keys.add(contract_key)
        metric_key = (row["evaluation_dataset"], row["operational_id"])
        _require(metric_key not in metric_keys, f"cell_result_inputs.tsv:{line_number}: duplicate metric key {metric_key}")
        metric_keys.add(metric_key)
        canonical_path = _normalise_locator(row["result_path"])
        _require(canonical_path not in canonical_paths, f"cell_result_inputs.tsv:{line_number}: duplicate canonical path {canonical_path}")
        canonical_paths.add(canonical_path)
    expected_contract_keys = {
        (config_id, dataset_id)
        for config_id in config_by_id
        for dataset_id in dataset_by_evaluation.values()
    }
    _require(
        contract_keys == expected_contract_keys,
        f"{INPUT_MANIFEST_PATH}: key set is not the complete 20x14 product",
    )
    resolution_counts = {
        value: sum(row["path_resolution"] == value for row in rows)
        for value in PATH_RESOLUTIONS
    }
    _require(
        resolution_counts == {"direct": 268, "owner_confirmed_move": 12},
        f"{INPUT_MANIFEST_PATH}: unexpected path-resolution counts {resolution_counts}",
    )
    return rows, metric_keys


def _validate_results(keys: set[tuple[str, str]]) -> None:
    _require(
        _sha256(RESULTS_PATH) == EXPECTED_METRIC_BASELINE_SHA256,
        f"{RESULTS_PATH}: file bytes differ from the frozen metric baseline",
    )
    rows = _read_table(RESULTS_PATH, ",", RESULT_COLUMNS)
    _require(len(rows) == EXPECTED_ROW_COUNT, f"{RESULTS_PATH}: expected {EXPECTED_ROW_COUNT} rows, found {len(rows)}")
    result_keys: set[tuple[str, str]] = set()
    for line_number, row in enumerate(rows, start=2):
        key = (row["dataset"], row["method"])
        _require(key not in result_keys, f"all_metrics.csv:{line_number}: duplicate unique key {key}")
        result_keys.add(key)
        for metric in ("consistency", "cbdir", "iccoh", "batc"):
            try:
                value = float(row[metric])
            except ValueError as exc:
                raise ValidationError(f"all_metrics.csv:{line_number}: {metric} is not numeric") from exc
            _require(math.isfinite(value), f"all_metrics.csv:{line_number}: {metric} must be finite")
    _require(result_keys == keys, f"{RESULTS_PATH}: dataset/method key set differs from cell_result_inputs.tsv; missing={sorted(keys - result_keys)}, extra={sorted(result_keys - keys)}")


def _normalise_locator(value: str) -> str:
    """Return a path relative to the audited project root."""
    if value.startswith("velocity:"):
        value = "velocity/" + value[len("velocity:") :].lstrip("/")
    return value


def _resolve_under_root(value: str, working_root: Path) -> Path:
    _reject_private_path(value, "result path")
    relative = _normalise_locator(value)
    _require(not relative.startswith(("project:", "paper:")), f"result_path uses unsupported locator: {value!r}")
    relative_path = Path(relative)
    _require(not relative_path.is_absolute(), f"result_path must be relative: {value!r}")
    root = working_root.resolve()
    candidate = (root / relative_path).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ValidationError(f"result_path escapes working root: {value!r}") from exc
    return candidate


def _load_render_aliases() -> dict[str, str]:
    if not PATH_ALIASES_PATH.is_file():
        return {}
    aliases: dict[str, str] = {}
    for row in _read_tsv_columns(PATH_ALIASES_PATH, ("artifact_type", "historical_path", "current_path")):
        if row["artifact_type"] == "render":
            aliases[_normalise_locator(row["historical_path"])] = _normalise_locator(
                row["current_path"]
            )
    return aliases


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _check_files(rows: list[dict[str, str]], working_root: Path, verify_hashes: bool) -> None:
    _require(working_root.is_dir(), f"--working-root is not a directory: {working_root}")
    for line_number, row in enumerate(rows, start=2):
        candidate = _resolve_under_root(row["result_path"], working_root)
        _require(candidate.is_file(), f"cell_result_inputs.tsv:{line_number}: result file does not exist: {candidate}")
        expected_size = int(row["file_size_bytes"])
        actual_size = candidate.stat().st_size
        _require(actual_size == expected_size, f"cell_result_inputs.tsv:{line_number}: size mismatch for {candidate}: expected {expected_size}, found {actual_size}")
        if verify_hashes:
            expected_hash = row["sha256"].lower()
            actual_hash = _sha256(candidate)
            _require(actual_hash == expected_hash, f"cell_result_inputs.tsv:{line_number}: sha256 mismatch for {candidate}")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate the Phase 3A cell-level contract using only the standard library.")
    parser.add_argument(
        "--working-root",
        type=Path,
        default=None,
        help="Optional local root used to resolve result_path values and check files/sizes.",
    )
    parser.add_argument(
        "--verify-hashes",
        action="store_true",
        help="With --working-root, compute SHA-256 for every listed result file.",
    )
    parser.add_argument(
        "--require-artifacts",
        action="store_true",
        help="Require the frozen cell_result_inputs.tsv and all_metrics.csv files to exist.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        if args.verify_hashes:
            _require(args.working_root is not None, "--verify-hashes requires --working-root")
        config = _load_json(CONFIG_PATH)
        dataset_manifest = _read_tsv_columns(DATASET_MANIFEST_PATH, ("dataset_id", "role"))
        method_manifest = _read_tsv_columns(
            METHOD_MANIFEST_PATH, ("config_id", "operational_id")
        )
        _validate_config(config, dataset_manifest, method_manifest)

        input_present = INPUT_MANIFEST_PATH.is_file()
        results_present = RESULTS_PATH.is_file()
        artifact_present = input_present or results_present
        if args.require_artifacts or artifact_present:
            _require(input_present, f"missing required artifact: {INPUT_MANIFEST_PATH}")
            _require(results_present, f"missing required artifact: {RESULTS_PATH}")
            input_rows, keys = _validate_cell_result_inputs(config, method_manifest)
            _validate_results(keys)
            if args.working_root is not None:
                _check_files(input_rows, args.working_root, args.verify_hashes)
        elif args.working_root is not None:
            _require(args.working_root.is_dir(), f"--working-root is not a directory: {args.working_root}")
            print("SKIP: cell result artifacts are not present; no file checks were run.")

        print("Cell Phase 3A contract validation passed.")
        return 0
    except (ValidationError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
