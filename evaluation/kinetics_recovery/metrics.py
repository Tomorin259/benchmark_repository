#!/usr/bin/env python3
"""Evaluate kinetic-parameter recovery from saved method-result H5AD files.

This module owns the kinetic experiment from method-result H5AD input through
parameter extraction, dataset-specific trajectory processing, numerical
scoring, and compact CSV output. It does not prepare datasets or train
methods. Method-specific H5AD layouts are handled explicitly at this input
boundary because they are part of the scientific result contract.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MODULE_DIR / "config.json"


def source_path(artifact_root: Path, relative_path: str) -> Path:
    """Resolve one config path while keeping the public config relocatable."""
    path = Path(relative_path)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Config paths must be relative to artifact root: {relative_path}")
    resolved = artifact_root / path
    if not resolved.is_file():
        raise FileNotFoundError(f"Required source does not exist: {resolved}")
    return resolved


def read_rpe1_reference(path: Path, parameter: str) -> pd.DataFrame:
    table = pd.read_csv(path, skiprows=1).set_index("gene_symbol")
    prefix = "norm_kappa_" if parameter == "alpha" else "norm_gamma_"
    columns = [name for name in table.columns if name.startswith(prefix)]
    values = table.loc[:, columns]
    values.columns = range(values.shape[1])
    return values.apply(pd.to_numeric, errors="coerce")


def read_organoid_reference(path: Path, branch: str, parameter: str) -> pd.DataFrame:
    table = pd.read_csv(path, skiprows=1).set_index("gene_symbol")
    stem = "kappa" if parameter == "alpha" else "gamma"
    columns = sorted(
        name for name in table.columns if f"norm_{stem}_branch_{branch}" in name
    )
    values = table.loc[:, columns]
    values.columns = range(values.shape[1])
    return values.apply(pd.to_numeric, errors="coerce")


def calculate_metrics(
    prediction: pd.DataFrame, reference: pd.DataFrame
) -> pd.DataFrame:
    """Calculate per-gene Spearman and cosine scores after alignment."""
    genes = reference.index.intersection(prediction.index, sort=False)
    timepoints = reference.columns.intersection(prediction.columns, sort=False)
    prediction = prediction.loc[genes, timepoints]
    reference = reference.loc[genes, timepoints]

    rows = []
    for gene in genes:
        predicted = prediction.loc[gene].to_numpy(dtype=float)
        observed = reference.loc[gene].to_numpy(dtype=float)
        valid = np.isfinite(predicted) & np.isfinite(observed)
        count = int(valid.sum())
        if count < 2:
            correlation = np.nan
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                correlation = spearmanr(predicted[valid], observed[valid]).statistic
        if count == 0:
            cosine = np.nan
        else:
            denominator = np.linalg.norm(predicted[valid]) * np.linalg.norm(observed[valid])
            cosine = (
                np.nan
                if denominator == 0
                else np.dot(predicted[valid], observed[valid]) / denominator
            )
        rows.append(
            {
                "Gene": gene,
                "Spearman_Corr": correlation,
                "Cosine_Similarity": cosine,
                "Valid_Samples": count,
            }
        )
    return pd.DataFrame(rows)


def validate_standard_result(
    adata: ad.AnnData, dataset: str, settings: dict, source: dict
) -> None:
    """Validate only the fields consumed by direct layer extraction."""
    expected_cells = settings[dataset]["n_obs"]
    if adata.n_obs != expected_cells:
        raise ValueError(
            f"{dataset}: expected {expected_cells} cells, found {adata.n_obs}"
        )
    if adata.obs_names.has_duplicates or adata.var_names.has_duplicates:
        raise ValueError(f"{dataset}: obs_names and var_names must be unique")

    if dataset == "rpe1":
        required_obs = [settings["rpe1"]["obs_key"]]
    else:
        required_obs = [
            settings["organoid"]["branch_key"],
            settings["organoid"]["cluster_key"],
            settings["organoid"]["pseudotime_key"],
        ]
    missing_obs = [name for name in required_obs if name not in adata.obs]
    if missing_obs:
        raise KeyError(f"{dataset}: missing obs fields {missing_obs}")

    required_layers = list(source.get("layers", {}).values()) + source.get(
        "required_layers", []
    )
    missing_layers = [name for name in required_layers if name not in adata.layers]
    if missing_layers:
        raise KeyError(f"{dataset}: missing parameter layers {missing_layers}")


def layer_frame(adata: ad.AnnData, layer: str) -> pd.DataFrame:
    values = adata.layers[layer]
    if sparse.issparse(values):
        values = values.toarray()
    return pd.DataFrame(
        np.asarray(values), index=adata.obs_names, columns=adata.var_names
    )


def rpe1_raw_from_layer(
    adata: ad.AnnData, layer: str, settings: dict
) -> pd.DataFrame:
    """Average a cell-by-gene layer at each fixed RPE1 cycle position."""
    values = layer_frame(adata, layer)
    positions = pd.to_numeric(adata.obs[settings["obs_key"]], errors="raise")
    return values.groupby(positions).mean().reindex(range(settings["positions"])).T


def process_rpe1(raw: pd.DataFrame, settings: dict, clip_negative: bool) -> pd.DataFrame:
    """Apply RPE1 smoothing, median normalization, and log2 transformation."""
    if clip_negative:
        raw = raw.clip(lower=0)
    original_nan = raw.isna()
    smoothed = raw.T.rolling(
        settings["smooth_window"], center=True, min_periods=1
    ).mean().T
    medians = smoothed.median(axis=1).replace(0, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        processed = np.log2(smoothed.div(medians, axis=0) + settings["epsilon"])
    return processed.mask(original_nan)


def select_organoid_branch(
    adata: ad.AnnData, settings: dict, branch: str
) -> tuple[ad.AnnData, pd.Series]:
    branch_settings = settings["branches"][branch]
    branch_mask = adata.obs[settings["branch_key"]].astype(str).str.contains(branch)
    cluster_mask = adata.obs[settings["cluster_key"]].astype(str).isin(
        branch_settings["clusters"]
    )
    selected = adata[branch_mask & cluster_mask].copy()
    pseudotime = pd.to_numeric(
        selected.obs[settings["pseudotime_key"]], errors="raise"
    )
    if branch_settings["reverse_time"]:
        aligned_time = branch_settings["start"] - pseudotime
    else:
        aligned_time = pseudotime - branch_settings["start"]
    selected = selected[aligned_time >= 0].copy()
    return selected, aligned_time.loc[selected.obs_names]


def organoid_raw_from_layer(
    adata: ad.AnnData, aligned_time: pd.Series, layer: str, settings: dict
) -> pd.DataFrame:
    """Average a cell-by-gene layer in the fixed organoid time bins."""
    values = layer_frame(adata, layer)
    bins = pd.cut(
        aligned_time, bins=settings["bins"], labels=False, include_lowest=True
    )
    return values.groupby(bins).mean().reindex(range(settings["bins"])).T


def process_organoid(raw: pd.DataFrame, settings: dict) -> pd.DataFrame:
    """Apply organoid smoothing and early-bin baseline normalization."""
    smoothed = raw.T.rolling(
        settings["smooth_window"], center=True, min_periods=1
    ).mean().T
    baseline = smoothed.iloc[:, : settings["baseline_bins"]].mean(axis=1)
    return smoothed.div(baseline + settings["epsilon"], axis=0)


def dynamo_rpe1_raw(
    source: dict, reference: pd.DataFrame, settings: dict, artifact_root: Path
) -> pd.DataFrame:
    """Extract Dynamo alpha from its labelled sparse RPE1 H5AD."""
    path = source_path(artifact_root, source["result"])
    with h5py.File(path, "r") as handle:
        layer_group = handle["layers"][source["layer"]]
        if layer_group.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{path}: layers/{source['layer']} is not stored as CSR")
        shape = tuple(int(value) for value in layer_group.attrs["shape"])
        matrix = sparse.csr_matrix(
            (
                np.asarray(layer_group["data"]),
                np.asarray(layer_group["indices"]),
                np.asarray(layer_group["indptr"]),
            ),
            shape=shape,
        )
        var_names = pd.Index(
            [
                value.decode("utf-8") if isinstance(value, bytes) else str(value)
                for value in np.asarray(handle["var"]["_index"])
            ]
        )
        positions = pd.Series(
            np.asarray(handle["obs"][source["position_key"]]),
            dtype=float,
            name=source["position_key"],
        )
    if matrix.shape[1] != len(var_names):
        raise ValueError(f"{path}: layer and var index have incompatible dimensions")
    if matrix.shape[0] != settings["n_obs"]:
        raise ValueError(
            f"Dynamo RPE1: expected {settings['n_obs']} cells, found {matrix.shape[0]}"
        )
    if positions.shape[0] != matrix.shape[0]:
        raise ValueError("Dynamo RPE1: position vector does not match cell count")
    if not np.isfinite(positions).all():
        raise ValueError("Dynamo RPE1: cycle positions contain non-finite values")

    genes = var_names.intersection(reference.index, sort=False)
    column_indices = var_names.get_indexer(genes)
    values = pd.DataFrame(
        matrix[:, column_indices].toarray(), columns=genes, dtype=float
    )
    bins = pd.cut(
        positions,
        bins=range(settings["positions"] + 1),
        labels=False,
        include_lowest=True,
    )
    raw = values.groupby(bins, observed=True).mean().T
    # Dynamo's labelled alpha export contains genes with no observed alpha.
    return raw.loc[~(raw == 0).all(axis=1)]


def tsvelo_rpe1_raw(
    source: dict, reference: pd.DataFrame, settings: dict, artifact_root: Path
) -> pd.DataFrame:
    """Extract and merge alpha from the two lineage-specific TSvelo H5ADs."""
    alpha_frames = []
    positions = []
    for relative_path in source["results"]:
        path = source_path(artifact_root, relative_path)
        adata = ad.read_h5ad(path)
        if adata.obs_names.has_duplicates or adata.var_names.has_duplicates:
            raise ValueError(f"{path}: obs_names and var_names must be unique")
        if source["layer"] not in adata.layers:
            raise KeyError(f"{path}: missing layers/{source['layer']}")
        if source["position_key"] not in adata.obs:
            raise KeyError(f"{path}: missing obs/{source['position_key']}")
        alpha_frames.append(layer_frame(adata, source["layer"]))
        positions.append(
            pd.to_numeric(
                adata.obs[source["position_key"]], errors="raise"
            ).rename(source["position_key"])
        )

    combined = pd.concat(alpha_frames)
    # The two lineage files overlap in cells; the benchmark uses their mean.
    merged_alpha = combined.groupby(level=0).mean()
    merged_positions = pd.concat(positions)
    merged_positions = merged_positions[
        ~merged_positions.index.duplicated(keep="first")
    ]
    merged = merged_alpha.join(merged_positions, how="left")
    if merged[source["position_key"]].isna().any():
        raise ValueError("TSvelo RPE1: a merged cell is missing its cycle position")

    binned = merged.groupby(source["position_key"]).mean()
    binned = binned.reindex(pd.Index(range(settings["positions"]), dtype=float))
    raw = binned.T
    genes = reference.index.intersection(raw.index, sort=False)
    return raw.loc[genes]


def score_component(
    method: str,
    dataset: str,
    branch: str | None,
    parameter: str,
    prediction: pd.DataFrame,
    reference: pd.DataFrame,
) -> tuple[dict, pd.DataFrame]:
    metrics = calculate_metrics(prediction, reference)
    metrics.insert(0, "Parameter", parameter)
    metrics.insert(0, "Branch", branch if branch is not None else pd.NA)
    metrics.insert(0, "Dataset", dataset)
    metrics.insert(0, "Method", method)
    component = {
        "Method": method,
        "Dataset": dataset,
        "Branch": branch if branch is not None else pd.NA,
        "Parameter": parameter,
        "Component": (
            f"rpe1_{parameter}"
            if dataset == "rpe1"
            else f"organoid_branch{branch}_{parameter}"
        ),
        "Mean": metrics["Spearman_Corr"].mean(),
        "N": metrics["Spearman_Corr"].notna().sum(),
    }
    return component, metrics


def save_parameter_matrices(
    output_dir: Path | None,
    dataset: str,
    branch: str | None,
    method: str,
    parameter: str,
    raw: pd.DataFrame,
    processed: pd.DataFrame,
) -> None:
    if output_dir is None:
        return
    matrix_directory = output_dir / "matrices" / dataset
    trajectory_directory = output_dir / "trajectories" / dataset
    if branch is not None:
        matrix_directory /= f"branch{branch}"
        trajectory_directory /= f"branch{branch}"
    matrix_directory.mkdir(parents=True, exist_ok=True)
    trajectory_directory.mkdir(parents=True, exist_ok=True)
    raw.to_csv(matrix_directory / f"{method}_{parameter}_raw.csv")
    processed.to_csv(trajectory_directory / f"{method}_{parameter}.csv")


def evaluate_rpe1(
    method: str,
    source: dict,
    settings: dict,
    references: dict[str, pd.DataFrame],
    artifact_root: Path,
    output_dir: Path | None,
) -> tuple[list[dict], list[pd.DataFrame]]:
    """Extract and score the configured RPE1 parameters for one method."""
    reader = source["reader"]
    if reader == "layers":
        adata = ad.read_h5ad(source_path(artifact_root, source["result"]))
        validate_standard_result(adata, "rpe1", settings, source)
        raw_by_parameter = {
            parameter: rpe1_raw_from_layer(adata, layer, settings["rpe1"])
            for parameter, layer in source["layers"].items()
        }
    elif reader == "dynamo_csr":
        raw_by_parameter = {
            "alpha": dynamo_rpe1_raw(
                source, references["alpha"], settings["rpe1"], artifact_root
            )
        }
    elif reader == "tsvelo_lineages":
        raw_by_parameter = {
            "alpha": tsvelo_rpe1_raw(
                source, references["alpha"], settings["rpe1"], artifact_root
            )
        }
    else:
        raise ValueError(f"{method} RPE1: unsupported reader {reader!r}")

    components = []
    gene_metrics = []
    for parameter, raw in raw_by_parameter.items():
        processed = process_rpe1(
            raw, settings["rpe1"], source.get("clip_negative", False)
        )
        component, metrics = score_component(
            method, "rpe1", None, parameter, processed, references[parameter]
        )
        components.append(component)
        gene_metrics.append(metrics)
        save_parameter_matrices(
            output_dir, "rpe1", None, method, parameter, raw, processed
        )
    return components, gene_metrics


def evaluate_organoid(
    method: str,
    source: dict,
    settings: dict,
    references: dict[str, dict[str, pd.DataFrame]],
    artifact_root: Path,
    output_dir: Path | None,
) -> tuple[list[dict], list[pd.DataFrame]]:
    if source["reader"] != "layers":
        raise ValueError(f"{method} organoid: unsupported reader {source['reader']!r}")
    adata = ad.read_h5ad(source_path(artifact_root, source["result"]))
    validate_standard_result(adata, "organoid", settings, source)

    components = []
    gene_metrics = []
    for branch in settings["organoid"]["branches"]:
        selected, aligned_time = select_organoid_branch(
            adata, settings["organoid"], branch
        )
        for parameter, layer in source["layers"].items():
            raw = organoid_raw_from_layer(
                selected, aligned_time, layer, settings["organoid"]
            )
            processed = process_organoid(raw, settings["organoid"])
            component, metrics = score_component(
                method,
                "organoid",
                branch,
                parameter,
                processed,
                references[branch][parameter],
            )
            components.append(component)
            gene_metrics.append(metrics)
            save_parameter_matrices(
                output_dir, "organoid", branch, method, parameter, raw, processed
            )
    return components, gene_metrics


def summarize(components: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    system_rows = []
    for (method, parameter, dataset), group in components.groupby(
        ["Method", "Parameter", "Dataset"], dropna=False
    ):
        system_rows.append(
            {
                "Method": method,
                "Parameter": parameter,
                "System": dataset,
                "Mean": group["Mean"].mean(),
                "N_components": len(group),
                "N_genes": group["N"].sum(),
            }
        )
    systems = pd.DataFrame(system_rows)
    methods = (
        systems.groupby(["Method", "Parameter"], as_index=False)
        .agg(
            Mean=("Mean", "mean"),
            N_systems=("System", "count"),
            N_genes=("N_genes", "sum"),
        )
        .sort_values(["Parameter", "Method"])
    )
    return systems.sort_values(["Parameter", "Method", "System"]), methods


def ensure_output(path: Path) -> None:
    if path.exists() and not path.is_dir():
        raise NotADirectoryError(path)
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Output directory must be empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def load_references(
    artifact_root: Path, settings: dict
) -> tuple[dict[str, pd.DataFrame], dict[str, dict[str, pd.DataFrame]]]:
    rpe1_reference_path = source_path(artifact_root, settings["rpe1"]["reference"])
    full_rpe1_alpha = read_rpe1_reference(rpe1_reference_path, "alpha")
    alpha_gene_path = source_path(
        artifact_root, settings["rpe1"]["alpha_gene_subset"]
    )
    alpha_gene_subset = pd.Index(pd.read_csv(alpha_gene_path, index_col=0).index)
    if alpha_gene_subset.has_duplicates:
        raise ValueError(f"{alpha_gene_path}: the gene subset contains duplicate names")
    missing_alpha_genes = alpha_gene_subset.difference(full_rpe1_alpha.index)
    if not missing_alpha_genes.empty:
        raise ValueError(
            "RPE1 alpha gene subset is not contained in the reference table: "
            f"{', '.join(missing_alpha_genes[:5])}"
        )
    rpe1_references = {
        "alpha": full_rpe1_alpha.loc[alpha_gene_subset],
        "gamma": read_rpe1_reference(rpe1_reference_path, "gamma"),
    }

    organoid_reference_path = source_path(
        artifact_root, settings["organoid"]["reference"]
    )
    organoid_references = {
        branch: {
            parameter: read_organoid_reference(
                organoid_reference_path, branch, parameter
            )
            for parameter in ("alpha", "gamma")
        }
        for branch in settings["organoid"]["branches"]
    }
    return rpe1_references, organoid_references


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact-root",
        type=Path,
        required=True,
        help="Root containing the relative reference and method-result paths in config.json.",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--methods", nargs="*")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.check_only == (args.output is not None):
        parser.error("use exactly one of --output or --check-only")

    with args.config.open(encoding="utf-8") as handle:
        settings = json.load(handle)
    methods = args.methods or list(settings["methods"])
    unknown = [method for method in methods if method not in settings["methods"]]
    if unknown:
        parser.error(f"unknown methods: {', '.join(unknown)}")

    artifact_root = args.artifact_root.resolve()
    output = None if args.check_only else args.output.resolve()
    if output is not None:
        ensure_output(output)

    rpe1_references, organoid_references = load_references(
        artifact_root, settings["datasets"]
    )
    if output is not None:
        rpe1_directory = output / "matrices" / "rpe1"
        rpe1_directory.mkdir(parents=True, exist_ok=True)
        for parameter, values in rpe1_references.items():
            values.to_csv(rpe1_directory / f"reference_{parameter}.csv")
        for branch, branch_values in organoid_references.items():
            directory = output / "matrices" / "organoid" / f"branch{branch}"
            directory.mkdir(parents=True, exist_ok=True)
            for parameter, values in branch_values.items():
                values.to_csv(directory / f"reference_{parameter}.csv")

    components = []
    gene_metrics = []
    for method in methods:
        method_sources = settings["methods"][method]
        if "rpe1" in method_sources:
            rows, metrics = evaluate_rpe1(
                method,
                method_sources["rpe1"],
                settings["datasets"],
                rpe1_references,
                artifact_root,
                output,
            )
            components.extend(rows)
            gene_metrics.extend(metrics)
        if "organoid" in method_sources:
            rows, metrics = evaluate_organoid(
                method,
                method_sources["organoid"],
                settings["datasets"],
                organoid_references,
                artifact_root,
                output,
            )
            components.extend(rows)
            gene_metrics.extend(metrics)

    component_frame = pd.DataFrame(components).sort_values(
        ["Parameter", "Method", "Dataset", "Branch"], na_position="first"
    )
    if args.check_only:
        print(component_frame.to_string(index=False))
        return

    metrics_frame = pd.concat(gene_metrics, ignore_index=True)
    systems, methods_summary = summarize(component_frame)
    component_frame.to_csv(output / "component_means.csv", index=False)
    metrics_frame.to_csv(output / "gene_metrics.csv", index=False)
    systems.to_csv(output / "system_means.csv", index=False)
    methods_summary.to_csv(output / "method_summary.csv", index=False)
    print(component_frame.to_string(index=False))
    print(f"Wrote kinetic-parameter metrics to {output}")


if __name__ == "__main__":
    main()
