#!/usr/bin/env python3
"""Rebuild and evaluate the Case 3 kinetic-parameter experiment.

The script consumes saved method result H5AD objects; it neither prepares data
nor trains a velocity model. Most methods expose parameter layers directly.
Dynamo and TSvelo are explicit exceptions: their historical alpha CSVs were
derived from method-specific H5AD layouts, so this script reproduces those
extractions before applying the shared RPE1 evaluation transform.
"""

from __future__ import annotations

import argparse
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MODULE_DIR / "data" / "case3.json"


def dense(matrix):
    """Return a dense NumPy array for a dense or sparse AnnData layer."""
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


def read_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def source_path(working_root: Path, relative_path: str) -> Path:
    path = working_root / relative_path
    if not path.is_file():
        raise FileNotFoundError(f"Required source does not exist: {path}")
    return path


def read_rpe1_reference(path: Path, parameter: str) -> pd.DataFrame:
    table = pd.read_csv(path, skiprows=1).set_index("gene_symbol")
    prefix = "norm_kappa_" if parameter == "alpha" else "norm_gamma_"
    values = table.loc[:, [name for name in table.columns if name.startswith(prefix)]]
    values.columns = range(values.shape[1])
    return values.apply(pd.to_numeric, errors="coerce")


def read_gene_subset(path: Path) -> pd.Index:
    """Read the legacy RPE1 alpha gene selection, not its derived values."""
    values = pd.read_csv(path, index_col=0)
    if values.index.has_duplicates:
        raise ValueError(f"{path}: the gene subset contains duplicate names")
    return pd.Index(values.index)


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
    """Calculate a per-gene score after aligning genes and timepoints."""
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
    """Check only the fields used by the direct-layer extraction path."""
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

    required_layers = list(source["layers"].values()) + source.get(
        "required_layers", []
    )
    missing_layers = [name for name in required_layers if name not in adata.layers]
    if missing_layers:
        raise KeyError(f"{dataset}: missing parameter layers {missing_layers}")


def layer_frame(adata: ad.AnnData, layer: str) -> pd.DataFrame:
    return pd.DataFrame(
        dense(adata.layers[layer]), index=adata.obs_names, columns=adata.var_names
    )


def rpe1_raw_from_layer(
    adata: ad.AnnData, layer: str, settings: dict
) -> pd.DataFrame:
    """Average a direct cell-by-gene layer at each fixed RPE1 position."""
    values = layer_frame(adata, layer)
    positions = pd.to_numeric(adata.obs[settings["obs_key"]], errors="raise")
    return values.groupby(positions).mean().reindex(range(settings["positions"])).T


def process_rpe1(raw: pd.DataFrame, settings: dict, clip_negative: bool) -> pd.DataFrame:
    """Apply the historical shared RPE1 smoothing, normalization, and log2 step."""
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
    """Average a direct organoid layer in the fixed 200 pseudotime bins."""
    values = layer_frame(adata, layer)
    bins = pd.cut(
        aligned_time, bins=settings["bins"], labels=False, include_lowest=True
    )
    return values.groupby(bins).mean().reindex(range(settings["bins"])).T


def process_organoid(raw: pd.DataFrame, settings: dict) -> pd.DataFrame:
    """Apply the historical organoid smoothing and early-bin normalization."""
    smoothed = raw.T.rolling(
        settings["smooth_window"], center=True, min_periods=1
    ).mean().T
    baseline = smoothed.iloc[:, : settings["baseline_bins"]].mean(axis=1)
    return smoothed.div(baseline + settings["epsilon"], axis=0)


def decode_h5_strings(values) -> pd.Index:
    """Decode an HDF5 string vector without loading the complete AnnData object."""
    decoded = [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in np.asarray(values)
    ]
    return pd.Index(decoded)


def read_csr_layer(
    path: Path, layer: str, position_key: str
) -> tuple[sparse.csr_matrix, pd.Index, pd.Series]:
    """Read Dynamo's sparse cell-wise alpha layer and its labels directly."""
    with h5py.File(path, "r") as handle:
        layer_group = handle["layers"][layer]
        if layer_group.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{path}: layers/{layer} is not stored as CSR")
        shape = tuple(int(value) for value in layer_group.attrs["shape"])
        matrix = sparse.csr_matrix(
            (
                np.asarray(layer_group["data"]),
                np.asarray(layer_group["indices"]),
                np.asarray(layer_group["indptr"]),
            ),
            shape=shape,
        )
        var_names = decode_h5_strings(handle["var"]["_index"])
        positions = pd.Series(
            np.asarray(handle["obs"][position_key]), dtype=float, name=position_key
        )
    if matrix.shape[1] != len(var_names):
        raise ValueError(f"{path}: layer and var index have incompatible dimensions")
    return matrix, var_names, positions


def dynamo_rpe1_raw(
    source: dict, reference: pd.DataFrame, settings: dict, working_root: Path
) -> pd.DataFrame:
    """Rebuild the historical Dynamo alpha matrix from its labelled H5AD.

    The historical CSV used the sparse cell-wise alpha layer, binned positions
    with pd.cut(range(0, 302)), retained reference genes, and dropped genes that
    were zero at every observed bin. This executable path replaces that CSV.
    """
    path = source_path(working_root, source["result"])
    matrix, var_names, positions = read_csr_layer(
        path, source["layer"], source["position_key"]
    )
    if matrix.shape[0] != settings["n_obs"]:
        raise ValueError(
            f"Dynamo RPE1: expected {settings['n_obs']} cells, found {matrix.shape[0]}"
        )
    if positions.shape[0] != matrix.shape[0]:
        raise ValueError("Dynamo RPE1: position vector does not match cell count")
    if not np.isfinite(positions).all():
        raise ValueError("Dynamo RPE1: Cell_cycle_possition contains non-finite values")

    # Preserve the source H5AD's gene order, matching the original CSV script.
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
    if source.get("drop_all_zero_genes", False):
        raw = raw.loc[~(raw == 0).all(axis=1)]
    return raw


def tsvelo_rpe1_raw(
    source: dict, reference: pd.DataFrame, settings: dict, working_root: Path
) -> pd.DataFrame:
    """Rebuild the historical TSvelo alpha matrix from its two lineage H5ADs."""
    alpha_frames = []
    positions = []
    for relative_path in source["results"]:
        path = source_path(working_root, relative_path)
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

    # Repeated cells belong to both lineages: average alpha and retain first label.
    combined = pd.concat(alpha_frames)
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


def write_csv(
    output_dir: Path,
    dataset: str,
    branch: str | None,
    filename: str,
    values: pd.DataFrame,
    kind: str,
) -> str:
    """Persist a plot-ready matrix, separating raw extraction from processing."""
    path = output_dir / kind / dataset
    if branch is not None:
        path = path / f"branch{branch}"
    path.mkdir(parents=True, exist_ok=True)
    file_path = path / filename
    values.to_csv(file_path)
    return str(file_path.relative_to(output_dir))


def write_rpe1_references(
    output_dir: Path, references: dict[str, pd.DataFrame]
) -> dict[str, str]:
    return {
        parameter: write_csv(
            output_dir,
            "rpe1",
            None,
            f"reference_{parameter}.csv",
            values,
            "matrices",
        )
        for parameter, values in references.items()
    }


def write_organoid_references(
    output_dir: Path, references: dict[str, dict[str, pd.DataFrame]]
) -> dict[str, dict[str, str]]:
    return {
        branch: {
            parameter: write_csv(
                output_dir,
                "organoid",
                branch,
                f"reference_{parameter}.csv",
                values,
                "matrices",
            )
            for parameter, values in branch_values.items()
        }
        for branch, branch_values in references.items()
    }


def record_matrices(
    output_dir: Path | None,
    dataset: str,
    branch: str | None,
    method: str,
    parameter: str,
    raw: pd.DataFrame,
    processed: pd.DataFrame,
) -> dict | None:
    if output_dir is None:
        return None
    return {
        "method": method,
        "dataset": dataset,
        "branch": branch,
        "parameter": parameter,
        "raw": write_csv(
            output_dir,
            dataset,
            branch,
            f"{method}_{parameter}_raw.csv",
            raw,
            "matrices",
        ),
        "processed": write_csv(
            output_dir,
            dataset,
            branch,
            f"{method}_{parameter}.csv",
            processed,
            "trajectories",
        ),
    }


def evaluate_rpe1(
    method: str,
    source: dict,
    settings: dict,
    references: dict[str, pd.DataFrame],
    working_root: Path,
    output_dir: Path | None,
) -> tuple[list[dict], list[pd.DataFrame], list[dict]]:
    """Extract RPE1 parameters from each method's actual result layout."""
    adapter = source["adapter"]
    raw_by_parameter: dict[str, pd.DataFrame]
    if adapter == "h5ad_layers":
        adata = ad.read_h5ad(source_path(working_root, source["result"]))
        validate_standard_result(adata, "rpe1", settings, source)
        raw_by_parameter = {
            parameter: rpe1_raw_from_layer(adata, layer, settings["rpe1"])
            for parameter, layer in source["layers"].items()
        }
    elif adapter == "dynamo_labelled_alpha":
        raw_by_parameter = {
            "alpha": dynamo_rpe1_raw(
                source, references["alpha"], settings["rpe1"], working_root
            )
        }
    elif adapter == "tsvelo_lineage_alpha":
        raw_by_parameter = {
            "alpha": tsvelo_rpe1_raw(
                source, references["alpha"], settings["rpe1"], working_root
            )
        }
    else:
        raise ValueError(f"{method} RPE1: unsupported adapter {adapter!r}")

    components = []
    gene_metrics = []
    matrices = []
    for parameter, raw in raw_by_parameter.items():
        processed = process_rpe1(
            raw, settings["rpe1"], source.get("clip_negative", False)
        )
        component, metrics = score_component(
            method, "rpe1", None, parameter, processed, references[parameter]
        )
        components.append(component)
        gene_metrics.append(metrics)
        record = record_matrices(
            output_dir, "rpe1", None, method, parameter, raw, processed
        )
        if record is not None:
            matrices.append(record)
    return components, gene_metrics, matrices


def evaluate_organoid(
    method: str,
    source: dict,
    settings: dict,
    references: dict[str, dict[str, pd.DataFrame]],
    working_root: Path,
    output_dir: Path | None,
) -> tuple[list[dict], list[pd.DataFrame], list[dict]]:
    if source["adapter"] != "h5ad_layers":
        raise ValueError(f"{method} organoid: unsupported adapter {source['adapter']!r}")
    adata = ad.read_h5ad(source_path(working_root, source["result"]))
    validate_standard_result(adata, "organoid", settings, source)

    components = []
    gene_metrics = []
    matrices = []
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
            record = record_matrices(
                output_dir,
                "organoid",
                branch,
                method,
                parameter,
                raw,
                processed,
            )
            if record is not None:
                matrices.append(record)
    return components, gene_metrics, matrices


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--working-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--methods", nargs="*")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.check_only == (args.output is not None):
        parser.error("use exactly one of --output or --check-only")

    settings = read_config(args.config)
    methods = args.methods or list(settings["methods"])
    unknown = [method for method in methods if method not in settings["methods"]]
    if unknown:
        parser.error(f"unknown methods: {', '.join(unknown)}")

    output = None if args.check_only else args.output.resolve()
    if output is not None:
        ensure_output(output)
    working_root = args.working_root.resolve()

    rpe1_reference_path = source_path(working_root, settings["rpe1"]["reference"])
    full_rpe1_alpha = read_rpe1_reference(rpe1_reference_path, "alpha")
    alpha_gene_subset = read_gene_subset(
        source_path(working_root, settings["rpe1"]["alpha_gene_subset"])
    )
    missing_alpha_genes = alpha_gene_subset.difference(full_rpe1_alpha.index)
    if not missing_alpha_genes.empty:
        raise ValueError(
            "RPE1 alpha gene subset is not contained in Table S1: "
            f"{', '.join(missing_alpha_genes[:5])}"
        )
    rpe1_references = {
        "alpha": full_rpe1_alpha.loc[alpha_gene_subset],
        "gamma": read_rpe1_reference(rpe1_reference_path, "gamma"),
    }
    organoid_reference_path = source_path(
        working_root, settings["organoid"]["reference"]
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

    reference_files = {}
    if output is not None:
        reference_files = {
            "rpe1": write_rpe1_references(output, rpe1_references),
            "organoid": write_organoid_references(output, organoid_references),
        }

    components = []
    gene_metrics = []
    matrices = []
    for method in methods:
        method_sources = settings["methods"][method]
        if "rpe1" in method_sources:
            rows, metrics, files = evaluate_rpe1(
                method,
                method_sources["rpe1"],
                settings,
                rpe1_references,
                working_root,
                output,
            )
            components.extend(rows)
            gene_metrics.extend(metrics)
            matrices.extend(files)
        if "organoid" in method_sources:
            rows, metrics, files = evaluate_organoid(
                method,
                method_sources["organoid"],
                settings,
                organoid_references,
                working_root,
                output,
            )
            components.extend(rows)
            gene_metrics.extend(metrics)
            matrices.extend(files)

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
    (output / "manifest.json").write_text(
        json.dumps(
            {
                "created_at": datetime.now(timezone.utc).isoformat(),
                "config": str(args.config.resolve()),
                "methods": methods,
                "method_sources": {
                    method: settings["methods"][method] for method in methods
                },
                "reference_matrices": reference_files,
                "parameter_matrices": matrices,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(component_frame.to_string(index=False))
    print(f"Wrote Case 3 metrics to {output}")


if __name__ == "__main__":
    main()
