#!/usr/bin/env python3
"""Run the four-method signed-transition ensemble on selected datasets."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import anndata as ad

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from ensemble.core import (
        METHOD_NAMES,
        VARIANT_NAMES,
        VARIANT_SPECS,
        aggregate_offdiagonal,
        align_method,
        align_method_values,
        complete_signed_operator,
        compute_signed_transition,
        DEGENERATE_ROW_TOLERANCE,
        compute_velocity_graph,
        install_fixed_neighbors,
        load_config,
        load_fixed_neighbor_asset,
        offdiagonal,
        project_with_scvelo,
        protocol_metadata,
        save_stream_plot,
        sha256,
        write_h5ad,
    )
    from ensemble.metrics import score_global_metrics, write_rankings
else:
    from .core import (
        METHOD_NAMES,
        VARIANT_NAMES,
        VARIANT_SPECS,
        aggregate_offdiagonal,
        align_method,
        align_method_values,
        complete_signed_operator,
        compute_signed_transition,
        DEGENERATE_ROW_TOLERANCE,
        compute_velocity_graph,
        install_fixed_neighbors,
        load_config,
        load_fixed_neighbor_asset,
        offdiagonal,
        project_with_scvelo,
        protocol_metadata,
        save_stream_plot,
        sha256,
        write_h5ad,
    )
    from .metrics import score_global_metrics, write_rankings


MODULE_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = MODULE_DIR.parent
DEFAULT_CONFIG = MODULE_DIR / "config.json"
DEFAULT_REFERENCE_NEIGHBORS = (
    REPOSITORY_ROOT / "evaluation" / "cell" / "reference_neighbors"
)
DEFAULT_CELL_METRICS = REPOSITORY_ROOT / "evaluation" / "cell" / "results" / "cell_metrics.csv"


def resolve_artifact(artifact_root: Path, locator: str) -> Path:
    """Resolve a config locator below --artifact-root without escaping it."""
    relative = Path(locator)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"artifact locator must be relative: {locator}")
    root = Path(artifact_root).resolve()
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise ValueError(
            f"artifact locator escapes artifact root: {locator}"
        ) from error
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    return candidate


def resolve_neighbor_asset(root: Path, filename: str) -> Path:
    """Resolve one repository-local NPZ named by the dataset config."""
    relative = Path(filename)
    if relative.is_absolute() or ".." in relative.parts or relative.name != filename:
        raise ValueError(f"invalid reference-neighbor filename: {filename}")
    path = Path(root) / relative
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def logical_input(locator: str, path: Path, layer: str) -> dict:
    return {
        "locator": f"velocity:{locator}",
        "sha256": sha256(path),
        "layer": layer,
    }


def write_scores(path: Path, rows: list[dict]) -> None:
    fields = ["dataset", "variant", "label", "cbdir", "batc"]
    output_rows = [
        {
            "dataset": row["dataset"],
            "variant": row["variant"],
            "label": row.get("label", row.get("variant_label")),
            "cbdir": row["cbdir"],
            "batc": row["batc"],
        }
        for row in rows
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)


def run_dataset(
    dataset: str,
    dataset_config: dict,
    *,
    artifact_root: Path,
    reference_neighbors_root: Path,
    n_jobs: int,
    variants: tuple[str, ...],
    output_root: Path,
    h5ad_root: Path | None,
    save_streams: bool,
) -> tuple[list[dict], dict]:
    reference_locator = dataset_config["reference"]
    reference_path = resolve_artifact(artifact_root, reference_locator)
    reference = ad.read_h5ad(reference_path)

    asset_path = resolve_neighbor_asset(
        reference_neighbors_root, dataset_config["reference_neighbor_file"]
    )
    asset = load_fixed_neighbor_asset(asset_path)
    reference_neighbor_detail = install_fixed_neighbors(reference, asset)

    logical_inputs = {
        "reference": logical_input(reference_locator, reference_path, "Ms"),
        "methods": {},
    }
    transitions = []
    method_metadata = {}
    for method in METHOD_NAMES:
        locator = dataset_config["methods"][method]
        path = resolve_artifact(artifact_root, locator)
        logical_inputs["methods"][method] = logical_input(
            locator, path, "velocity"
        )
        method_adata = ad.read_h5ad(path)
        aligned = align_method(reference, method_adata, method)
        positive, negative = compute_velocity_graph(aligned, n_jobs=n_jobs)
        transition, fallback_rows = compute_signed_transition(
            aligned, positive, negative, return_fallback=True
        )
        transitions.append(offdiagonal(transition))
        method_metadata[method] = {
            "n_fallback_rows": int(fallback_rows.sum()),
            "fallback_tolerance": DEGENERATE_ROW_TOLERANCE,
        }

    template, values = align_method_values(transitions)
    rows = []
    variant_metadata = {}
    for variant in variants:
        aggregate, aggregate_detail = aggregate_offdiagonal(template, values, variant)
        operator, completion_detail = complete_signed_operator(template, aggregate)
        result, invalid = project_with_scvelo(
            reference,
            operator,
            variant=variant,
            logical_inputs=logical_inputs,
        )
        metrics = score_global_metrics(
            result,
            dataset,
            dataset_config,
            operator=operator,
        )
        spec = VARIANT_SPECS[variant]
        output_locator = f"ensemble:{spec['directory']}/{dataset}"
        if h5ad_root is not None:
            h5ad_path = Path(h5ad_root) / spec["directory"] / f"{dataset}.h5ad"
            write_h5ad(result, h5ad_path)
        if save_streams:
            stream_path = (
                output_root
                / "figures"
                / "streams"
                / spec["directory"]
                / f"{dataset}.png"
            )
            save_stream_plot(
                result,
                stream_path,
                label_key=dataset_config["label"],
                title=f"{spec['label']} - {dataset}",
            )
        completion = dict(completion_detail)
        completion["invalid_embedded_values"] = invalid
        variant_metadata[variant] = {
            "label": spec["label"],
            "definition": spec["definition"],
            "aggregate": aggregate_detail,
            "completion": completion,
            "output_locator": output_locator,
            "metrics": {
                "cbdir": metrics["cbdir"],
                "batc": metrics["batc"],
            },
        }
        rows.append(
            {
                "dataset": dataset,
                "variant": variant,
                "variant_label": spec["label"],
                "cbdir": metrics["cbdir"],
                "batc": metrics["batc"],
                "n_obs": reference.n_obs,
                "n_common_edges": completion["n_edges"],
                "directional_mass_mean": completion["directional_mass_mean"],
                "neutral_mass_mean": completion["neutral_mass_mean"],
                "invalid_embedded_values": invalid,
            }
        )
    return rows, {
        "dataset": dataset,
        "reference_neighbor": reference_neighbor_detail,
        "inputs": logical_inputs,
        "methods": method_metadata,
        "variants": variant_metadata,
    }


def write_summary(path: Path, rows: list[dict], datasets: list[str]) -> None:
    lines = [
        "# Signed-T ensemble run",
        "",
        f"Datasets: {', '.join(datasets)}",
        "",
        "| Dataset | Variant | CBDir | BATC |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['dataset']} | {row.get('label', row.get('variant_label'))} | "
            f"{row['cbdir']:.8g} | {row['batc']:.8g} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--reference-neighbor-root",
        type=Path,
        default=DEFAULT_REFERENCE_NEIGHBORS,
    )
    parser.add_argument(
        "--dataset",
        action="append",
        help="dataset key; repeat to select more than one (default: all configured)",
    )
    parser.add_argument(
        "--variant",
        action="append",
        choices=VARIANT_NAMES,
        help="variant key; repeat to select more than one (default: both)",
    )
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument(
        "--h5ad-root",
        type=Path,
        help="optional external staging root for H5AD outputs",
    )
    parser.add_argument("--save-streams", action="store_true")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    datasets = args.dataset or list(config["datasets"])
    unknown = sorted(set(datasets).difference(config["datasets"]))
    if unknown:
        raise ValueError(f"unknown dataset keys: {unknown}")
    variants = tuple(args.variant or VARIANT_NAMES)
    if args.n_jobs < 1:
        raise ValueError("--n-jobs must be positive")

    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    rows = []
    runs = []
    for dataset in datasets:
        dataset_rows, metadata = run_dataset(
            dataset,
            config["datasets"][dataset],
            artifact_root=args.artifact_root,
            reference_neighbors_root=args.reference_neighbor_root,
            n_jobs=args.n_jobs,
            variants=variants,
            output_root=output_root,
            h5ad_root=args.h5ad_root,
            save_streams=args.save_streams,
        )
        rows.extend(dataset_rows)
        runs.append(metadata)

    write_scores(output_root / "scores.csv", rows)
    (output_root / "metadata.json").write_text(
        json.dumps(
            {
                "schema_version": 6,
                "status": "selected-datasets",
                "datasets": datasets,
                "variants": list(variants),
                "n_jobs": args.n_jobs,
                "protocol": protocol_metadata(),
                "runs": runs,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    write_summary(output_root / "summary.md", rows, datasets)
    if set(datasets) == set(config["datasets"]) and DEFAULT_CELL_METRICS.is_file():
        write_rankings(output_root / "scores.csv", DEFAULT_CELL_METRICS, args.config, output_root / "rankings_vs_20.csv", output_root / "ranking_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
