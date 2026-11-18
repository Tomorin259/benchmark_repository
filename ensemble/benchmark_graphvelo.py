#!/usr/bin/env python3
"""Compare native-operator ensembles with GraphVelo and optional TSvelo."""

from __future__ import annotations

import argparse
import csv
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys

import anndata as ad
from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ensemble.benchmark import resolve_artifact, write_scores
from ensemble.core import (
    DEFAULT_CONFIG, METHOD_NAMES, VARIANT_NAMES, aggregate_offdiagonal,
    align_method, align_method_values, complete_signed_operator,
    compute_signed_transition, compute_velocity_graph, dense_float64, install_fixed_neighbors,
    load_config, load_fixed_neighbor_asset, offdiagonal, protocol_metadata,
    row_l1_norm, sha256,
)
from ensemble.geometric_median import row_geometric_median
from ensemble.metrics import score_global_metrics
from ensemble.reproject import project_graphvelo
from evaluation.cell.reference_embeddings import load_forebrain_frames


ROOT = Path(__file__).resolve().parents[1]
METHODS = (*METHOD_NAMES, "GraphVelo_scv")
ZERO_TOLERANCE = 1e-12
TSVELO_CONTRASTS = (
    ("add_tsvelo_to_four", "4", "5_tsvelo"),
    ("add_tsvelo_to_graphvelo", "5", "6"),
    ("add_graphvelo_to_tsvelo", "5_tsvelo", "6"),
)


def member_selections(replace_each=False, include_tsvelo=False):
    """Define the original, added-member, and four replacement configurations."""
    selections = {"four": tuple(range(4)), "five": tuple(range(5))}
    if include_tsvelo:
        selections["tsvelo_five"] = (0, 1, 2, 3, 5)
        selections["tsvelo_graphvelo_six"] = tuple(range(6))
    if replace_each:
        for removed in METHOD_NAMES:
            selections["replace_" + removed] = tuple(
                index for index, method in enumerate(METHODS) if method != removed
            )
    return selections


def load_graphvelo(row, artifact_root, phi_root, reference):
    """Read the current native Phi, verify its output, and normalize signed mass."""
    path = resolve_artifact(artifact_root, row["result_path"].split(":", 1)[1])
    if sha256(path) != row["sha256"]:
        raise ValueError("GraphVelo source differs from the registered input")
    with h5py.File(path, "r") as handle:
        names = read_elem(handle["obs"]).index
        coordinates = read_elem(handle["obsm/X_umap"])
        expected = read_elem(handle["obsm/gv_umap"])
        if "gv" in handle["obsp"]:
            phi = sparse.csr_matrix(read_elem(handle["obsp/gv"]))
            provenance = {"kind": "stored_obsp_gv"}
        else:
            folder = phi_root / row["evaluation_dataset"]
            metadata_path = folder / "metadata.json"
            metadata = json.loads(metadata_path.read_text())
            if (metadata["source_sha256"] != row["sha256"]
                    or not metadata["summary"]["all_outputs_entrywise_equal"]):
                raise ValueError("Phi sidecar does not verify the current GraphVelo source")
            for name in ("phi.npz", "cell_ids.npy"):
                if sha256(folder / name) != metadata["outputs_sha256"][name]:
                    raise ValueError("Phi sidecar hash differs")
            phi = sparse.load_npz(folder / "phi.npz")
            sidecar_names = np.load(folder / "cell_ids.npy", allow_pickle=False)
            if not np.array_equal(sidecar_names, np.asarray(names, dtype=str)):
                raise ValueError("Phi sidecar cell order differs from the method H5AD")
            provenance = {
                "kind": "verified_sidecar",
                "phi_locator": row["evaluation_dataset"] + "/phi.npz",
                "phi_sha256": sha256(folder / "phi.npz"),
                "metadata_sha256": sha256(metadata_path),
            }
    if phi.shape != (len(names), len(names)) or not np.isfinite(phi.data).all():
        raise ValueError("native Phi has invalid axes or values")
    replay = np.zeros(coordinates.shape, dtype=np.float64)
    for cell in range(phi.shape[0]):
        start, end = phi.indptr[cell:cell + 2]
        neighbors = phi.indices[start:end]
        replay[cell] = phi.data[start:end] @ (coordinates[neighbors] - coordinates[cell])
    replay_error = float(np.max(np.abs(replay - expected)))
    if not np.allclose(replay, expected, rtol=1e-12, atol=1e-11):
        raise ValueError("native Phi no longer reproduces its stored embedding")
    positions = pd.Index(names).get_indexer(reference.obs_names)
    if (len(names) != reference.n_obs or np.any(positions < 0)
            or not pd.Index(names).is_unique):
        raise ValueError("GraphVelo and ensemble reference cell sets differ")
    phi = offdiagonal(phi[positions][:, positions])
    mass = row_l1_norm(phi)
    scale = np.divide(1.0, mass, out=np.zeros_like(mass), where=mass > ZERO_TOLERANCE)
    normalized = phi.multiply(scale[:, None]).tocsr()
    normalized.eliminate_zeros()
    return normalized, {
        **provenance, "source_locator": row["result_path"], "source_sha256": row["sha256"],
        "native_embedding_replay_max_abs_error": replay_error,
        "zero_rows": int(np.count_nonzero(mass <= ZERO_TOLERANCE)),
        "native_mass_median": float(np.median(mass)),
        "native_mass_mean": float(mass.mean()),
        "normalized_mass_max": float(row_l1_norm(normalized).max()),
    }


def run_dataset(dataset, spec, graphvelo_row, args, release, tsvelo_row=None):
    reference_path = resolve_artifact(args.artifact_root, spec["reference"])
    expected = release["inputs"][f"input:{dataset}:reference"]
    if (reference_path.stat().st_size != expected["size_bytes"]
            or sha256(reference_path) != expected["sha256"]):
        raise ValueError("reference differs from the accepted four-member input")
    reference = ad.read_h5ad(reference_path)
    asset_path = ROOT / "evaluation/cell/reference_neighbors" / spec["reference_neighbor_file"]
    asset = load_fixed_neighbor_asset(asset_path)
    install_fixed_neighbors(reference, asset)
    if dataset == "forebrain":
        frames = load_forebrain_frames()
        order = pd.Index(frames["obs_names"]).get_indexer(reference.obs_names)
        if np.any(order < 0):
            raise ValueError("Forebrain cells differ from the frozen coordinate assets")
        for basis in ("pca", "umap"):
            reference.obsm["X_" + basis] = frames["X_" + basis][order].copy()

    transitions, sources, fallback_counts = [], {}, {}
    run_methods = (*METHOD_NAMES, "TSvelo") if args.include_tsvelo else METHOD_NAMES
    methods = (*METHODS, "TSvelo") if args.include_tsvelo else METHODS
    for method in run_methods:
        if method == "TSvelo":
            locator = tsvelo_row["result_path"].split(":", 1)[1]
            expected = {"sha256": tsvelo_row["sha256"],
                        "size_bytes": int(tsvelo_row["file_size_bytes"])}
        else:
            locator = spec["methods"][method]
            expected = release["inputs"][f"input:{dataset}:method:{method}"]
        path = resolve_artifact(args.artifact_root, locator)
        actual = sha256(path)
        if path.stat().st_size != expected["size_bytes"] or actual != expected["sha256"]:
            raise ValueError("method input differs from the accepted ensemble: " + method)
        print(dataset + ": constructing " + method + " signed operator", flush=True)
        native = ad.read_h5ad(path)
        aligned = align_method(reference, native, method)
        positive, negative = compute_velocity_graph(aligned, n_jobs=args.n_jobs)
        transition, fallback = compute_signed_transition(
            aligned, positive, negative, return_fallback=True
        )
        transitions.append(offdiagonal(transition))
        sources[method] = {"locator": locator, "sha256": actual}
        if method == "TSvelo":
            sources[method]["native_gene_count"] = native.n_vars
            sources[method]["selected_velocity_genes"] = int(aligned.var["velocity_genes"].sum())
            sources[method]["nonfinite_native_velocity_entries"] = int(
                np.count_nonzero(~np.isfinite(dense_float64(native.layers["velocity"])))
            )
        fallback_counts[method] = int(fallback.sum())
        del native, aligned, positive, negative, transition

    graphvelo, phi_metadata = load_graphvelo(
        graphvelo_row, args.artifact_root, args.phi_root, reference
    )
    transitions.insert(4, graphvelo)
    template, values = align_method_values(transitions, n_methods=len(methods))
    folder = args.output_root / dataset
    folder.mkdir()
    if args.include_tsvelo:
        for method, transition in zip(methods, transitions):
            sparse.save_npz(folder / ("operator_input_" + method + ".npz"), transition)
    if not args.replace_each and not args.include_tsvelo:
        sparse.save_npz(folder / "graphvelo_signed.npz", graphvelo)
    bases = ("pca", "umap") if dataset == "forebrain" else ("umap",)
    plot = ad.AnnData(obs=reference.obs[[spec["label"]]].copy())
    plot.uns["neighbors"] = {"indices": reference.uns["neighbors"]["indices"].copy()}
    arrays = {
        "cell_ids": np.asarray(reference.obs_names, dtype=str),
        "labels": np.asarray(reference.obs[spec["label"]].astype(str)),
        "neighbors": plot.uns["neighbors"]["indices"],
    }
    for basis in bases:
        plot.obsm["X_" + basis] = reference.obsm["X_" + basis]
        arrays["X_" + basis] = reference.obsm["X_" + basis]
    rows, diagnostics = [], {}
    baseline_table = {
        r["variant"]: r for r in csv.DictReader(
            (ROOT / "ensemble/results/graphvelo_projection/scores.csv").open()
        ) if r["dataset"] == dataset
    }
    for selection, member_ids in member_selections(
        args.replace_each, args.include_tsvelo
    ).items():
        count = len(member_ids)
        selected = values[list(member_ids)]
        for aggregation in ("mean", "geomedian"):
            variant = f"{aggregation}{count}"
            if selection.startswith("replace_"):
                variant += "_" + selection
            elif selection == "tsvelo_five":
                variant += "_tsvelo"
            legacy = VARIANT_NAMES[0 if aggregation == "mean" else 1]
            if count == 4:
                aggregate, detail = aggregate_offdiagonal(template, selected, legacy)
            elif aggregation == "mean":
                aggregate = selected.mean(axis=0)
                detail = {"n_methods": count, "rule": "equal arithmetic mean"}
            else:
                aggregate, detail, cell_detail = row_geometric_median(
                    template, selected, return_cell_diagnostics=True
                )
                if not args.replace_each and not (args.include_tsvelo and selection == "five"):
                    np.savez_compressed(
                        folder / (variant + "_weights.npz"),
                        cell_ids=arrays["cell_ids"],
                        methods=np.asarray([methods[index] for index in member_ids]),
                        coefficients=cell_detail["coefficients"],
                    )
            operator, completion = complete_signed_operator(template, aggregate)
            matrix_error = 0.0
            if selection == "four":
                entry = release["outputs"][f"output:{legacy}:{dataset}:h5ad"]
                frozen = resolve_artifact(args.frozen_ensemble_root, entry["relative_path"])
                if frozen.stat().st_size != entry["size_bytes"] or sha256(frozen) != entry["sha256"]:
                    raise ValueError("frozen baseline H5AD hash differs")
                with h5py.File(frozen, "r") as handle:
                    previous = read_elem(handle["obsp/ensemble_transition_signed"])
                difference = operator - previous
                difference.eliminate_zeros()
                matrix_error = float(np.max(np.abs(difference.data))) if difference.nnz else 0.0
                if matrix_error > 1e-9:
                    raise ValueError(f"four-member operator changed: {dataset}/{variant}/{matrix_error}")
            elif selection == "five" and (args.replace_each or args.include_tsvelo):
                baseline_folder = ROOT / "ensemble/results/graphvelo_five" / dataset
                baseline_metadata = json.loads((baseline_folder / "metadata.json").read_text())
                name = f"{variant}_operator.npz"
                if sha256(baseline_folder / name) != baseline_metadata["outputs_sha256"][name]:
                    raise ValueError("five-member control operator changed")
                previous = sparse.load_npz(baseline_folder / name)
                difference = operator - previous
                difference.eliminate_zeros()
                matrix_error = float(np.max(np.abs(difference.data))) if difference.nnz else 0.0
                if matrix_error > 1e-9:
                    raise ValueError("five-member control no longer reproduces")
            elif not args.replace_each:
                sparse.save_npz(folder / f"{variant}_operator.npz", operator)
            for basis in bases:
                field = project_graphvelo(operator, plot.obsm["X_" + basis])
                plot.obsm["ensemble_velocity_" + basis] = field
                arrays[variant + "_" + basis] = field
            metrics = score_global_metrics(plot, dataset, spec)
            if selection == "four" and any(
                abs(metrics[m] - float(baseline_table[legacy][m])) > 5e-7
                for m in ("cbdir", "batc")
            ):
                raise ValueError("four-member projected scores changed")
            rows.append({
                "dataset": dataset, "variant": variant,
                "label": (("Mean-" if aggregation == "mean" else "GeoMedian-") + str(count)
                          + (" | " + selection if selection.startswith("replace_")
                             else " | TSvelo" if selection == "tsvelo_five"
                             else " | TSvelo + GraphVelo" if count == 6 else "")),
                "selection": selection,
                "members": [methods[index] for index in member_ids],
                "cbdir": metrics["cbdir"], "batc": metrics["batc"],
            })
            diagnostics[variant] = {
                "aggregation": detail, "completion": completion,
                "metrics": metrics, "frozen_operator_max_abs_error": matrix_error,
                "control_operator_checked": selection == "four" or (
                    selection == "five" and (args.replace_each or args.include_tsvelo)
                ),
            }
            print(f"{dataset}/{variant}: CBDir {metrics['cbdir']:.6f}; BATC {metrics['batc']:.6f}",
                  flush=True)
    np.savez_compressed(folder / "projections.npz", **arrays)
    output_hashes = {p.name: sha256(p) for p in folder.glob("*.npz")}
    return {
        "dataset": dataset, "n_cells": reference.n_obs, "rows": rows,
        "reference": {"locator": spec["reference"], "sha256": sha256(reference_path)},
        "reference_neighbor_sha256": asset["sha256"],
        "methods": sources, "fallback_counts": fallback_counts,
        "graphvelo": phi_metadata, "diagnostics": diagnostics,
        "outputs_sha256": output_hashes,
    }


def write_results(output, runs, protocol):
    rows = [row for run in runs for row in run["rows"]]
    write_scores(output / "scores.csv", rows)
    if protocol.get("replacement_ablation", False):
        comparison, summaries = [], []
        for run in runs:
            by_variant = {row["variant"]: row for row in run["rows"]}
            for removed in METHOD_NAMES:
                for aggregation in ("mean", "geomedian"):
                    four = by_variant[aggregation + "4"]
                    five = by_variant[aggregation + "5"]
                    candidate = by_variant[aggregation + "4_replace_" + removed]
                    record = {"dataset": run["dataset"], "removed_method": removed,
                              "aggregation": aggregation}
                    for metric in ("cbdir", "batc"):
                        record.update({
                            "four_" + metric: four[metric],
                            "replacement_" + metric: candidate[metric],
                            "five_" + metric: five[metric],
                            "delta_" + metric: candidate[metric] - four[metric],
                            "delta_vs_five_" + metric: candidate[metric] - five[metric],
                        })
                    comparison.append(record)
        for removed in METHOD_NAMES:
            for aggregation in ("mean", "geomedian"):
                selected = [row for row in comparison
                            if row["removed_method"] == removed and row["aggregation"] == aggregation]
                for metric in ("cbdir", "batc"):
                    delta = np.asarray([row["delta_" + metric] for row in selected])
                    delta_five = np.asarray([row["delta_vs_five_" + metric] for row in selected])
                    summaries.append({
                        "removed_method": removed, "aggregation": aggregation,
                        "metric": metric, "n_datasets": len(selected),
                        "four_mean": float(np.mean([row["four_" + metric] for row in selected])),
                        "replacement_mean": float(np.mean([row["replacement_" + metric] for row in selected])),
                        "five_mean": float(np.mean([row["five_" + metric] for row in selected])),
                        "delta_mean": float(delta.mean()),
                        "delta_vs_five_mean": float(delta_five.mean()),
                        "improved": int(np.count_nonzero(delta > 1e-12)),
                        "reduced": int(np.count_nonzero(delta < -1e-12)),
                        "improved_vs_five": int(np.count_nonzero(delta_five > 1e-12)),
                        "reduced_vs_five": int(np.count_nonzero(delta_five < -1e-12)),
                    })
    elif protocol.get("tsvelo_included", False):
        comparison, summaries = [], []
        for run in runs:
            by_variant = {row["variant"]: row for row in run["rows"]}
            for contrast, baseline_suffix, candidate_suffix in TSVELO_CONTRASTS:
                for aggregation in ("mean", "geomedian"):
                    baseline = by_variant[aggregation + baseline_suffix]
                    candidate = by_variant[aggregation + candidate_suffix]
                    record = {"dataset": run["dataset"], "contrast": contrast,
                              "aggregation": aggregation}
                    for metric in ("cbdir", "batc"):
                        record.update({
                            "baseline_" + metric: baseline[metric],
                            "candidate_" + metric: candidate[metric],
                            "delta_" + metric: candidate[metric] - baseline[metric],
                        })
                    comparison.append(record)
        for contrast, _, _ in TSVELO_CONTRASTS:
            for aggregation in ("mean", "geomedian"):
                selected = [row for row in comparison
                            if row["contrast"] == contrast and row["aggregation"] == aggregation]
                for metric in ("cbdir", "batc"):
                    delta = np.asarray([row["delta_" + metric] for row in selected])
                    summaries.append({
                        "contrast": contrast, "aggregation": aggregation, "metric": metric,
                        "n_datasets": len(selected),
                        "baseline_mean": float(np.mean([row["baseline_" + metric] for row in selected])),
                        "candidate_mean": float(np.mean([row["candidate_" + metric] for row in selected])),
                        "delta_mean": float(delta.mean()), "median_delta": float(np.median(delta)),
                        "improved": int(np.count_nonzero(delta > 1e-12)),
                        "reduced": int(np.count_nonzero(delta < -1e-12)),
                    })
    else:
        comparison = []
        for run in runs:
            by_variant = {r["variant"]: r for r in run["rows"]}
            for aggregation in ("mean", "geomedian"):
                old, new = by_variant[aggregation + "4"], by_variant[aggregation + "5"]
                comparison.append({
                    "dataset": run["dataset"], "aggregation": aggregation,
                    "four_cbdir": old["cbdir"], "five_cbdir": new["cbdir"],
                    "delta_cbdir": new["cbdir"] - old["cbdir"],
                    "four_batc": old["batc"], "five_batc": new["batc"],
                    "delta_batc": new["batc"] - old["batc"],
                })
        summaries = []
        for aggregation in ("mean", "geomedian"):
            selected = [r for r in comparison if r["aggregation"] == aggregation]
            for metric in ("cbdir", "batc"):
                delta = np.asarray([r["delta_" + metric] for r in selected])
                summaries.append({
                    "aggregation": aggregation, "metric": metric, "n_datasets": len(selected),
                    "four_mean": float(np.mean([r["four_" + metric] for r in selected])),
                    "five_mean": float(np.mean([r["five_" + metric] for r in selected])),
                    "delta_mean": float(delta.mean()),
                    "improved": int(np.count_nonzero(delta > 1e-12)),
                    "reduced": int(np.count_nonzero(delta < -1e-12)),
                })
    for name, records in (("comparison.csv", comparison), ("summary.csv", summaries)):
        with (output / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
    (output / "metadata.json").write_text(json.dumps({
        "protocol": protocol, "completed_datasets": [run["dataset"] for run in runs],
        "dataset_metadata_sha256": {
            run["dataset"]: sha256(output / run["dataset"] / "metadata.json") for run in runs
        },
        "table_sha256": {
            name: sha256(output / name) for name in ("scores.csv", "comparison.csv", "summary.csv")
        },
    }, indent=2, sort_keys=True) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--phi-root", type=Path, required=True)
    parser.add_argument("--frozen-ensemble-root", type=Path, required=True)
    parser.add_argument("--release-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--dataset", action="append")
    parser.add_argument("--n-jobs", type=int, default=8)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--replace-each", action="store_true",
        help="also replace each original method with GraphVelo in four-member ensembles",
    )
    mode.add_argument(
        "--include-tsvelo", action="store_true",
        help="also compare original four plus TSvelo, with and without native GraphVelo",
    )
    args = parser.parse_args(argv)
    config = load_config(DEFAULT_CONFIG)
    selected = args.dataset or list(config["datasets"])
    if len(selected) != len(set(selected)) or set(selected).difference(config["datasets"]):
        raise ValueError("select unique configured dataset keys")
    if args.n_jobs < 1:
        raise ValueError("--n-jobs must be positive")
    release = json.loads(args.release_manifest.read_text())
    input_path = ROOT / "evaluation/cell/inputs.tsv"
    with input_path.open() as handle:
        registered = list(csv.DictReader(handle, delimiter="\t"))
    graphvelo = {r["evaluation_dataset"]: r for r in registered
                 if r["operational_id"] == "graphvelo_scv"}
    tsvelo = {r["evaluation_dataset"]: r for r in registered
              if r["operational_id"] == "TSvelo"}
    methods = (*METHODS, "TSvelo") if args.include_tsvelo else METHODS
    protocol = {
        "status": "exploratory", "methods": list(methods),
        "replacement_ablation": args.replace_each, "tsvelo_included": args.include_tsvelo,
        "member_selections": {name: [methods[index] for index in indices]
                              for name, indices in member_selections(
                                  args.replace_each, args.include_tsvelo
                              ).items()},
        "repository_head": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True
        ).strip(),
        "code_sha256": {
            name: sha256(ROOT / name) for name in
            ("ensemble/benchmark_graphvelo.py", "ensemble/core.py",
             "ensemble/geometric_median.py", "ensemble/reproject.py", "ensemble/metrics.py")
        },
        "config_sha256": sha256(DEFAULT_CONFIG), "inputs_sha256": sha256(input_path),
        "release_manifest_sha256": sha256(args.release_manifest),
        "four_method_protocol": protocol_metadata(), "n_jobs": args.n_jobs,
        "graphvelo_normalization": "offdiagonal Phi / row signed-L1; near-zero rows give zero",
        "four_method_normalization": "retain each native signed-T offdiagonal directional mass",
        "projection": "sum_j W_ij (Y_j-Y_i); fixed frozen coordinates; no training",
        "versions": {name: version(name) for name in
                     ("numpy", "scipy", "scvelo", "anndata", "scanpy", "scikit-learn")},
    }
    args.output_root.mkdir(parents=True, exist_ok=True)
    runs = []
    for dataset in selected:
        folder = args.output_root / dataset
        if folder.exists():
            run = json.loads((folder / "metadata.json").read_text())
            if run["protocol"] != protocol:
                raise ValueError("existing dataset uses a different protocol: " + dataset)
            for name, expected in run["outputs_sha256"].items():
                if sha256(folder / name) != expected:
                    raise ValueError("existing output hash differs: " + dataset)
            print(dataset + ": resumed verified result", flush=True)
        else:
            spec = config["datasets"][dataset]
            key = Path(spec["reference_neighbor_file"]).stem
            row = graphvelo[key]
            tsvelo_row = tsvelo[key] if args.include_tsvelo else None
            run = run_dataset(dataset, spec, row, args, release, tsvelo_row=tsvelo_row)
            run["protocol"] = protocol
            (folder / "metadata.json").write_text(json.dumps(run, indent=2, sort_keys=True) + "\n")
        runs.append(run)
        write_results(args.output_root, runs, protocol)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
