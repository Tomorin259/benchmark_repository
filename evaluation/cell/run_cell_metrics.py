#!/usr/bin/env python3
"""Run the fixed-reference cell evaluation one method/dataset pair at a time.

The controller deliberately invokes two environments.  scvelo 0.2.5 prepares
the fixed-graph AnnData and calculates the four non-TiVelo metrics; TiVelo
0.1.4 then calculates the two official TiVelo metrics from the same saved
AnnData.  A copied TiVelo formula runs in the first environment only as a
per-pair parity check.  Source result files are always read-only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = MODULE_DIR.parents[1]
DEFAULT_CONFIG = MODULE_DIR / "config.json"
DEFAULT_INPUTS = MODULE_DIR / "inputs.tsv"


DETAIL_COLUMNS = [
    "dataset", "method", "metric", "source", "target", "score",
    "evaluable_source_cells",
]
QA_COLUMNS = ["dataset", "method", "metric", "raw_values_match", "aggregate_abs_difference"]
EXTENDED_DETAIL_COLUMNS = [
    "dataset", "method", "metric", "source", "target", "score",
    "evaluable_source_cells", "raw_nonfinite_source_cells",
    "branch_winner_source_cells",
]
FOREBRAIN_DETAIL_COLUMNS = [
    "dataset", "config_id", "method", "source", "target", "score",
    "evaluable_source_cells", "included_in_mean",
]


def evaluation_protocol(row: dict, config: dict, file_hashes: dict) -> dict:
    """Bind a pair cache to source bytes and the complete scientific protocol."""
    dataset = row["evaluation_dataset"]
    cbdir = config["datasets"][dataset].get("cbdir", {})
    space = cbdir.get("embedding", "umap")
    protocol = {
        "config_id": row["config_id"],
        "source_sha256": row["sha256"],
        "velocity_layer_key": row["velocity_layer_key"],
        "dataset": dataset,
        "cbdir_space": space,
        "cbdir_dimensions": cbdir.get("dimensions", 2),
        "frame_sha256": cbdir.get("coordinate_sha256", "registered_reference_umap"),
        "files": file_hashes,
    }
    protocol["sha256"] = hashlib.sha256(
        json.dumps(protocol, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return protocol


def pair_is_current(payload: dict, protocol: dict) -> bool:
    return (
        protocol is not None
        and payload.get("status") == "ok"
        and "extended_edge_details" in payload
        and payload.get("evaluation_protocol") == protocol
    )


def read_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"operational_id", "evaluation_dataset", "result_path"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path}: missing required manifest columns")
    seen = set()
    for row in rows:
        key = (row["evaluation_dataset"], row["operational_id"])
        if key in seen:
            raise ValueError(f"duplicate manifest input: {key}")
        seen.add(key)
    return rows


def numeric_sequences_match(left, right, tolerance: float = 1e-12) -> bool:
    if len(left) != len(right):
        return False
    return all(abs(float(a) - float(b)) <= tolerance for a, b in zip(left, right))


def compare_tivelo(indicators: dict, official: dict) -> list[dict]:
    result = []
    for metric in ("highdim_cbdir", "transcosine"):
        raw_key = f"{metric}_raw"
        left = indicators["metrics"][raw_key]
        right = official[raw_key]
        keys_match = set(left) == set(right)
        raw_match = keys_match and all(
            numeric_sequences_match(left[key], right[key]) for key in left
        )
        frozen_aggregate = indicators["metrics"][metric]
        official_aggregate = official[metric]
        if frozen_aggregate is None and official_aggregate is None:
            aggregate_difference = 0.0
        elif frozen_aggregate is None or official_aggregate is None:
            raise ValueError(
                f"{metric}: frozen/official mismatch "
                f"(one aggregate is None, other={official_aggregate!r})"
            )
        else:
            aggregate_difference = abs(
                float(frozen_aggregate) - float(official_aggregate)
            )
        # Raw TiVelo outputs must match. Aggregate precision differs because
        # TiVelo's numpy mean preserves a float32 intermediate for TransCosine.
        if not raw_match or aggregate_difference > 1e-7:
            raise ValueError(
                f"{metric}: frozen/official mismatch "
                f"(raw_match={raw_match}, aggregate_difference={aggregate_difference})"
            )
        result.append(
            {
                "metric": metric,
                "raw_values_match": raw_match,
                "aggregate_abs_difference": aggregate_difference,
            }
        )
    return result


def edge_records(dataset: str, method: str, official: dict) -> list[dict]:
    records = []
    for metric in ("highdim_cbdir", "transcosine"):
        for edge, values in official[f"{metric}_raw"].items():
            source, target = edge.split("->", 1)
            records.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "metric": metric,
                    "source": source,
                    "target": target,
                    "score": (sum(values) / len(values)) if values else None,
                    "evaluable_source_cells": len(values),
                }
            )
    return records


def extended_edge_records(dataset: str, method: str, indicators: dict) -> list[dict]:
    records = []
    for edge, values in indicators["metrics"]["cbdir_edge_raw"].items():
        source, target = edge.split("->", 1)
        nonfinite = [
            value is None or not math.isfinite(float(value)) for value in values
        ]
        records.append(
            {
                "dataset": dataset,
                "method": method,
                "metric": "cbdir",
                "source": source,
                "target": target,
                "score": (
                    sum(values) / len(values)
                    if values and not any(nonfinite) else None
                ),
                "evaluable_source_cells": len(values),
                "raw_nonfinite_source_cells": sum(nonfinite),
                "branch_winner_source_cells": "",
            }
        )
    for item in indicators["metrics"]["batc_edge_details"]:
        records.append(
            {
                "dataset": dataset,
                "method": method,
                "metric": "batc_candidate_source_mean",
                "source": item["source"],
                "target": item["target"],
                "score": item["candidate_mean_source_cells"],
                "evaluable_source_cells": item["source_cells"],
                "raw_nonfinite_source_cells": item["candidate_nonfinite_source_cells"],
                "branch_winner_source_cells": item["branch_winner_source_cells"],
            }
        )
    return records


def build_pair_payload(dataset: str, method: str, indicators: dict, official: dict, qa: list[dict]) -> dict:
    metadata = indicators["metadata"]
    metric = indicators["metrics"]
    trans_raw = official["transcosine_raw"]
    highdim_raw = official["highdim_cbdir_raw"]
    return {
        "dataset": dataset,
        "method": method,
        "status": "ok",
        "metrics": {
            "consistency": metric["consistency"],
            "iccoh": metric["iccoh"],
            "cbdir": metric["cbdir"],
            "transcosine": official["transcosine"],
            "batc": metric["batc"],
            "highdim_cbdir": official["highdim_cbdir"],
        },
        "counts": {
            "transcosine_evaluable_edges": sum(bool(values) for values in trans_raw.values()),
            "transcosine_evaluable_source_cells": sum(len(values) for values in trans_raw.values()),
            "highdim_cbdir_evaluable_edges": sum(bool(values) for values in highdim_raw.values()),
            "highdim_cbdir_evaluable_source_cells": sum(len(values) for values in highdim_raw.values()),
        },
        "metadata": metadata,
        "tivelo_parity": qa,
        "edge_details": edge_records(dataset, method, official),
        "extended_edge_details": extended_edge_records(dataset, method, indicators),
    }


def run(command: list[str], log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}); see {log_path}")


def write_forebrain_cbdir_audit(pairs: list[dict], result_root: Path, config: dict) -> None:
    """Bind a complete Forebrain rerun's edge scores to its current inputs."""
    selected = [pair for pair in pairs if pair["dataset"] == "forebrain"]
    if not selected or {pair["method"] for pair in selected} != set(config.get("method_order", [])):
        for name in ("forebrain_cbdir_pca_edge_details.csv", "forebrain_cbdir_pca_manifest.json"):
            path = result_root / name
            if path.is_file():
                path.unlink()
        return
    if len(selected) != len(config["method_order"]):
        raise ValueError("Forebrain audit requires unique method pairs")
    frames = selected[0]["metadata"]["metric_frames"]
    files = selected[0]["evaluation_protocol"]["files"]
    cbdir = config["datasets"]["forebrain"]["cbdir"]
    if (frames["cbdir"]["space"], frames["cbdir"]["dimensions"],
        frames["cbdir"]["coordinate_sha256"]) != (
        "pca", 30, cbdir["coordinate_sha256"],
    ) or (frames["batc"]["space"], frames["batc"]["dimensions"]) != ("umap", 2):
        raise ValueError("Forebrain audit requires the frozen PCA30 and UMAP2 frames")
    edges = config["datasets"]["forebrain"]["cluster_edges"]
    expected = {(str(source), str(target)) for source, target in edges}
    records, inputs = [], []
    for pair in sorted(selected, key=lambda item: item["evaluation_protocol"]["config_id"]):
        metadata, protocol = pair["metadata"], pair["evaluation_protocol"]
        if metadata["metric_frames"] != frames or protocol["files"] != files:
            raise ValueError("Forebrain pairs have different frames or evaluation inputs")
        if metadata["result_sha256"] != protocol["source_sha256"]:
            raise ValueError("Forebrain source bytes differ from the evaluated protocol")
        current = [row for row in pair["extended_edge_details"] if row["metric"] == "cbdir"]
        if len(current) != len(expected) or {
            (str(row["source"]), str(row["target"])) for row in current
        } != expected:
            raise ValueError("Forebrain audit requires each configured CBDir edge once")
        values = []
        for row in current:
            score = row["score"]
            included = score is not None and score != 0.0
            if included:
                if not math.isfinite(float(score)) or row["evaluable_source_cells"] <= 0:
                    raise ValueError("Forebrain CBDir cannot include an undefined edge")
                values.append(float(score))
            records.append({"dataset": "forebrain", "config_id": protocol["config_id"],
                "method": pair["method"], "source": str(row["source"]),
                "target": str(row["target"]), "score": score,
                "evaluable_source_cells": row["evaluable_source_cells"],
                "included_in_mean": included})
        if not values or not math.isclose(sum(values) / len(values),
            pair["metrics"]["cbdir"], rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("Forebrain edge means differ from the primary CBDir score")
        inputs.append({"method": pair["method"], "config_id": protocol["config_id"],
            "locator": metadata["result_path"], "sha256": metadata["result_sha256"],
            "velocity_layer_key": metadata["velocity_layer_key"]})
    edge_path = result_root / "forebrain_cbdir_pca_edge_details.csv"
    with edge_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FOREBRAIN_DETAIL_COLUMNS)
        writer.writeheader()
        writer.writerows(records)
    manifest = {
        "schema_version": 1, "dataset": "forebrain", "status": "reproduced_pca30_cbdir",
        "metric": "cbdir", "coordinate_space": "pca", "dimensions": 30,
        "pca_coordinate_sha256": frames["cbdir"]["coordinate_sha256"],
        "umap_coordinate_sha256": frames["batc"]["coordinate_sha256"],
        "projection": {"basis": "pca", "direct_pca_projection": False,
                       "autoscale": True, "all_comps": True},
        "evaluation_files": files,
        "scvelo_versions": sorted({pair["metadata"]["scvelo_version"] for pair in selected}),
        "inputs": inputs,
        "cbdir_aggregation": "UnitVelo edge means; exact-zero and empty edges excluded; remaining edges equally weighted",
        "edge_details_file": edge_path.name,
        "edge_details_sha256": hashlib.sha256(edge_path.read_bytes()).hexdigest(),
        "edge_rows": len(records), "primary_table": "cell_metrics.csv",
        "primary_table_sha256": hashlib.sha256((result_root / "cell_metrics.csv").read_bytes()).hexdigest(),
        "guard_tables": {name: hashlib.sha256((result_root / name).read_bytes()).hexdigest()
            for name in ("cell_metrics_supplementary.csv", "cell_metrics_audit.csv",
                         "cell_metrics_edge_details.csv", "tivelo_parity.csv")},
        "randomness_policy": "Frozen coordinate and neighbor assets; no PCA/UMAP construction or training",
    }
    (result_root / "forebrain_cbdir_pca_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def write_tables(pair_root: Path, result_root: Path, config: dict, protocols: dict = None) -> None:
    pairs = []
    for path in sorted(pair_root.glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") == "ok":
            key = (payload["dataset"], payload["method"])
            if protocols is not None and not pair_is_current(payload, protocols.get(key)):
                continue
            pairs.append(payload)
    primary_rows = []
    supplementary_rows = []
    audit_rows = []
    details = []
    extended_details = []
    qa_rows = []
    for payload in pairs:
        metrics = payload["metrics"]
        counts = payload["counts"]
        metadata = payload["metadata"]
        row = {
            "dataset": payload["dataset"], "method": payload["method"],
            "consistency": metrics["consistency"], "iccoh": metrics["iccoh"],
            "cbdir": metrics["cbdir"], "transcosine": metrics["transcosine"],
            "batc": metrics["batc"], "highdim_cbdir": metrics["highdim_cbdir"],
            **counts,
            "result_path": metadata["result_path"],
            "result_sha256": metadata["result_sha256"],
            "velocity_layer_key": metadata["velocity_layer_key"],
            "invalid_velocity_values_replaced_with_zero": metadata[
                "invalid_velocity_values_replaced_with_zero"
            ],
            "status": payload["status"],
        }
        primary_rows.append(
            {column: row[column] for column in config["result_tables"]["primary"]}
        )
        supplementary_rows.append(
            {column: row[column] for column in config["result_tables"]["supplementary"]}
        )
        audit_rows.append(
            {column: row[column] for column in config["result_tables"]["audit"]}
        )
        details.extend(payload["edge_details"])
        extended_details.extend(payload["extended_edge_details"])
        qa_rows.extend(
            {"dataset": payload["dataset"], "method": payload["method"], **item}
            for item in payload["tivelo_parity"]
        )
    result_root.mkdir(parents=True, exist_ok=True)
    for filename, fieldnames, rows in (
        ("cell_metrics.csv", config["result_tables"]["primary"], primary_rows),
        ("cell_metrics_supplementary.csv", config["result_tables"]["supplementary"], supplementary_rows),
        ("cell_metrics_audit.csv", config["result_tables"]["audit"], audit_rows),
        ("cell_metrics_edge_details.csv", DETAIL_COLUMNS, details),
        ("cell_metrics_edge_details_extended.csv", EXTENDED_DETAIL_COLUMNS, extended_details),
        ("tivelo_parity.csv", QA_COLUMNS, qa_rows),
    ):
        with (result_root / filename).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    write_forebrain_cbdir_audit(pairs, result_root, config)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument(
        "--inputs",
        type=Path,
        default=DEFAULT_INPUTS,
        help="logical Cell input manifest (default: evaluation/cell/inputs.tsv)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Cell evaluation config (default: evaluation/cell/config.json)",
    )
    parser.add_argument(
        "--reference-neighbors-root",
        type=Path,
        default=MODULE_DIR / "reference_neighbors",
        help="fixed reference-30 assets; defaults to the repository copy",
    )
    parser.add_argument("--reference-moments-root", type=Path, required=True)
    parser.add_argument("--reference-embeddings-root", type=Path, default=MODULE_DIR / "reference_embeddings")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path)
    parser.add_argument("--projection-root", type=Path)
    parser.add_argument("--indicators-python", type=Path, required=True)
    parser.add_argument("--tivelo-python", type=Path, required=True)
    parser.add_argument("--prepare-script", type=Path, required=True)
    parser.add_argument("--official-script", type=Path, required=True)
    parser.add_argument("--method", action="append")
    parser.add_argument("--dataset", action="append")
    parser.add_argument(
        "--n-jobs",
        type=int,
        help="velocity_graph workers; defaults to config.json",
    )
    parser.add_argument(
        "--pair-workers", type=int, default=1,
        help="number of dataset/method combinations to run concurrently",
    )
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    all_rows = read_rows(args.inputs)
    rows = all_rows
    config = json.loads(args.config.read_text(encoding="utf-8"))
    configured_n_jobs = int(config["metrics"]["velocity_graph_n_jobs"])
    n_jobs = args.n_jobs if args.n_jobs is not None else configured_n_jobs
    if args.method:
        rows = [row for row in rows if row["operational_id"] in set(args.method)]
    if args.dataset:
        rows = [row for row in rows if row["evaluation_dataset"] in set(args.dataset)]
    if not rows:
        raise ValueError("no rows selected")

    output_root = args.output_root.resolve()
    work_root = (
        args.work_root.resolve()
        if args.work_root is not None
        else output_root / ".work"
    )
    pair_root = work_root / "pairs"
    log_root = work_root / "logs"
    reference_umap_cache = work_root / "reference_umap_scvelo_deterministic"
    result_root = output_root / "results"
    if result_root == MODULE_DIR / "results":
        expected = {(dataset, method) for dataset in config["datasets"] for method in config["method_order"]}
        selected_keys = {(row["evaluation_dataset"], row["operational_id"]) for row in rows}
        if selected_keys != expected:
            raise ValueError("A partial run must use a separate --output-root, not the canonical Cell bundle")
    for directory in (pair_root, work_root, log_root, reference_umap_cache):
        directory.mkdir(parents=True, exist_ok=True)

    def file_sha(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    shared_hashes = {
        "config": file_sha(args.config), "inputs": file_sha(args.inputs),
        "prepare_script": file_sha(args.prepare_script),
        "official_script": file_sha(args.official_script),
        "frame_loader": file_sha(MODULE_DIR / "reference_embeddings.py"),
        "metrics": file_sha(MODULE_DIR / "metrics.py"),
        "batc": file_sha(MODULE_DIR / "batc.py"),
        "reference_neighbors_code": file_sha(MODULE_DIR / "reference_neighbors.py"),
    }
    dataset_hashes = {}
    for dataset in {row["evaluation_dataset"] for row in rows}:
        files = {
            **shared_hashes,
            "reference_neighbors": file_sha(args.reference_neighbors_root / f"{dataset}.npz"),
            "reference_moments": file_sha(args.reference_moments_root / f"{dataset}.npz"),
        }
        if dataset == "forebrain":
            files["frame_manifest"] = file_sha(args.reference_embeddings_root / "manifest.json")
            files["umap_frame"] = file_sha(args.reference_embeddings_root / "forebrain_umap.npz")
        dataset_hashes[dataset] = files
    protocols = {
        (row["evaluation_dataset"], row["operational_id"]): evaluation_protocol(
            row, config, dataset_hashes[row["evaluation_dataset"]]
        ) for row in rows
    }

    def process_row(position: int, row: dict) -> None:
        method = row["operational_id"]
        dataset = row["evaluation_dataset"]
        pair_path = pair_root / f"{dataset}__{method}.json"
        protocol = protocols[(dataset, method)]
        if args.resume and pair_path.is_file():
            payload = json.loads(pair_path.read_text(encoding="utf-8"))
            if pair_is_current(payload, protocol):
                print(f"[{position}/{len(rows)}] skip {dataset}/{method}", flush=True)
                return
        stem = f"{dataset}__{method}"
        prepared = work_root / f"{stem}.h5ad"
        indicators_json = work_root / f"{stem}.indicators.json"
        official_json = work_root / f"{stem}.official.json"
        try:
            print(f"[{position}/{len(rows)}] {dataset}/{method}", flush=True)
            run(
                [
                    str(args.indicators_python), str(args.prepare_script),
                    "--artifact-root", str(args.artifact_root), "--inputs", str(args.inputs),
                    "--config", str(args.config),
                    "--reference-neighbors-root", str(args.reference_neighbors_root),
                    "--reference-moments-root", str(args.reference_moments_root),
                    "--reference-embeddings-root", str(args.reference_embeddings_root),
                    "--reference-umap-cache", str(reference_umap_cache),
                    "--method", method, "--dataset", dataset,
                    "--prepared-output", str(prepared),
                    "--indicators-output", str(indicators_json),
                    "--n-jobs", str(n_jobs),
                    *(["--projection-output", str(args.projection_root / f"{dataset}__{method}.npz")]
                      if args.projection_root else []),
                ],
                log_root / f"{stem}.indicators.log",
            )
            run(
                [
                    str(args.tivelo_python), str(args.official_script),
                    "--prepared-input", str(prepared), "--config", str(args.config),
                    "--dataset", dataset, "--output", str(official_json),
                ],
                log_root / f"{stem}.tivelo.log",
            )
            indicators = json.loads(indicators_json.read_text(encoding="utf-8"))
            official = json.loads(official_json.read_text(encoding="utf-8"))
            qa = compare_tivelo(indicators, official)
            payload = build_pair_payload(dataset, method, indicators, official, qa)
            payload["evaluation_protocol"] = protocol
            pair_path.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            print(
                f"  ok: transcosine={payload['metrics']['transcosine']:.6g}, "
                f"highdim_cbdir={payload['metrics']['highdim_cbdir']:.6g}",
                flush=True,
            )
        except Exception as error:
            pair_path.write_text(
                json.dumps(
                    {"dataset": dataset, "method": method, "status": "error", "error": str(error)},
                    indent=2, sort_keys=True,
                ) + "\n",
                encoding="utf-8",
            )
            print(f"  ERROR: {error}", file=sys.stderr, flush=True)
            raise
        finally:
            for temporary_path in (prepared, indicators_json, official_json):
                if temporary_path.exists():
                    temporary_path.unlink()

    indexed_rows = list(enumerate(rows, start=1))
    if args.pair_workers <= 1:
        for position, row in indexed_rows:
            process_row(position, row)
    else:
        failures = []
        with ThreadPoolExecutor(max_workers=args.pair_workers) as executor:
            futures = {
                executor.submit(process_row, position, row): (position, row)
                for position, row in indexed_rows
            }
            for future in as_completed(futures):
                position, row = futures[future]
                try:
                    future.result()
                except Exception as error:
                    failures.append((position, row, error))
        if failures:
            position, row, error = sorted(failures, key=lambda item: item[0])[0]
            raise RuntimeError(
                f"parallel evaluation failed at [{position}] "
                f"{row['evaluation_dataset']}/{row['operational_id']}: {error}"
            ) from error
    write_tables(pair_root, result_root, config, protocols)
    print(f"Wrote {len(rows)} selected rows to {result_root}", flush=True)


if __name__ == "__main__":
    main()
