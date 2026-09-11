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
import json
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
    }


def run(command: list[str], log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}); see {log_path}")


def write_tables(pair_root: Path, result_root: Path, config: dict) -> None:
    pairs = []
    for path in sorted(pair_root.glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") == "ok":
            pairs.append(payload)
    primary_rows = []
    supplementary_rows = []
    audit_rows = []
    details = []
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
        ("tivelo_parity.csv", QA_COLUMNS, qa_rows),
    ):
        with (result_root / filename).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)


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
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path)
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

    rows = read_rows(args.inputs)
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
    for directory in (pair_root, work_root, log_root, reference_umap_cache):
        directory.mkdir(parents=True, exist_ok=True)

    def process_row(position: int, row: dict) -> None:
        method = row["operational_id"]
        dataset = row["evaluation_dataset"]
        pair_path = pair_root / f"{dataset}__{method}.json"
        if args.resume and pair_path.is_file():
            payload = json.loads(pair_path.read_text(encoding="utf-8"))
            if payload.get("status") == "ok":
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
                    "--reference-umap-cache", str(reference_umap_cache),
                    "--method", method, "--dataset", dataset,
                    "--prepared-output", str(prepared),
                    "--indicators-output", str(indicators_json),
                    "--n-jobs", str(n_jobs),
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
    write_tables(pair_root, result_root, config)
    print(f"Wrote {len(rows)} selected rows to {result_root}", flush=True)


if __name__ == "__main__":
    main()
