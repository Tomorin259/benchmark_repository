#!/usr/bin/env python3
"""Rank all Cell methods by paired-bootstrap inclusion in the Pareto front."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ensemble.core import sha256

ROOT = Path(__file__).resolve().parents[1]
METRICS = ("cbdir", "batc")
SOURCE_NAMES = ("evaluation/cell/results/cell_metrics.csv", "evaluation/cell/inputs.tsv",
                "evaluation/cell/results/cell_metrics_audit.csv", "methods/configurations.tsv")


def dominance(values):
    """Match the original five-objective analysis, retaining paired arrays."""
    values = np.asarray(values, dtype=np.float64)
    left, right = values[..., :, None, :], values[..., None, :, :]
    return np.all(left >= right, axis=-1) & np.any(left > right, axis=-1)


def load_formal_scores():
    """Read the registered score grid and verify its source audit."""
    def read_rows(name):
        with (ROOT / name).open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t" if name.endswith(".tsv") else ","))
    inputs = read_rows(SOURCE_NAMES[1])
    registry = {row["config_id"]: row for row in read_rows(SOURCE_NAMES[3])}
    configuration = {row["operational_id"]: row["config_id"] for row in inputs}
    methods = sorted(configuration, key=configuration.get)
    rows = read_rows(SOURCE_NAMES[0])
    datasets = sorted({row["dataset"] for row in rows})
    by_pair = {(row["dataset"], row["method"]): row for row in rows}
    expected = {(row["evaluation_dataset"], row["operational_id"]) for row in inputs}
    if (len(rows) != 280 or len(by_pair) != 280 or set(by_pair) != expected
            or len(methods) != 20 or len(datasets) != 14):
        raise ValueError("use the complete formal 20-method by 14-dataset score grid")
    audits = {(row["dataset"], row["method"]): row for row in read_rows(SOURCE_NAMES[2])}
    for row in inputs:
        audit = audits[row["evaluation_dataset"], row["operational_id"]]
        if (audit["status"] != "ok" or audit["result_sha256"] != row["sha256"]
                or audit["result_path"] != row["result_path"]):
            raise ValueError("score audit and registered source differ")
    configs = [
        {"config_id": configuration[method], "method": method,
         "manuscript_name": registry[configuration[method]]["manuscript_name"]}
        for method in methods
    ]
    return rows, configs


def bootstrap_rankings(rows, configs, resamples, seed):
    """Resample paired dataset scores and rank the full method roster."""
    datasets = sorted({row["dataset"] for row in rows})
    by_pair = {(row["dataset"], row["method"]): row for row in rows}
    methods = [config["method"] for config in configs]
    values = np.asarray([[[float(by_pair[dataset, method][metric]) for metric in METRICS]
                          for method in methods] for dataset in datasets])
    if not np.isfinite(values).all():
        raise ValueError("all scores must be finite")
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(datasets), size=(resamples, len(datasets)))
    membership = np.empty((resamples, len(methods)), dtype=bool)
    for start in range(0, resamples, 250):
        end = min(start + 250, resamples)
        means = values[draws[start:end]].mean(axis=1)
        membership[start:end] = ~dominance(means).any(axis=1)
    original_means = values.mean(axis=0)
    original_front = ~dominance(original_means).any(axis=0)
    counts = membership.sum(axis=0)
    order = sorted(range(20), key=lambda index: (-counts[index], configs[index]["config_id"]))
    ranked = []
    for index in order:
        frequency = counts[index] / resamples
        ranked.append({
            "rank": 1 + int((counts > counts[index]).sum()),
            **configs[index],
            "bootstrap_front_count": int(counts[index]),
            "bootstrap_front_percent": float(100 * frequency),
            "bootstrap_frequency_mc_se_percent": float(100 * np.sqrt(frequency * (1 - frequency) / resamples)),
            "original_mean_front": bool(original_front[index]),
            "cbdir_mean": float(original_means[index, 0]), "batc_mean": float(original_means[index, 1]),
        })
    return ranked, draws, membership, datasets, values


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--resamples", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261109)
    args = parser.parse_args(argv)
    if args.resamples < 1 or args.output_root.exists():
        raise ValueError("use a positive resample count and a new output directory")
    sources = {name: sha256(ROOT / name) for name in SOURCE_NAMES}
    rows, configs = load_formal_scores()
    ranked, draws, membership, datasets, values = bootstrap_rankings(rows, configs, args.resamples, args.seed)
    methods = [config["method"] for config in configs]
    output = args.output_root
    output.mkdir(parents=True)
    with (output / "ranking.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ranked[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(ranked)
    np.savez_compressed(output / "bootstrap.npz", dataset_indices=draws, front_membership=membership,
                        datasets=np.asarray(datasets), methods=np.asarray(methods), scores=values)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(8.4, 7.5), layout="constrained")
    percents = [row["bootstrap_front_percent"] for row in ranked]
    ax.barh(np.arange(20), percents,
            color=["#0072B2" if row["original_mean_front"] else "#8D99A5" for row in ranked])
    ax.set_yticks(np.arange(20), [row["manuscript_name"] for row in ranked])
    ax.invert_yaxis()
    ax.set_xlim(0, max(percents) * 1.17)
    for index, value in enumerate(percents):
        ax.text(value + max(percents) * .015, index, f"{value:.2f}%", va="center", fontsize=9)
    ax.set_xlabel("Bootstrap inclusion in the CBDir–BATC Pareto front (%)")
    ax.set_title(f"20 configurations | 14 datasets | {args.resamples:,} paired draws", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for extension in ("png", "pdf"):
        fig.savefig(output / ("ranking." + extension), dpi=180, facecolor="white")
    plt.close(fig)
    for name, digest in sources.items():
        if sha256(ROOT / name) != digest:
            raise ValueError("formal source changed during the analysis")
    manifest = {
        "status": "complete", "metrics": list(METRICS), "resamples": args.resamples, "seed": args.seed,
        "bootstrap_unit": "dataset", "draw_size": 14, "replacement": True,
        "pairing": "every draw uses the same dataset indices for all 20 methods and both metrics",
        "statistic": "percentage of bootstrap macro-mean pairs on the nondominated front",
        "dataset_weighting": "equal within each draw", "directions": "maximize both",
        "dominance": "all A >= B and at least one A > B; full precision; epsilon=0",
        "interpretation": "empirical dataset-reweighting stability, not probability of being correct",
        "score_space": "current formal method-native Cell scores, including updated GraphVelo bone marrow",
        "repository_head": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "generator_sha256": sha256(Path(__file__)), "sources_sha256": sources,
        "numpy_version": np.__version__, "rng": "default_rng / PCG64", "batch_size": 250,
        "front_size_mean": float(membership.sum(axis=1).mean()),
        "caption": "Methods ranked by inclusion frequency in the Pareto front of macro-mean CBDir and BATC "
                   "over 10,000 paired dataset-bootstrap draws. Each draw samples 14 of the observed datasets "
                   "with replacement and retains all methods and both metrics together. Blue bars denote "
                   "methods on the original mean front. Percentages need not sum to 100 because multiple "
                   "methods can share a front. Frequencies describe stability on these observed datasets.",
        "outputs_sha256": {path.name: sha256(path) for path in output.iterdir() if path.is_file()},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    for row in ranked:
        print(f"{row['rank']:2d} {row['manuscript_name']:<24} {row['bootstrap_front_percent']:.2f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

