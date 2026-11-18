"""Reproduce S9 with the registered Gene-matched Cell output for each fit."""

import argparse
import csv
import errno
import os
from decimal import Decimal
from itertools import combinations_with_replacement
from math import factorial, isfinite, prod
from pathlib import Path
from statistics import mean
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from evaluation.gene.config import METHOD_ORDER, REAL_CONFIGS
from evaluation.gene.refresh_common_gene_summary import common_gene_tables

CELL_FIELDS = ("consistency", "iccoh", "cbdir", "batc", "transcosine")
OBJECTIVES = ("gene_cbdir", "cell_cbdir", "batc")
MATCHED = {"C04", "C05", "C12", "C16"}


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t" if path.suffix == ".tsv" else ","))


def index(rows, *fields):
    indexed = {tuple(row[field].casefold() for field in fields): row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"Duplicate key in {fields}")
    return indexed


def matched_scores(root):
    gene = root / "evaluation/gene/results/gene_cbdir_reference30_all.csv"
    tables = common_gene_tables(gene)
    selected = {(r["dataset"], r["gene"]) for r in tables["membership.csv"]}
    sources = {}
    for row in read(gene):
        if (row["dataset"], row["gene"]) in selected:
            key = (row["dataset"], row["method"])
            path = row["source_h5ad"].replace("artifact://velocity/", "velocity:", 1)
            if sources.setdefault(key, path) != path or int(row["n_evaluable_edges"]) < 1:
                raise ValueError(f"Inconsistent Gene input for {key}")
    registry = index(read(root / "methods/configurations.tsv"), "manuscript_name")
    cell = root / "evaluation/cell"
    routes = {}
    for route, directory in (("primary_cell", cell), ("gene_matched", cell / "gene_matched")):
        inputs = read(directory / "inputs.tsv")
        if route == "gene_matched":
            inputs = [r for r in inputs if r["evaluation_role"] == "gene_matched_cell"]
        routes[route] = (index(inputs, "dataset_id", "config_id"),
                         index(read(directory / "results/cell_metrics.csv"), "dataset", "method"),
                         index(read(directory / "results/cell_metrics_audit.csv"), "dataset", "method"))
    native = index([r for r in read(cell / "native_outputs.tsv") if r["scope"] == "real"],
                   "evaluation_dataset")
    scores = []
    for common in tables["dataset_summary.csv"]:
        dataset, label = common["dataset"], common["method"]
        method = registry[("multivelo (rna_only)" if label == "MultiVelo" else label.casefold(),)]
        cid = method["config_id"]
        route = "gene_matched" if cid in MATCHED else "primary_cell"
        inputs, metrics, audit = routes[route]
        source = inputs[dataset.casefold(), cid.casefold()]
        cell_dataset = source["evaluation_dataset"]
        key = cell_dataset.casefold(), method["operational_id"].casefold()
        score, record = metrics[key], audit[key]
        if (record["status"] != "ok" or record["result_path"] != source["result_path"]
                or record["result_sha256"] != source["sha256"]
                or record["velocity_layer_key"] != source["velocity_layer_key"]):
            raise ValueError(f"Cell audit differs for {key}")
        gene_path = sources[dataset, label]
        relation = "same_result_file"
        if cid == "C04":
            bridge = native[cell_dataset.casefold(),]
            if (bridge["fit_result_path"] != gene_path or bridge["result_path"] != source["result_path"]
                    or bridge["sha256"] != source["sha256"] or source["velocity_layer_key"] != "velo_s_norm"):
                raise ValueError(f"Missing MultiVelo fit bridge for {dataset}")
            relation = "native_fit_to_normalized_cell_output"
        elif gene_path != source["result_path"]:
            raise ValueError(f"Different Gene and Cell fits for {key}")
        if not all(isfinite(float(score[f])) for f in CELL_FIELDS):
            raise ValueError(f"Missing Cell score for {key}")
        scores.append(dict(dataset=dataset, cell_dataset=cell_dataset, n_common_genes=common["n_common_genes"],
                           method=method["manuscript_name"], config_id=cid, operational_id=method["operational_id"],
                           gene_cbdir=common["mean"], cell_cbdir=score["cbdir"], batc=score["batc"],
                           cell_score_route=route, gene_result_path=gene_path, cell_result_path=source["result_path"],
                           cell_result_sha256=source["sha256"], fit_relation=relation,
                           consistency=score["consistency"], iccoh=score["iccoh"], transcosine=score["transcosine"]))
    return scores


def front(vectors):
    return [not any(all(x >= y for x, y in zip(other, vector)) and other != vector
                    for other in vectors) for vector in vectors]


def resampled_fronts(matrix, counts):
    # The positive support supplies the Gene objectives. Joint objectives are
    # paired weighted sums; the common positive denominator cannot affect dominance.
    gene, joint = [], []
    for rows in matrix:
        gene.append(tuple(r[0] for r, c in zip(rows, counts) if c))
        joint.append(tuple(sum(c * r[j] for c, r in zip(counts, rows)) for j in range(3)))
    return front(gene), front(joint)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=ROOT)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    if os.path.lexists(args.output_root):
        parser.error(f"--output-root must be a new directory: {args.output_root}")
    scores = matched_scores(args.repository_root)
    datasets = list(REAL_CONFIGS)
    methods = list(dict.fromkeys(r["method"] for r in scores))
    by_key = {(r["dataset"], r["method"]): r for r in scores}
    if len(datasets) != 7 or len(methods) != len(METHOD_ORDER) or len(by_key) != 63 or len(scores) != 63:
        raise ValueError("S9 requires seven datasets and nine methods")
    data = [[[Decimal(str(by_key[d, m][f])) for f in OBJECTIVES] for d in datasets] for m in methods]
    # Integerize the recorded decimals exactly to preserve equality on resampling.
    scale = 10 ** max(-x.as_tuple().exponent for rows in data for row in rows for x in row)
    matrix = [[[int(x * scale) for x in row] for row in rows] for rows in data]
    n = len(datasets)
    observed = resampled_fronts(matrix, [1] * n)
    hits = [[0] * len(methods) for _ in observed]
    total = 0
    for draw in combinations_with_replacement(range(n), n):
        counts = [draw.count(i) for i in range(n)]
        weight = factorial(n) // prod(factorial(c) for c in counts)
        total += weight
        for tally, flags in zip(hits, resampled_fronts(matrix, counts)):
            for i, included in enumerate(flags):
                tally[i] += weight * included
    if total != n ** n:
        raise ValueError("Incomplete enumeration of paired dataset resamples")
    names = ["gene", "joint"]
    frontier, omitted, plotted, means = [], [], [], []
    display = index(read(args.repository_root / "data/datasets.tsv"), "dataset_id")
    for j, d in enumerate(datasets):
        flags = resampled_fronts(matrix, [int(i != j) for i in range(n)])
        for k, setting in enumerate(names):
            for i, m in enumerate(methods):
                omitted.append(dict(setting=setting, omitted_dataset=d, method=m, on_front=flags[k][i]))
    for i, m in enumerate(methods):
        for panel, field in zip("abcdef", (*CELL_FIELDS, "gene_cbdir")):
            source = "cell_cbdir" if field == "cbdir" else field
            values = [float(by_key[d, m][source]) for d in datasets]
            means.append(dict(panel=panel, method=m, field=field, mean=mean(values)))
            for d, value in zip(datasets, values):
                label = d if panel == "f" else display[d.casefold(),]["paper_name"]
                plotted.append(dict(panel=panel, method=m, dataset=label, field=field, value=value))
        for f in OBJECTIVES:
            value = mean(float(by_key[d, m][f]) for d in datasets)
            plotted.append(dict(panel="g", method=m, dataset="macro", field="mean_" + f, value=value))
        for k, setting in enumerate(names):
            freq = hits[k][i] / total
            frontier.append(dict(setting=setting, method=m, on_front=observed[k][i],
                                 bootstrap_front_hits=hits[k][i], bootstrap_ordered_draws=total,
                                 bootstrap_frequency=freq))
            plotted.append(dict(panel="h", method=m, dataset="bootstrap", field=setting + "_front_percent",
                                value=100 * freq))
    args.output_root.parent.mkdir(exist_ok=True, parents=True)
    with TemporaryDirectory(dir=args.output_root.parent, prefix=".s9-") as temp:
        staged = Path(temp) / "tables"
        staged.mkdir()
        for name, rows in (("dataset_scores", scores), ("frontiers", frontier), ("leave_one_out", omitted),
                           ("plot_data", plotted), ("upper_means", means)):
            with (staged / (name + ".csv")).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        if os.path.lexists(args.output_root):
            parser.error(f"--output-root must be a new directory: {args.output_root}")
        try:
            staged.rename(args.output_root)
        except OSError as exc:
            # Another replay can publish its populated directory between the check and rename.
            if exc.errno in (errno.ENOTEMPTY, errno.EEXIST):
                parser.error(f"--output-root must be a new directory: {args.output_root}")
            raise


if __name__ == "__main__":
    main()
