"""Refresh Figure-summary Real-CBDir values from the reference-30 gene results."""

import argparse
import csv
import io
import math
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "evaluation/summary/results/figure6_summary"


def read_csv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def main(write):
    with (ROOT / "data/datasets.tsv").open(newline="") as handle:
        dataset_registry = list(csv.DictReader(handle, delimiter="\t"))
    dataset_names = {row["paper_name"].casefold(): row["dataset_id"].casefold()
                     for row in dataset_registry}

    wps_path = RESULTS / "fig2_gene_metric_wps_values.csv"
    wps_fields, wps_rows = read_csv(wps_path)
    method_keys = {row["Figure2_Method"]: row["Source_Method_Key"].casefold()
                   for row in wps_rows if row["Source_Method_Key"]}

    source_path = ROOT / "evaluation/gene/results/gene_cbdir_reference30_all.csv"
    _, source_rows = read_csv(source_path)
    scores = defaultdict(list)
    for row in source_rows:
        if row["status"] == "ok":
            scores[(row["dataset"].casefold(), row["method"].casefold())].append(
                float(row["new_cbdir"])
            )

    _, canonical_summary = read_csv(
        ROOT / "evaluation/gene/results/gene_cbdir_summary_by_dataset_method.csv"
    )
    canonical = {}
    for row in canonical_summary:
        key = (row["dataset"].casefold(), row["method"].casefold())
        values = scores[key]
        expected = float(row["new_mean"])
        if len(values) != int(row["n_genes"]) or not math.isclose(
            statistics.mean(values), expected, rel_tol=0, abs_tol=1e-12
        ):
            raise ValueError(f"Reference-30 source mismatch: {key}")
        canonical[key] = (expected, len(values))

    dataset_path = RESULTS / "fig2_gene_metric_dataset_means.csv"
    dataset_fields, dataset_rows = read_csv(dataset_path)
    by_method = defaultdict(list)
    count = 0
    for row in dataset_rows:
        if row["Metric"] != "Real-CBDir":
            continue
        key = (dataset_names[row["Dataset"].casefold()], method_keys[row["Method"]])
        mean, n_genes = canonical[key]
        if int(float(row["N"])) != n_genes:
            raise ValueError(f"Gene count mismatch: {key}")
        row["Mean"] = repr(mean)
        by_method[row["Method"]].append((mean, n_genes))
        count += 1
    if count != 63 or len(by_method) != 9 or any(len(v) != 7 for v in by_method.values()):
        raise ValueError("Expected exactly 9 methods by 7 real datasets")

    summary_path = RESULTS / "fig2_gene_metric_summary.csv"
    summary_fields, summary_rows = read_csv(summary_path)
    for row in summary_rows:
        if row["Metric"] != "Real-CBDir":
            continue
        values = by_method[row["Method"]]
        if int(row["N_dataset_or_system"]) != len(values) or int(
            row["N_raw_genes_or_rows"]
        ) != sum(n for _, n in values):
            raise ValueError(f"Summary count mismatch: {row['Method']}")
        row["Mean"] = repr(statistics.mean(x for x, _ in values))
        row["Std"] = repr(statistics.stdev(x for x, _ in values))

    for row in wps_rows:
        if row["Figure2_Method"] not in by_method:
            continue
        values = by_method[row["Figure2_Method"]]
        mean = statistics.mean(x for x, _ in values)
        std = statistics.stdev(x for x, _ in values)
        row["Real_CBDir_mean_plus_1"] = f"{mean + 1:.4f}"
        row["Real_CBDir_SD"] = f"{std:.4f}"
        if int(float(row["Real_CBDir_n_datasets"])) != len(values) or int(
            float(row["Real_CBDir_n_gene_rows"])
        ) != sum(n for _, n in values):
            raise ValueError(f"Figure count mismatch: {row['Figure2_Method']}")

    updates = [(dataset_path, dataset_fields, dataset_rows),
               (summary_path, summary_fields, summary_rows),
               (wps_path, wps_fields, wps_rows)]
    for path, fields, rows in updates:
        output = io.StringIO(newline="")
        writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        if write:
            path.write_text(output.getvalue())
        elif path.read_text() != output.getvalue():
            raise SystemExit(f"Outdated Real-CBDir summary: {path}")
    print(f"Reference-30 Real-CBDir: {count} dataset means, 9 method summaries")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Update the three summary CSVs")
    main(parser.parse_args().write)
