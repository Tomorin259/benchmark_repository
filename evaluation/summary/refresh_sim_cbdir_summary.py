"""Refresh Figure 7 simulated CBDir-gene from the frozen reference-30 scores."""

import argparse
import csv
import hashlib
import io
import math
import statistics
from collections import defaultdict
from pathlib import Path

try:
    from .refresh_real_cbdir_summary import check_summary, read_csv
except ImportError:
    from refresh_real_cbdir_summary import check_summary, read_csv


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "evaluation/summary/results/figure7_summary"
SCENARIOS = {"Normal + MURK": ("murk", 10),
             "Normal + branching": ("branching", 20),
             "Hybrid (normal + MURK + branching)": ("all", 40)}


def main(write):
    wps_path = RESULTS / "fig7_gene_metric_wps_values.csv"
    wps_fields, wps_rows = read_csv(wps_path)
    method_keys = {row["Figure7_Method"]: row["Source_Method_Key"].casefold()
                   for row in wps_rows if row["Sim_CBDir_mean_plus_1"]}
    _, source = read_csv(ROOT / "evaluation/gene/results/gene_cbdir_reference30_all.csv")
    scores, seen = defaultdict(list), set()
    for row in source:
        if row["dataset"] not in {s for s, _ in SCENARIOS.values()}:
            continue
        key = (row["dataset"], row["method"].casefold())
        identity = (*key, row["gene"])
        value = float(row["new_cbdir"] or "nan")
        if identity in seen or row["status"] != "ok" or not math.isfinite(value):
            raise ValueError(f"Invalid or duplicate simulation score: {identity}")
        seen.add(identity)
        scores[key].append(value)
    if len(method_keys) != 9 or len(seen) != 630:
        raise ValueError("Expected nine methods and 630 frozen special-gene scores")

    dataset_path = RESULTS / "fig7_gene_metric_dataset_means.csv"
    dataset_fields, dataset_rows = read_csv(dataset_path)
    by_method, matched = defaultdict(list), set()
    for row in dataset_rows:
        if row["Metric"] != "Sim-CBDir":
            continue
        scenario, expected_n = SCENARIOS[row["Dataset"]]
        key = (scenario, method_keys[row["Method"]])
        values = scores[key]
        if key in matched or len(values) != expected_n or int(float(row["N"])) != expected_n:
            raise ValueError(f"Simulation gene coverage differs: {key}")
        matched.add(key)
        row["Mean"] = repr(statistics.mean(values))
        by_method[row["Method"]].append((float(row["Mean"]), expected_n))
    if matched != set(scores):
        raise ValueError("Figure 7 scenario/method coverage differs from reference-30")

    summary_path = RESULTS / "fig7_gene_metric_summary.csv"
    summary_fields, summary_rows = read_csv(summary_path)
    for row in summary_rows:
        if row["Metric"] != "Sim-CBDir":
            continue
        values = by_method[row["Method"]]
        if len(values) != 3 or int(row["N_dataset_or_system"]) != 3 or int(row["N_raw_genes_or_rows"]) != 70:
            raise ValueError(f"Simulation summary coverage differs: {row['Method']}")
        row["Mean"] = repr(statistics.mean(v for v, _ in values))
        row["Std"] = repr(statistics.stdev(v for v, _ in values))
    for row in wps_rows:
        if row["Figure7_Method"] in by_method:
            row["Sim_CBDir_mean_plus_1"] = f"{statistics.mean(v for v, _ in by_method[row['Figure7_Method']]) + 1:.4f}"

    rendered = {}
    for path, fields, rows in [(dataset_path, dataset_fields, dataset_rows),
                               (summary_path, summary_fields, summary_rows),
                               (wps_path, wps_fields, wps_rows)]:
        output = io.StringIO(newline="")
        writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        rendered[path] = output.getvalue()
        if not write:
            check_summary(path, fields, rows, metric="Sim-CBDir")
    with (ROOT / "artifacts.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, registry = reader.fieldnames, list(reader)
    hashes = {str(path.relative_to(ROOT)): hashlib.sha256(text.encode()).hexdigest()
              for path, text in rendered.items()}
    registered = set()
    for row in registry:
        if row["path"] in hashes:
            if not write and row["sha256"] != hashlib.sha256((ROOT / row["path"]).read_bytes()).hexdigest():
                raise SystemExit(f"Outdated Figure 7 artifact hash: {row['path']}")
            row["sha256"] = hashes[row["path"]]
            registered.add(row["path"])
    if registered != set(hashes):
        raise ValueError("Figure 7 summary artifact registration is incomplete")
    if write:
        for path, text in rendered.items():
            path.write_text(text)
        with (ROOT / "artifacts.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(registry)
    print("Reference-30 Sim-CBDir: 27 scenario means, 9 equal-scenario means; CosSim unchanged")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    main(parser.parse_args().write)
