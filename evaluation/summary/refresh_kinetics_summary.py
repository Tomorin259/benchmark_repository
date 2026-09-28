"""Refresh integrated kinetic summaries from the current full-gene scores."""

import argparse
import csv
import hashlib
import io
import math
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = "evaluation/kinetics_recovery/results/central_results"
RESULTS = ROOT / "evaluation/summary/results/figure6_summary"


def read_csv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def main(write):
    _, genes = read_csv(ROOT / SOURCE / "gene_metrics.csv")
    scores = defaultdict(list)
    seen = set()
    for row in genes:
        key = (row["Method"], row["Dataset"], row["Branch"], row["Parameter"])
        identity = (*key, row["Gene"])
        if identity in seen:
            raise ValueError(f"Duplicate gene score: {identity}")
        seen.add(identity)
        value = float(row["Spearman_Corr"] or "nan")
        if math.isfinite(value):
            scores[key].append(value)

    # Validate the same aggregation boundaries used by the evaluator.
    _, components = read_csv(ROOT / SOURCE / "component_means.csv")
    systems = defaultdict(list)
    component_rows = []
    component_keys = set()
    for row in components:
        key = (row["Method"], row["Dataset"], row["Branch"], row["Parameter"])
        values = scores[key]
        mean = statistics.mean(values)
        if key in component_keys or len(values) != int(row["N"]) or not math.isclose(
            mean, float(row["Mean"]), rel_tol=0, abs_tol=1e-12
        ):
            raise ValueError(f"Canonical component mismatch: {key}")
        component_keys.add(key)
        systems[(row["Method"], row["Parameter"], row["Dataset"])].append(
            (mean, len(values))
        )
        component_rows.append(dict(Method=row["Method"], Component=row["Component"],
                                   Mean=repr(mean), N=str(len(values)),
                                   Path=f"{SOURCE}/gene_metrics.csv"))
    if component_keys != set(scores):
        raise ValueError("Gene/component coverage differs")

    _, canonical_systems = read_csv(ROOT / SOURCE / "system_means.csv")
    methods = defaultdict(list)
    system_rows = []
    system_keys = set()
    for row in canonical_systems:
        key = (row["Method"], row["Parameter"], row["System"])
        values = systems[key]
        mean = statistics.mean(value for value, _ in values)
        count = sum(n for _, n in values)
        if key in system_keys or len(values) != int(row["N_components"]) or count != int(
            row["N_genes"]
        ) or not math.isclose(mean, float(row["Mean"]), rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"Canonical system mismatch: {key}")
        system_keys.add(key)
        methods[key[:2]].append((mean, count))
        system_rows.append(dict(Metric=f"Kinetic-Spearman-{row['Parameter']}",
                               Method=row["Method"], Dataset=row["System"],
                               Mean=repr(mean), N=str(len(values)), Unit="components",
                               N_raw_genes_or_rows=str(count),
                               Aggregation="within_parameter_components"))
    if system_keys != set(systems):
        raise ValueError("Component/system coverage differs")

    _, canonical_methods = read_csv(ROOT / SOURCE / "method_summary.csv")
    method_rows = []
    method_keys = set()
    for row in canonical_methods:
        key = (row["Method"], row["Parameter"])
        values = methods[key]
        mean = statistics.mean(value for value, _ in values)
        count = sum(n for _, n in values)
        if key in method_keys or len(values) != int(row["N_systems"]) or count != int(
            row["N_genes"]
        ) or not math.isclose(mean, float(row["Mean"]), rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"Canonical method mismatch: {key}")
        method_keys.add(key)
        std = statistics.stdev(value for value, _ in values) if len(values) > 1 else None
        method_rows.append(dict(Metric=f"Kinetic-Spearman-{row['Parameter']}",
                               Method=row["Method"], Mean=repr(mean),
                               Std="" if std is None else repr(std),
                               N_dataset_or_system=str(len(values)),
                               N_raw_genes_or_rows=str(count)))
    if method_keys != set(methods):
        raise ValueError("System/method coverage differs")

    updates = []
    for name, replacement in [("kinetic_component_means", component_rows),
                              ("gene_metric_dataset_means", system_rows),
                              ("gene_metric_summary", method_rows)]:
        path = RESULTS / f"fig2_{name}.csv"
        fields, rows = read_csv(path)
        retained = [] if name == "kinetic_component_means" else [
            row for row in rows if not row["Metric"].startswith("Kinetic-Spearman-")
        ]
        updates.append((path, fields, retained + replacement))

    path = RESULTS / "fig2_gene_metric_wps_values.csv"
    fields, rows = read_csv(path)
    mapped = set()
    for row in rows:
        for parameter in ("alpha", "gamma"):
            key = (row["Figure2_Method"], parameter)
            prefix = f"Kinetic_Spearman_{parameter}"
            if key in methods:
                values = methods[key]
                row[f"{prefix}_mean_plus_1"] = repr(statistics.mean(x for x, _ in values) + 1)
                row[f"{prefix}_n_systems"] = str(len(values))
                row[f"{prefix}_n_gene_rows"] = str(sum(n for _, n in values))
                mapped.add(key)
            else:
                for suffix in ("mean_plus_1", "n_systems", "n_gene_rows"):
                    row[f"{prefix}_{suffix}"] = ""
    if mapped != set(methods):
        raise ValueError("WPS method coverage differs from canonical results")
    updates.append((path, fields, rows))

    rendered = {}
    for path, fields, rows in updates:
        output = io.StringIO(newline="")
        writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        rendered[path] = output.getvalue()

    registry_path = ROOT / "artifacts.tsv"
    with registry_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        registry_fields, registry_rows = reader.fieldnames, list(reader)
    hashes = {str(path.relative_to(ROOT)): hashlib.sha256(text.encode()).hexdigest()
              for path, text in rendered.items()}
    registered = []
    for row in registry_rows:
        if row["path"] in hashes:
            row["sha256"] = hashes[row["path"]]
            registered.append(row["path"])
    if len(registered) != len(hashes) or set(registered) != set(hashes):
        raise ValueError("Summary artifact registration is missing or duplicated")
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=registry_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(registry_rows)
    rendered[registry_path] = output.getvalue()

    for path, text in rendered.items():
        if write:
            path.write_text(text)
        elif path.read_text() != text:
            raise SystemExit(f"Outdated kinetics summary: {path}")
    print(f"Current kinetics: {len(component_rows)} components, "
          f"{len(system_rows)} systems, {len(method_rows)} method/parameter summaries")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Update the four summary CSVs")
    main(parser.parse_args().write)
