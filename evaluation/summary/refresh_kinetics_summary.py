"""Refresh integrated kinetic summaries from the expanded-reference coverage table."""

import argparse
import csv
import hashlib
import io
import math
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = "evaluation/kinetics_recovery/results/expanded_reference/coverage_summary.csv"
RESULTS = ROOT / "evaluation/summary/results/figure7_summary"


def read_csv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def csv_values_equivalent(current_text, refreshed_text, tolerance=1e-15):
    """Compare CSV content while ignoring harmless float serialization drift."""
    current = csv.DictReader(io.StringIO(current_text))
    refreshed = csv.DictReader(io.StringIO(refreshed_text))
    if current.fieldnames != refreshed.fieldnames:
        return False
    current_rows, refreshed_rows = list(current), list(refreshed)
    if len(current_rows) != len(refreshed_rows):
        return False
    for left, right in zip(current_rows, refreshed_rows):
        for field in current.fieldnames:
            if left[field] == right[field]:
                continue
            try:
                if math.isclose(float(left[field]), float(right[field]),
                                rel_tol=0.0, abs_tol=tolerance):
                    continue
            except ValueError:
                pass
            return False
    return True


def main(write):
    _, coverage = read_csv(ROOT / SOURCE)
    systems = defaultdict(dict)
    component_rows = []
    for row in sorted(coverage, key=lambda r: (r["Parameter"], r["Method"], r["Panel"])):
        method, parameter, system = row["Method"], row["Parameter"], row["Dataset"]
        if system == "rpe1" and not row["Branch"]:
            component = f"rpe1_{parameter}"
            expected_panel = f"rpe1_all_{parameter}"
        elif system == "organoid" and row["Branch"] in ("1.0", "2.0", "1", "2"):
            component = f"organoid_branch{int(float(row['Branch']))}_{parameter}"
            expected_panel = component
        else:
            raise ValueError(f"Unknown system or branch in {row['Panel']}")
        group = systems[(parameter, method, system)]
        count = int(row["Valid_N"])
        mean = float(row["Mean"])
        if (parameter not in ("alpha", "gamma") or row["Panel"] != expected_panel
                or component in group or not math.isfinite(mean)
                or not 0 < count <= int(row["Candidate_N"])):
            raise ValueError(f"Invalid coverage for {method}, {row['Panel']}")
        group[component] = (mean, count)
        component_rows.append(dict(Method=method, Component=component,
                                   Mean=row["Mean"], N=str(count), Path=SOURCE))

    # Each branch has equal weight, irrespective of its number of finite gene scores.
    methods = defaultdict(dict)
    system_rows = []
    for (parameter, method, system), components in sorted(systems.items()):
        expected = ({f"rpe1_{parameter}"} if system == "rpe1" else
                    {f"organoid_branch1_{parameter}", f"organoid_branch2_{parameter}"})
        if set(components) != expected:
            raise ValueError(f"Incomplete {system} components for {method}, {parameter}")
        mean = statistics.mean(value for value, _ in components.values())
        count = sum(n for _, n in components.values())
        methods[(method, parameter)][system] = (mean, count)
        system_rows.append(dict(Metric=f"Kinetic-Spearman-{parameter}", Method=method,
                                Dataset=system, Mean=repr(mean), N=str(len(components)),
                                Unit="components", N_raw_genes_or_rows=str(count),
                                Aggregation="within_parameter_components"))

    method_rows = []
    for (method, parameter), system_values in methods.items():
        if set(system_values) != {"rpe1", "organoid"}:
            raise ValueError(f"Missing biological system for {method}, {parameter}")
        values = [mean for mean, _ in system_values.values()]
        count = sum(n for _, n in system_values.values())
        # Round the sample variance to a float before taking its square root on every Python version.
        sd = math.sqrt(statistics.variance(values))
        method_rows.append(dict(Metric=f"Kinetic-Spearman-{parameter}", Method=method,
                                Mean=repr(statistics.mean(values)), Std=repr(sd),
                                N_dataset_or_system=str(len(values)), N_raw_genes_or_rows=str(count)))

    updates = []
    for name, replacement in [("kinetic_component_means", component_rows),
                              ("gene_metric_dataset_means", system_rows),
                              ("gene_metric_summary", method_rows)]:
        path = RESULTS / f"fig7_{name}.csv"
        fields, rows = read_csv(path)
        retained = [] if name == "kinetic_component_means" else [
            row for row in rows if not row["Metric"].startswith("Kinetic-Spearman-")
        ]
        updates.append((path, fields, retained + replacement))

    path = RESULTS / "fig7_gene_metric_wps_values.csv"
    fields, rows = read_csv(path)
    mapped = set()
    for row in rows:
        for parameter in ("alpha", "gamma"):
            key = (row["Figure7_Method"], parameter)
            prefix = f"Kinetic_Spearman_{parameter}"
            if key in methods:
                values = methods[key].values()
                row[f"{prefix}_mean_plus_1"] = repr(statistics.mean(x for x, _ in values) + 1)
                row[f"{prefix}_n_systems"] = str(len(values))
                row[f"{prefix}_n_gene_rows"] = str(sum(n for _, n in values))
                mapped.add(key)
            else:
                for suffix in ("mean_plus_1", "n_systems", "n_gene_rows"):
                    row[f"{prefix}_{suffix}"] = ""
    if mapped != set(methods):
        raise ValueError("WPS method coverage differs from expanded-reference results")
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
    equivalent_csvs = {
        path for path, text in rendered.items()
        if path.suffix == ".csv" and path.exists()
        and csv_values_equivalent(path.read_text(), text)
    }
    hashes = {
        str(path.relative_to(ROOT)): hashlib.sha256(text.encode()).hexdigest()
        for path, text in rendered.items()
        if path not in equivalent_csvs
    }
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
        if path.read_text() == text or path in equivalent_csvs:
            continue
        if write:
            path.write_text(text)
        else:
            raise SystemExit(f"Outdated kinetics summary: {path}")
    print(f"Current kinetics: {len(component_rows)} components, "
          f"{len(system_rows)} systems, {len(method_rows)} method/parameter summaries")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Update the four summary CSVs")
    main(parser.parse_args().write)
