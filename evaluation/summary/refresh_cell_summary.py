"""Derive the Figure 7 Cell statistics from the complete frozen Cell table."""

import argparse
import csv
import io
import math
from pathlib import Path
import statistics


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / "evaluation/cell/results/cell_metrics.csv"
DEFAULT_OUTPUT = ROOT / "evaluation/summary/results/figure7_summary/fig7_cell_metric_summary.csv"
METRICS = ("consistency", "cbdir", "iccoh", "transcosine", "batc")
FIELDS = ["config_id", "method", "manuscript_name", "metric", "n_datasets",
          "mean", "sample_sd", "median", "q1", "q3", "iqr",
          "q3_minus_mean", "mean_minus_q1", "mean_plus_1"]


def quantile(values, probability):
    """Linear interpolation, matching NumPy percentile and Excel QUARTILE.INC."""
    values = sorted(values)
    position = (len(values) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    return values[lower] + (values[upper] - values[lower]) * (position - lower)


def build_summary(source, repository_root=ROOT):
    with (repository_root / "methods/configurations.tsv").open(newline="") as handle:
        registry = {r["config_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    with (repository_root / "evaluation/cell/inputs.tsv").open(newline="") as handle:
        inputs = list(csv.DictReader(handle, delimiter="\t"))
    config_by_method = {}
    for row in inputs:
        method = row["operational_id"]
        if config_by_method.setdefault(method, row["config_id"]) != row["config_id"]:
            raise ValueError("Cell input configuration identity is inconsistent")
    with Path(source).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    identities = {(r["dataset"], r["method"]) for r in rows}
    expected = {(r["evaluation_dataset"], r["operational_id"]) for r in inputs}
    if len(rows) != 280 or len(identities) != 280 or identities != expected:
        raise ValueError("Cell summary requires the complete registered 14 by 20 table")
    result = []
    for method, config_id in sorted(config_by_method.items(), key=lambda item: item[1]):
        current = [r for r in rows if r["method"] == method]
        if len(current) != 14:
            raise ValueError("Each method must give every dataset equal weight")
        for metric in METRICS:
            values = [float(r[metric]) for r in current]
            if not all(math.isfinite(v) for v in values):
                raise ValueError("Nonfinite primary Cell metric in summary")
            mean = statistics.mean(values)
            q1, q3 = quantile(values, .25), quantile(values, .75)
            result.append(dict(config_id=config_id, method=method,
                manuscript_name=registry[config_id]["manuscript_name"], metric=metric,
                n_datasets=14, mean=mean, sample_sd=statistics.stdev(values),
                median=statistics.median(values), q1=q1, q3=q3, iqr=q3-q1,
                q3_minus_mean=q3-mean, mean_minus_q1=mean-q1, mean_plus_1=mean+1))
    return result


def refresh(source=DEFAULT_SOURCE, output=DEFAULT_OUTPUT, write=False):
    expected = build_summary(source)
    if write:
        content = io.StringIO(newline="")
        writer = csv.DictWriter(content, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(expected)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(content.getvalue())
    else:
        with output.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != FIELDS:
                raise ValueError("Cell summary schema differs")
            actual = list(reader)
        if len(actual) != len(expected):
            raise ValueError("Cell summary coverage differs")
        for old, new in zip(actual, expected):
            for key in FIELDS[:5]:
                if str(old[key]) != str(new[key]):
                    raise ValueError("Cell summary identity differs")
            for key in FIELDS[5:]:
                if not math.isclose(float(old[key]), new[key], rel_tol=0, abs_tol=1e-12):
                    raise ValueError("Outdated Cell summary: " + new["method"] + "/" + new["metric"] + "/" + key)
    print("Cell summary: 100 method/metric rows, 14 equal-weight datasets, sample SD")
    return expected


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    refresh(args.metrics, args.output, args.write)
