"""Reproduce Table S10A/B and descriptive BATC dataset-sensitivity results."""

import argparse
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from scipy.stats import f, friedmanchisquare


CELL = Path(__file__).resolve().parent
ROOT = CELL.parents[1]
sys.path.insert(0, str(ROOT))
from evaluation.gene.reproduce_crosslevel import matched_scores

MATCHED = CELL.parent / "gene/results/crosslevel/dataset_scores.csv"
CELL_METRICS = ("consistency", "iccoh", "cbdir", "transcosine", "batc")
MATCHED_METRICS = ("gene_cbdir", "cell_cbdir", "batc")


def validate_identities(cell, matched):
    """Check the complete Cell input map and the registered Gene-matched fit routes."""
    inputs = pd.read_csv(CELL / "inputs.tsv", sep="\t")
    expected_cell = inputs.rename(columns={"evaluation_dataset": "dataset", "operational_id": "method"})
    identity = ["dataset", "method", "config_id", "operational_id", "cell_dataset",
                "cell_score_route", "gene_result_path", "cell_result_path", "cell_result_sha256", "fit_relation"]
    for scores, expected, fields in (
        (cell, expected_cell, ["dataset", "method"]),
        (matched, pd.DataFrame(matched_scores(ROOT)), identity),
    ):
        observed = scores[fields]
        if observed.isna().any().any() or observed.duplicated(["dataset", "method"]).any():
            raise ValueError("Missing or duplicate dataset/method identity")
        if set(observed.itertuples(index=False, name=None)) != set(expected[fields].itertuples(index=False, name=None)):
            raise ValueError(f"Score identities differ from registered inputs for {fields}")


def omnibus(scores, metrics, shape):
    """Block by dataset, then correct Iman–Davenport p-values within the panel."""
    rows = []
    for metric in metrics:
        matrix = scores.pivot(index="dataset", columns="method", values=metric)
        if matrix.shape != shape or not np.isfinite(matrix.to_numpy()).all():
            raise ValueError(f"{metric} requires a complete, finite {shape} score matrix")
        n, k = shape
        max_q = n * (k - 1)
        ranks = matrix.rank(axis=1).to_numpy()
        if (ranks == ranks[:, :1]).all():
            q, statistic = 0.0, 0.0
        elif (ranks == ranks[0]).all():
            q, statistic = float(max_q), np.inf
        else:
            q = friedmanchisquare(*matrix.to_numpy().T).statistic
            statistic = (n - 1) * q / (max_q - q)
        df1, df2 = k - 1, (k - 1) * (n - 1)
        rows.append(dict(metric=metric, n_datasets=n, n_configurations=k,
                         friedman_q=q, kendall_w=q / max_q, iman_davenport_f=statistic,
                         df_numerator=df1, df_denominator=df2, p_raw=f.sf(statistic, df1, df2)))
    table = pd.DataFrame(rows)
    ordered = table.p_raw.sort_values()
    adjusted = (ordered * np.arange(len(metrics), 0, -1)).cummax().clip(upper=1)
    table["p_holm"] = adjusted
    return table


def batc_ranks(scores):
    """Rank all configurations by dataset and by equal-dataset mean after each omission."""
    batc = scores.pivot(index="dataset", columns="method", values="batc")
    per_dataset = batc.stack().rename("batc").to_frame()
    per_dataset["rank"] = batc.rank(axis=1, ascending=False, method="average").stack()
    omissions = []
    for dataset in batc.index:
        means = batc.drop(index=dataset).mean().rename("mean_batc").to_frame()
        means["rank"] = means.mean_batc.rank(ascending=False, method="average")
        means["omitted_dataset"] = dataset
        omissions.append(means.reset_index())
    return per_dataset.reset_index(), pd.concat(omissions, ignore_index=True)[
        ["omitted_dataset", "method", "mean_batc", "rank"]
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell", type=Path, default=CELL / "results/cell_metrics.csv")
    parser.add_argument("--matched", type=Path, default=MATCHED)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        parser.error("--output requires a new directory")
    cell = pd.read_csv(args.cell, float_precision="round_trip")
    matched = pd.read_csv(args.matched, float_precision="round_trip")
    validate_identities(cell, matched)
    primary_tests = omnibus(cell, CELL_METRICS, (14, 20))
    matched_tests = omnibus(matched, MATCHED_METRICS, (7, 9))
    ranks, omissions = batc_ranks(cell)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.output.parent) as temporary:
        bundle = Path(temporary) / "statistics"
        bundle.mkdir()
        for name, table in (("s10a", primary_tests), ("s10b", matched_tests),
                            ("batc_dataset_ranks", ranks), ("batc_leave_one_out", omissions)):
            table.to_csv(bundle / f"{name}.csv", index=False)
        if args.output.exists() or args.output.is_symlink():
            raise FileExistsError(args.output)
        bundle.rename(args.output)


if __name__ == "__main__":
    main()
