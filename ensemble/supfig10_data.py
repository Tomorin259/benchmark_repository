#!/usr/bin/env python3
"""Recompute the three Supplementary Figure 10 source tables."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ensemble.core import DEFAULT_CONFIG, load_config, sha256
from ensemble.pareto_bootstrap import (
    METRICS,
    ROOT,
    SOURCE_NAMES,
    bootstrap_rankings,
    dominance,
    load_formal_scores,
)

ENSEMBLE = "GeoMedian-5"
CONSTITUENT_CONFIGS = {"C05", "C08", "C12", "C16", "C19"}
N_RESAMPLES = 10_000
BOOTSTRAP_SEED = 20261109


def load_full5(path: Path, datasets: tuple[str, ...]) -> pd.DataFrame:
    """Read the S5B dataset table, the S5A detail, or the normalization scores."""
    selected = pd.read_csv(path)
    columns = ["dataset", *METRICS]
    if missing := set(columns) - set(selected.columns):
        raise ValueError(f"{path} is missing dataset score columns: {sorted(missing)}")

    if "id" in selected:
        selected = selected.loc[selected.id.eq("full5"), columns]
    elif {"condition", "n_members"} <= set(selected.columns):
        selected = selected.loc[selected.condition.eq("normalized") & selected.n_members.eq(5), columns]
    elif "variant" in selected:
        if not selected.variant.eq("normalized_geomedian5").all():
            raise ValueError(f"{path} must identify normalized full5 scores")
        selected = selected.loc[:, columns]
    elif set(selected.columns) != {"dataset", "display_name", *METRICS}:
        raise ValueError(f"{path} must contain Table S5B or dataset scores identifying normalized full5")

    summaries = selected.loc[selected.dataset.eq("macro_average"), list(METRICS)]
    selected = selected.loc[~selected.dataset.eq("macro_average"), columns]
    if len(selected) != len(datasets) or set(selected.dataset) != set(datasets):
        raise ValueError(f"{path} must contain one full5 score pair for each of the 14 datasets")
    selected = selected.set_index("dataset").loc[list(datasets)].reset_index()
    selected[list(METRICS)] = selected[list(METRICS)].apply(pd.to_numeric)
    if not np.isfinite(selected[list(METRICS)]).all().all():
        raise ValueError(f"{path} contains nonfinite full5 scores")
    if len(summaries) > 1 or (len(summaries) == 1 and not np.allclose(
        summaries.to_numpy(dtype=float)[0], selected[list(METRICS)].mean(), rtol=0, atol=1e-12
    )):
        raise ValueError(f"{path} macro_average does not match its dataset scores")
    return selected.assign(method=ENSEMBLE)


def write_panel_data(selected_score_path: Path, output_dir: Path) -> None:
    """Write the three source CSVs in a new directory given dataset-level full5 scores."""
    selected_score_path = Path(selected_score_path)
    output_dir = Path(output_dir).resolve()
    if output_dir.is_relative_to(ROOT / "figures/assembly"):
        raise ValueError("write panel data to a new directory outside the figure assembly")
    if output_dir.exists():
        raise FileExistsError(f"panel output directory already exists: {output_dir}")

    sources = [*(ROOT / name for name in SOURCE_NAMES), DEFAULT_CONFIG, selected_score_path]
    hashes = {path: sha256(path) for path in sources}
    rows, configs = load_formal_scores()
    ranking, _, _, evaluation_datasets, _ = bootstrap_rankings(rows, configs, N_RESAMPLES, BOOTSTRAP_SEED)
    panel_a = pd.DataFrame(ranking)

    # The Cell benchmark uses the stems of the frozen neighbor files as dataset IDs.
    # S5 uses the ensemble config keys, which differ for the La Manno dentate gyrus.
    dataset_keys = {Path(spec["reference_neighbor_file"]).stem: key
                    for key, spec in load_config()["datasets"].items()}
    if set(dataset_keys) != set(evaluation_datasets):
        raise ValueError("the ensemble and Cell benchmark use different dataset grids")
    selected = load_full5(selected_score_path, tuple(dataset_keys.values()))
    native = pd.read_csv(ROOT / SOURCE_NAMES[0], usecols=["dataset", "method", *METRICS])
    native["dataset"] = native.dataset.map(dataset_keys)

    labels = {config["method"]: config["manuscript_name"] for config in configs}
    constituents = {config["method"] for config in configs if config["config_id"] in CONSTITUENT_CONFIGS}
    labels[ENSEMBLE] = ENSEMBLE
    wide_scores = pd.concat([native, selected], ignore_index=True)
    wide_scores["method_label"] = wide_scores.method.map(labels)
    means = wide_scores.groupby("method", sort=False)[list(METRICS)].mean()
    orders = {metric: means[metric].sort_values(ascending=False, kind="stable").index for metric in METRICS}

    long_scores = wide_scores.melt(
        id_vars=["dataset", "method", "method_label"], value_vars=METRICS,
        var_name="metric", value_name="score",
    )
    long_scores["constituent"] = long_scores.method.isin(constituents)
    plot_rows = {(metric, method): index for metric, order in orders.items() for index, method in enumerate(order)}
    long_scores["plot_row"] = [plot_rows[metric, method] for metric, method in zip(
        long_scores.metric, long_scores.method
    )]

    macro = means.loc[orders["cbdir"]].reset_index()
    macro["method_label"] = macro.method.map(labels)
    macro["single_method_front"] = macro.method.isin(panel_a.loc[panel_a.original_mean_front, "method"])
    macro["overall_front"] = ~dominance(macro[list(METRICS)].to_numpy()).any(axis=0)

    if any(sha256(path) != digest for path, digest in hashes.items()):
        raise ValueError("a score or configuration source changed during S10 generation")
    output_dir.mkdir(parents=True)
    panel_a.to_csv(output_dir / "panel_a_data.csv", index=False)
    long_scores.to_csv(output_dir / "plot_data.csv", index=False, lineterminator="\r\n")
    macro.to_csv(output_dir / "macro_means.csv", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-scores", type=Path, required=True,
                        help="Table_S5_B.csv or dataset scores identifying normalized full5")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    write_panel_data(args.selected_scores, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
