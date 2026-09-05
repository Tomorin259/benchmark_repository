#!/usr/bin/env python3
"""Render kinetic-parameter recovery figures from one metrics output bundle.

The plotting layer consumes only the CSV files written by metrics.py. It does
not reopen method-result H5AD files or recompute benchmark metrics.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MODULE_DIR / "config.json"

METHOD_ORDER = ["Dynamo", "DeepVelo", "cellDancer", "TSvelo", "DeepKINET"]
COLORS = {
    "Dynamo": "#2ca02c",
    "DeepVelo": "#1f77b4",
    "cellDancer": "#ff7f0e",
    "TSvelo": "#e420ca",
    "DeepKINET": "#B5BD61",
}


def ensure_output(path: Path) -> None:
    if path.exists() and not path.is_dir():
        raise NotADirectoryError(path)
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Output directory must be empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def read_inputs(metrics_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_path = metrics_dir / "gene_metrics.csv"
    summary_path = metrics_dir / "method_summary.csv"
    required = [path for path in (gene_path, summary_path) if not path.is_file()]
    if required:
        raise FileNotFoundError(
            "Expected gene_metrics.csv and method_summary.csv in "
            f"{metrics_dir}; run metrics.py first."
        )
    genes = pd.read_csv(gene_path)
    summary = pd.read_csv(summary_path)
    required_genes = {"Method", "Dataset", "Parameter", "Spearman_Corr"}
    required_summary = {"Method", "Parameter", "Mean"}
    if missing := required_genes - set(genes.columns):
        raise ValueError(f"gene_metrics.csv is missing {sorted(missing)}")
    if missing := required_summary - set(summary.columns):
        raise ValueError(f"method_summary.csv is missing {sorted(missing)}")
    return genes, summary


def read_matrix(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required matrix does not exist: {path}")
    return pd.read_csv(path, index_col=0).apply(pd.to_numeric, errors="coerce")


def ordered_methods(values) -> list[str]:
    present = list(pd.Index(values).dropna().unique())
    return [method for method in METHOD_ORDER if method in present] + sorted(
        set(present) - set(METHOD_ORDER)
    )


def zscore(values) -> np.ndarray:
    series = pd.Series(values, dtype=float)
    std = series.std()
    if pd.isna(std) or std == 0:
        return (series - series.mean()).to_numpy()
    return ((series - series.mean()) / std).to_numpy()


def plot_rpe1_trajectories(
    metrics_dir: Path,
    gene_metrics: pd.DataFrame,
    figure_genes: list[str],
    parameter: str,
    output: Path,
) -> pd.DataFrame:
    """Plot configured RPE1 genes with display-only per-curve z-scores."""
    truth = read_matrix(
        metrics_dir / "matrices" / "rpe1" / f"reference_{parameter}.csv"
    )
    methods = ordered_methods(
        gene_metrics.loc[
            gene_metrics["Dataset"].eq("rpe1")
            & gene_metrics["Parameter"].eq(parameter),
            "Method",
        ]
    )
    method_frames = {}
    for method in methods:
        path = metrics_dir / "trajectories" / "rpe1" / f"{method}_{parameter}.csv"
        if path.is_file():
            method_frames[method] = read_matrix(path)
    methods = ordered_methods(method_frames)
    genes = [gene for gene in figure_genes if gene in truth.index]
    if not genes:
        raise ValueError(f"No configured RPE1 trajectory genes occur in {parameter}")

    ncols = 5
    nrows = math.ceil((len(genes) + 1) / ncols)
    figure, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 2.6 * nrows))
    axes = np.asarray(axes).reshape(-1)
    tidy_rows = []
    for index, axis in enumerate(axes):
        if index >= len(genes):
            axis.axis("off")
            continue
        gene = genes[index]
        truth_values = zscore(truth.loc[gene])
        truth_x = np.linspace(0, 1, len(truth_values))
        axis.plot(
            truth_x,
            truth_values,
            color="black",
            linewidth=2.2,
            zorder=10,
            label="scEU-seq",
        )
        tidy_rows.extend(
            {
                "Gene": gene,
                "Source": "scEU-seq",
                "Relative_Position": position,
                "Display_Z": value,
            }
            for position, value in zip(truth_x, truth_values)
        )
        for method in methods:
            frame = method_frames[method]
            if gene not in frame.index:
                continue
            values = zscore(frame.loc[gene])
            positions = np.linspace(0, 1, len(values))
            axis.plot(
                positions,
                values,
                color=COLORS.get(method, "0.4"),
                linewidth=1.4,
                alpha=0.85,
                label=method,
            )
            tidy_rows.extend(
                {
                    "Gene": gene,
                    "Source": method,
                    "Relative_Position": position,
                    "Display_Z": value,
                }
                for position, value in zip(positions, values)
            )
        axis.set_title(gene, fontsize=10, fontweight="bold")
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)

    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="lower center",
        ncol=min(len(labels), 6),
        frameon=False,
        bbox_to_anchor=(0.5, -0.01),
    )
    figure.tight_layout(rect=(0, 0.06, 1, 1))
    figure.savefig(
        output / f"rpe1_{parameter}_trajectories.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)
    return pd.DataFrame(tidy_rows)


def heatmap_matrix(
    gene_metrics: pd.DataFrame, parameter: str, genes: list[str]
) -> pd.DataFrame:
    subset = gene_metrics[
        gene_metrics["Dataset"].eq("rpe1")
        & gene_metrics["Parameter"].eq(parameter)
        & gene_metrics["Gene"].isin(genes)
    ]
    matrix = subset.pivot_table(
        index="Method", columns="Gene", values="Spearman_Corr", aggfunc="first"
    )
    columns = [gene for gene in genes if gene in matrix.columns]
    matrix = matrix.reindex(columns=columns)
    methods = ordered_methods(matrix.index)
    return matrix.reindex(index=methods)


def plot_rpe1_heatmap(
    gene_metrics: pd.DataFrame, figure_genes: list[str], output: Path
) -> dict[str, pd.DataFrame]:
    matrices = {
        parameter: heatmap_matrix(gene_metrics, parameter, figure_genes)
        for parameter in ("alpha", "gamma")
    }
    figure, axes = plt.subplots(1, 2, figsize=(24, 5.2))
    image = None
    for axis, parameter in zip(axes, ("alpha", "gamma")):
        matrix = matrices[parameter]
        values = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("#eeeeee")
        image = axis.imshow(values, cmap=cmap, vmin=-0.5, vmax=1, aspect="auto")
        for row, method in enumerate(matrix.index):
            for column, gene in enumerate(matrix.columns):
                value = matrix.loc[method, gene]
                if np.isfinite(value):
                    axis.text(
                        column,
                        row,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        fontsize=7,
                    )
        axis.set_title(f"RPE1 {parameter}", fontweight="bold")
        axis.set_yticks(range(len(matrix.index)), matrix.index)
        axis.set_xticks(
            range(len(matrix.columns)), matrix.columns, rotation=45, ha="right"
        )
        axis.set_xlabel("Gene")
        axis.set_ylabel("Method")
    figure.subplots_adjust(left=0.04, right=0.89, bottom=0.27, top=0.88, wspace=0.15)
    colorbar_axis = figure.add_axes([0.92, 0.27, 0.012, 0.55])
    figure.colorbar(image, cax=colorbar_axis, label="Spearman correlation")
    figure.savefig(
        output / "rpe1_alpha_gamma_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)
    return matrices


def cosine_classes(alpha: pd.DataFrame, gamma: pd.DataFrame) -> dict[str, set[str]]:
    genes = alpha.index.intersection(gamma.index, sort=False)
    columns = alpha.columns.intersection(gamma.columns, sort=False)
    classes = {"Cooperative": set(), "Neutral": set(), "Destabilizing": set()}
    for gene in genes:
        left = alpha.loc[gene, columns].to_numpy(dtype=float)
        right = gamma.loc[gene, columns].to_numpy(dtype=float)
        valid = np.isfinite(left) & np.isfinite(right)
        if not valid.any():
            continue
        denominator = np.linalg.norm(left[valid]) * np.linalg.norm(right[valid])
        if denominator == 0:
            continue
        similarity = float(np.dot(left[valid], right[valid]) / denominator)
        if similarity <= -0.5:
            classes["Cooperative"].add(gene)
        elif similarity >= 0.5:
            classes["Destabilizing"].add(gene)
        else:
            classes["Neutral"].add(gene)
    return classes


def draw_set_overlap(
    axis: plt.Axes,
    title: str,
    reference: set[str],
    prediction: set[str],
    prediction_label: str,
) -> dict[str, int]:
    """Draw a dependency-free two-set overlap panel."""
    reference_only = len(reference - prediction)
    prediction_only = len(prediction - reference)
    intersection = len(reference & prediction)
    left = Circle((0.42, 0.5), 0.28, facecolor="#80b1d3", alpha=0.6, edgecolor="0.4")
    right = Circle((0.58, 0.5), 0.28, facecolor="#fb8072", alpha=0.6, edgecolor="0.4")
    axis.add_patch(left)
    axis.add_patch(right)
    axis.text(0.30, 0.5, str(reference_only), ha="center", va="center", fontsize=14)
    axis.text(0.70, 0.5, str(prediction_only), ha="center", va="center", fontsize=14)
    axis.text(0.50, 0.5, str(intersection), ha="center", va="center", fontsize=14)
    axis.text(0.31, 0.17, "scEU-seq", ha="center", fontsize=9)
    axis.text(0.69, 0.17, prediction_label, ha="center", fontsize=9)
    axis.set_title(title, fontsize=10, fontweight="bold")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    return {
        "scEU_only": reference_only,
        f"{prediction_label}_only": prediction_only,
        "intersection": intersection,
    }


def plot_kinetic_classes(
    method: str,
    metrics_dir: Path,
    output: Path,
) -> pd.DataFrame:
    """Compare reference and method kinetic classes for one method."""
    truth = {
        parameter: read_matrix(
            metrics_dir / "matrices" / "rpe1" / f"reference_{parameter}.csv"
        )
        for parameter in ("alpha", "gamma")
    }
    paths = {
        parameter: metrics_dir
        / "trajectories"
        / "rpe1"
        / f"{method}_{parameter}.csv"
        for parameter in ("alpha", "gamma")
    }
    if any(not path.is_file() for path in paths.values()):
        return pd.DataFrame(
            columns=["Class", "scEU_only", f"{method}_only", "intersection"]
        )

    prediction = {parameter: read_matrix(path) for parameter, path in paths.items()}
    truth_classes = cosine_classes(truth["alpha"], truth["gamma"])
    predicted_classes = cosine_classes(prediction["alpha"], prediction["gamma"])

    figure, axes = plt.subplots(1, 3, figsize=(10, 3.5))
    rows = []
    for axis, name in zip(axes, ("Cooperative", "Neutral", "Destabilizing")):
        counts = draw_set_overlap(
            axis,
            name,
            truth_classes[name],
            predicted_classes[name],
            method,
        )
        rows.append({"Class": name, **counts})
    figure.tight_layout()
    figure.savefig(
        output / f"rpe1_{method.lower()}_kinetic_classes.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)
    return pd.DataFrame(rows)


def draw_boxplot(
    axis: plt.Axes, frame: pd.DataFrame, parameter: str
) -> None:
    subset = frame[frame["Parameter"].eq(parameter)].copy()
    methods = ordered_methods(subset["Method"])
    values = [
        subset.loc[subset["Method"].eq(method), "Spearman_Corr"].dropna().to_numpy()
        for method in methods
    ]
    if not methods:
        axis.set_visible(False)
        return
    box = axis.boxplot(
        values, labels=methods, patch_artist=True, showfliers=False
    )
    for patch, method in zip(box["boxes"], methods):
        patch.set_facecolor(COLORS.get(method, "0.7"))
        patch.set_alpha(0.55)
    for position, values_for_method in enumerate(values, start=1):
        if len(values_for_method):
            jitter = np.linspace(-0.08, 0.08, len(values_for_method))
            axis.scatter(
                position + jitter,
                values_for_method,
                s=13,
                color="black",
                alpha=0.55,
                zorder=3,
            )
            axis.text(
                position,
                np.nanmean(values_for_method) + 0.08,
                f"{np.nanmean(values_for_method):.2f}",
                ha="center",
                fontsize=8,
            )
    axis.axhline(0, color="0.6", linewidth=0.8, linestyle="--")
    axis.set_ylim(-1.1, 1.25)
    axis.set_title(f"{parameter} recovery", fontweight="bold")
    axis.set_ylabel("Gene-wise Spearman")
    axis.tick_params(axis="x", rotation=25)


def plot_organoid_branches(
    gene_metrics: pd.DataFrame, branch_genes: dict[str, list[str]], output: Path
) -> dict[str, pd.DataFrame]:
    """Render the curated organoid branch comparisons."""
    selected_tables = {}
    for branch, genes in branch_genes.items():
        subset = gene_metrics[
            gene_metrics["Dataset"].eq("organoid")
            & pd.to_numeric(gene_metrics["Branch"], errors="coerce").eq(int(branch))
            & gene_metrics["Gene"].isin(genes)
        ].copy()
        selected_tables[branch] = subset
        figure, axes = plt.subplots(1, 2, figsize=(10, 4.2))
        draw_boxplot(axes[0], subset, "alpha")
        draw_boxplot(axes[1], subset, "gamma")
        figure.suptitle(f"Organoid branch {branch}", fontweight="bold")
        figure.tight_layout()
        figure.savefig(
            output / f"organoid_branch{branch}_boxplot.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(figure)
    return selected_tables


def plot_kinetic_summary(summary: pd.DataFrame, output: Path) -> None:
    methods = ordered_methods(summary["Method"])
    parameters = [
        parameter
        for parameter in ("alpha", "gamma")
        if parameter in set(summary["Parameter"])
    ]
    positions = np.arange(len(methods))
    width = 0.7 / max(1, len(parameters))
    figure, axis = plt.subplots(figsize=(max(6, len(methods) * 1.2), 4.5))
    for index, parameter in enumerate(parameters):
        values = (
            summary[summary["Parameter"].eq(parameter)]
            .set_index("Method")["Mean"]
            .reindex(methods)
        )
        offset = (index - (len(parameters) - 1) / 2) * width
        axis.bar(positions + offset, values, width=width, label=parameter)
    axis.axhline(0, color="0.6", linewidth=0.8)
    axis.set_xticks(positions, methods, rotation=30)
    axis.set_ylabel("Equal-system mean Spearman")
    axis.set_title("Kinetics-recovery summary")
    axis.legend(title="Parameter")
    figure.tight_layout()
    figure.savefig(output / "kinetic_summary.png", dpi=300)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()

    metrics_dir = args.metrics.resolve()
    output = args.output.resolve()
    ensure_output(output)
    with args.config.open(encoding="utf-8") as handle:
        config = json.load(handle)
    gene_metrics, summary = read_inputs(metrics_dir)
    plot_data = output / "plot_data"
    plot_data.mkdir()

    for parameter in ("alpha", "gamma"):
        trajectories = plot_rpe1_trajectories(
            metrics_dir,
            gene_metrics,
            config["figures"]["trajectory_genes"],
            parameter,
            output,
        )
        trajectories.to_csv(
            plot_data / f"rpe1_{parameter}_trajectories.csv", index=False
        )

    heatmaps = plot_rpe1_heatmap(
        gene_metrics, config["figures"]["heatmap_genes"], output
    )
    for parameter, values in heatmaps.items():
        values.to_csv(plot_data / f"rpe1_{parameter}_heatmap.csv")

    for method in ("DeepVelo", "cellDancer"):
        classes = plot_kinetic_classes(method, metrics_dir, output)
        if not classes.empty:
            classes.to_csv(
                plot_data / f"rpe1_{method.lower()}_kinetic_classes.csv", index=False
            )

    branches = plot_organoid_branches(
        gene_metrics, config["figures"]["branch_genes"], output
    )
    for branch, values in branches.items():
        values.to_csv(
            plot_data / f"organoid_branch{branch}_gene_metrics.csv", index=False
        )

    plot_kinetic_summary(summary, output)
    print(f"Wrote kinetic-parameter figures to {output}")


if __name__ == "__main__":
    main()
