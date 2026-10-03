#!/usr/bin/env python3
"""Render Figure 6 quantitative panels from one kinetics evaluator output."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
try:
    from matplotlib_venn import venn2, venn2_circles
except ModuleNotFoundError:
    venn2 = venn2_circles = None
import numpy as np
import pandas as pd
import seaborn as sns


DISPLAY_METHOD_NAMES = {
    "DeepVelo": "DeepVelo",
    "Dynamo": "Dynamo(label data)",
    "cellDancer": "cellDancer",
    "DeepKINET": "DeepKINET",
    "TSvelo": "TSvelo",
}

TRAJECTORY_COLORS = {
    # Figure 6B follows the benchmark hues. Separate the two green methods
    # by lightness so their overlapping trajectories remain distinguishable.
    "DeepVelo": "#C9383B",
    "cellDancer": "#80C46F",
    "Dynamo(label data)": "#126B36",
    "TSvelo": "#A05D56",
    "DeepKINET": "#FA8F12",
}

BOXPLOT_METHOD_COLORS = {
    # Figure 6F follows evaluation.cell.plots.BENCHMARK_COLORS. The labelled
    # Dynamo control reuses the Dynamo(RNA_seq) method color.
    "DeepVelo": "#C9383B",
    "cellDancer": "#5BA658",
    "Dynamo(label data)": "#2CA02C",
    "DeepKINET": "#FA8F12",
}


def canonical_method_color(method: str) -> str:
    """Return the benchmark color for a Figure 6F method key."""
    display_name = DISPLAY_METHOD_NAMES.get(method, method)
    if display_name not in BOXPLOT_METHOD_COLORS:
        raise KeyError(f"No benchmark boxplot color registered for {method}")
    return BOXPLOT_METHOD_COLORS[display_name]

# Evaluator IDs and exported plot-data stay unchanged; this is a figure label.
ORGANOID_PLOT_METHOD_NAMES = {"Dynamo": "Dynamo"}
FIGURE_METHOD_NAMES = {"Dynamo": "Dynamo (scEU)", "Dynamo(label data)": "Dynamo (scEU)"}

# Display-only panels inherit the relative method order of Figure 7.  The
# metric-ranked heatmap remains sorted by its own correlation values.
PAPER_METHOD_ORDER = ("Dynamo", "DeepKINET", "TSvelo", "cellDancer", "DeepVelo")


def paper_method_order(paths: dict) -> dict:
    return {method: paths[method] for method in PAPER_METHOD_ORDER if method in paths}


def resolve(root: Path, relative_path: str) -> Path:
    path = Path(relative_path)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Figure 6 paths must be relative: {relative_path}")
    resolved = root / path
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    return resolved


def ensure_output(path: Path) -> None:
    if path.exists() and not path.is_dir():
        raise NotADirectoryError(path)
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Output directory must be empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def read_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.columns = pd.Index(range(frame.shape[1]))
    return frame.apply(pd.to_numeric, errors="coerce")


def read_metric(
    path: Path,
    method: str,
    parameter: str,
    genes: list[str] | None = None,
    dataset: str = "rpe1",
    branch: int | None = None,
) -> pd.DataFrame:
    frame = pd.read_csv(path)
    if {"Method", "Parameter", "Dataset", "Branch"}.issubset(frame.columns):
        selected = (
            frame["Method"].eq(DISPLAY_METHOD_NAMES.get(method, method))
            & frame["Parameter"].str.lower().eq(parameter.lower())
            & frame["Dataset"].eq(dataset)
        )
        selected &= (
            frame["Branch"].isna() if branch is None
            else pd.to_numeric(frame["Branch"], errors="raise").eq(branch)
        )
        frame = frame.loc[selected]
    if frame.empty or frame["Gene"].duplicated().any():
        raise ValueError(f"{path}: missing or duplicated rows for {method}/{dataset}/{branch}/{parameter}")
    correlation_column = next(
        column for column in frame.columns if "corr" in column.lower()
    )
    if genes is not None:
        frame = frame[frame["Gene"].isin(genes)]
    return pd.DataFrame(
        {
            "Gene": frame["Gene"],
            "Correlation": frame[correlation_column],
            "Method": DISPLAY_METHOD_NAMES.get(method, method),
            "Parameter": parameter,
        }
    )


def zscore(values):
    series = pd.Series(values, dtype=float)
    standard_deviation = series.std()
    if pd.isna(standard_deviation) or standard_deviation == 0:
        return series - series.mean()
    return (series - series.mean()) / standard_deviation


def plot_trajectory(
    parameter: str,
    reference_path: Path,
    method_paths: dict[str, Path],
    genes: list[str],
    output_path: Path,
) -> pd.DataFrame:
    """Preserve the 10-column trajectory layout and export plotted values."""
    truth = read_matrix(reference_path)
    methods = {name: read_matrix(path) for name, path in paper_method_order(method_paths).items()}
    valid_genes = [gene for gene in genes if gene in truth.index]
    plotted_rows = []
    rows = math.ceil(len(valid_genes) / 10)
    figure, axes = plt.subplots(rows, 10, figsize=(30, 3 * rows))
    axes = np.asarray(axes).reshape(-1)

    for index, axis in enumerate(axes):
        if index >= len(valid_genes):
            axis.axis("off")
            continue
        gene = valid_genes[index]
        truth_values = zscore(truth.loc[gene])
        plotted_rows.extend(
            {"Gene": gene, "Source": "scEU-seq", "Relative_Position": t, "Display_Z": v}
            for t, v in zip(np.linspace(0, 1, len(truth_values)), truth_values)
        )
        axis.plot(
            np.linspace(0, 1, len(truth_values)),
            truth_values,
            color="black",
            linewidth=2.5,
            zorder=10,
            label="Ground Truth",
        )
        for name, frame in methods.items():
            if gene not in frame.index:
                continue
            values = zscore(frame.loc[gene])
            plotted_rows.extend(
                {"Gene": gene, "Source": DISPLAY_METHOD_NAMES.get(name, name),
                 "Relative_Position": t, "Display_Z": v}
                for t, v in zip(np.linspace(0, 1, len(values)), values)
            )
            axis.plot(
                np.linspace(0, 1, len(values)),
                values,
                color=TRAJECTORY_COLORS[DISPLAY_METHOD_NAMES.get(name, name)],
                linewidth=1.5,
                alpha=0.8,
                label=FIGURE_METHOD_NAMES.get(name, name),
            )
        axis.grid(False)
        axis.set_xlabel("")
        axis.set_ylabel("")
        axis.set_title(gene, fontsize=12, fontweight="bold")

    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].axis("off")
    axes[-1].legend(handles, labels, loc="center", frameon=False, fontsize=12)
    plt.tight_layout()
    figure.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close(figure)
    return pd.DataFrame(plotted_rows)


def plot_heatmap(
    alpha_paths: dict[str, Path],
    gamma_paths: dict[str, Path],
    genes: list[str],
    output_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = {}
    for parameter, paths in (("Alpha", alpha_paths), ("Gamma", gamma_paths)):
        combined = pd.concat(
            [
                read_metric(path, method, parameter, genes)
                for method, path in paths.items()
            ],
            ignore_index=True,
        )
        matrix = combined.pivot(index="Method", columns="Gene", values="Correlation")
        ordered_genes = [gene for gene in genes if gene in matrix.columns]
        matrix = matrix[ordered_genes]
        data[parameter] = matrix.loc[
            matrix.mean(axis=1).sort_values(ascending=False).index
        ]

    width = max(24, max(len(frame.columns) for frame in data.values()) * 1.35)
    height = max(len(frame.index) for frame in data.values()) * 0.8 + 3
    figure, axes = plt.subplots(1, 2, figsize=(width, height))
    plt.subplots_adjust(left=0.12, wspace=0.24, right=0.895)
    colorbar_axis = figure.add_axes([0.955, 0.3, 0.007, 0.4])
    summary_columns = []

    for index, parameter in enumerate(("Alpha", "Gamma")):
        sns.heatmap(
            data[parameter],
            ax=axes[index],
            cmap="RdBu_r",
            center=0,
            vmin=-0.5,
            vmax=1,
            annot=True,
            fmt=".2f",
            linewidths=1,
            linecolor="white",
            square=True,
            cbar=index == 1,
            cbar_ax=colorbar_axis if index == 1 else None,
            cbar_kws={"label": "Spearman Correlation"} if index == 1 else None,
        )
        axes[index].set_anchor("N")
        axes[index].set_title(
            f"{parameter} Performance", fontsize=16, fontweight="bold", pad=15
        )
        axes[index].set_xlabel("Gene", fontsize=14, fontweight="bold")
        axes[index].set_ylabel("Method" if index == 0 else "", fontsize=14, fontweight="bold")
        axes[index].set_xticklabels(
            axes[index].get_xticklabels(), rotation=45, ha="right", fontsize=12
        )
        frame = data[parameter]
        labels = [
            FIGURE_METHOD_NAMES.get(method, method)
            for method in frame.index
        ]
        axes[index].set_yticklabels(labels, rotation=0, fontsize=12)
        # The ordering is descriptive: methods cover different display genes.
        header = axes[index].annotate(
            "Mean (n)", xy=(1.07, 1), xycoords="axes fraction",
            xytext=(0, 15), textcoords="offset points",
            fontsize=12, fontweight="bold", ha="center", va="baseline",
            annotation_clip=False,
        )
        mean_labels, count_labels = [], []
        for row, (_, scores) in enumerate(frame.iterrows()):
            mean_labels.append(axes[index].text(
                1.082, row + 0.5, f"{scores.mean():.2f}",
                transform=axes[index].get_yaxis_transform(),
                fontsize=12, fontfamily="DejaVu Sans", fontweight="bold",
                ha="right", va="center",
            ))
            count_labels.append(axes[index].text(
                1.096, row + 0.5, f"({scores.count()})",
                transform=axes[index].get_yaxis_transform(),
                fontsize=12, fontfamily="DejaVu Sans", fontweight="bold",
                ha="left", va="center",
            ))
        summary_columns.append((header, mean_labels, count_labels))

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    # Use the same physical column widths and padding in both panels.
    mean_width = max(
        label.get_window_extent(renderer).width
        for _, means, _ in summary_columns for label in means
    )
    count_width = max(
        label.get_window_extent(renderer).width
        for _, _, counts in summary_columns for label in counts
    )
    point = figure.dpi / 72
    for axis, (header, means, counts) in zip(axes, summary_columns):
        column_left = axis.bbox.x1 + 8 * point
        mean_right = column_left + mean_width
        count_left = mean_right + 4 * point
        inverse = axis.transAxes.inverted()
        for label in means:
            label.set_x(inverse.transform((mean_right, 0))[0])
        for label in counts:
            label.set_x(inverse.transform((count_left, 0))[0])
        header.xy = (inverse.transform(
            ((column_left + count_left + count_width) / 2, 0)
        )[0], 1)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    # Reserve only the measured width of the mean column and next row labels.
    # This keeps the gap compact without overlapping text at export size.
    mean_right = max(text.get_window_extent(renderer).x1 for text in axes[0].texts)
    gamma_label_width = max(
        label.get_window_extent(renderer).width for label in axes[1].get_yticklabels()
    )
    gap_pixels = 10 * figure.dpi / 72
    right_bounds = axes[1].get_position()
    right_start = (mean_right + gamma_label_width + gap_pixels) / figure.bbox.width
    axes[1].set_position([right_start, right_bounds.y0, right_bounds.width, right_bounds.height])
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    gamma_mean_right = max(text.get_window_extent(renderer).x1 for text in axes[1].texts)
    colorbar_left = (gamma_mean_right + gap_pixels) / figure.bbox.width
    heatmap_bounds = axes[0].get_position()
    colorbar_axis.set_position([colorbar_left, heatmap_bounds.y0, 0.007, heatmap_bounds.height])
    figure.legend(
        handles=[
            mpatches.Patch(
                facecolor="white",
                label="Not Calculated",
                edgecolor="gray",
                linewidth=1,
            )
        ],
        loc="upper left",
        bbox_to_anchor=(colorbar_left - 0.01, heatmap_bounds.y0 - 0.02),
        frameon=False,
        fontsize=10,
    )
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
    return data["Alpha"], data["Gamma"]


def cosine_classes(alpha_path: Path, gamma_path: Path) -> dict[str, set[str]]:
    alpha = read_matrix(alpha_path)
    gamma = read_matrix(gamma_path)
    genes = alpha.index.intersection(gamma.index)
    columns = alpha.columns.intersection(gamma.columns)
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


def plot_venn(
    reference_alpha: Path,
    reference_gamma: Path,
    prediction_alpha: Path,
    prediction_gamma: Path,
    output_path: Path,
    prediction_label: str,
    membership_output_path: Path,
) -> pd.DataFrame:
    if venn2 is None or venn2_circles is None:
        raise RuntimeError("matplotlib-venn is required when rendering Venn panels")
    truth = cosine_classes(reference_alpha, reference_gamma)
    prediction = cosine_classes(prediction_alpha, prediction_gamma)
    intervals = {
        "Cooperative": r"(cosine $\leq -0.5$)",
        "Neutral": r"($-0.5 <$ cosine $< 0.5$)",
        "Destabilizing": r"(cosine $\geq 0.5$)",
    }
    figure = plt.figure(figsize=(15, 5))
    rows = []
    for index, category in enumerate(
        ("Cooperative", "Neutral", "Destabilizing"), start=1
    ):
        left = truth[category]
        right = prediction[category]
        subsets = (len(left - right), len(right - left), len(left & right))
        axis = figure.add_subplot(1, 3, index)
        venn2(
            subsets=subsets,
            set_labels=("scEU", prediction_label),
            set_colors=("skyblue", "salmon"),
            alpha=0.6,
            ax=axis,
        )
        venn2_circles(subsets=subsets, linewidth=1.5, color="gray", ax=axis)
        # Figure coordinates keep both title lines aligned across all circles.
        figure.text((index - 0.5) / 3, 0.95, f"{category} Genes",
                    fontsize=12, fontweight="bold", ha="center", va="top")
        figure.text((index - 0.5) / 3, 0.89, intervals[category],
                    fontsize=11, ha="center", va="top")
        rows.append(
            {
                "Class": category,
                "scEU_only": subsets[0],
                f"{prediction_label}_only": subsets[1],
                "intersection": subsets[2],
            }
        )
    membership_rows = []
    for category in ("Cooperative", "Neutral", "Destabilizing"):
        for gene in sorted(truth[category] | prediction[category]):
            membership_rows.append(
                {
                    "Class": category,
                    "Gene": gene,
                    "scEU": gene in truth[category],
                    prediction_label: gene in prediction[category],
                }
            )
    pd.DataFrame(membership_rows).to_csv(membership_output_path, index=False)
    plt.tight_layout(rect=[0, 0.03, 1, 0.84])
    figure.savefig(output_path, dpi=300)
    plt.close(figure)
    return pd.DataFrame(rows)


def plot_branch(
    branch_id: int,
    method_paths: dict[str, dict[str, Path]],
    genes: list[str],
    output_path: Path,
) -> pd.DataFrame:
    method_paths = paper_method_order(method_paths)
    frames = []
    for method, paths in method_paths.items():
        display_name = ORGANOID_PLOT_METHOD_NAMES.get(
            method, DISPLAY_METHOD_NAMES.get(method, method)
        )
        display_method = f"{display_name}_branch{branch_id}"
        for parameter in ("Alpha", "Gamma"):
            path = paths.get(parameter)
            if path is not None:
                frame = read_metric(path, method, parameter, genes, dataset="organoid", branch=branch_id)
                frame["Method"] = display_method
                frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)

    figure, axes = plt.subplots(1, 2, figsize=(18, 8))
    palette = {
        f"{ORGANOID_PLOT_METHOD_NAMES.get(method, DISPLAY_METHOD_NAMES.get(method, method))}_branch{branch_id}": paths["color"]
        for method, paths in method_paths.items()
    }
    for axis, parameter in zip(axes, ("Alpha", "Gamma")):
        subset = combined[combined["Parameter"] == parameter]
        sns.boxplot(
            x="Method",
            y="Correlation",
            data=subset,
            ax=axis,
            palette=palette,
            width=0.6,
            boxprops=dict(alpha=0.5, linewidth=1.5),
            whiskerprops=dict(linewidth=1.5),
            capprops=dict(linewidth=1.5),
            medianprops=dict(color="black", linewidth=2),
            fliersize=0,
        )
        np.random.seed(0)  # Display jitter only; scores are unchanged.
        sns.stripplot(
            x="Method",
            y="Correlation",
            data=subset,
            ax=axis,
            color="black",
            size=7,
            alpha=0.7,
            jitter=True,
        )
        axis.set_title(
            f"Performance on {parameter}", fontsize=20, fontweight="bold", pad=15
        )
        axis.set_ylabel("Spearman Correlation (R)", fontsize=16, labelpad=10)
        axis.set_xlabel("")
        axis.tick_params(axis="x", labelsize=14, rotation=0)
        axis.set_xticklabels([
            tick.get_text().replace("Dynamo_", "Dynamo (scEU)_")
            for tick in axis.get_xticklabels()
        ])
        axis.tick_params(axis="y", labelsize=14)
        axis.set_ylim(-1.1, 1.25)
        axis.axhline(0, color="gray", linestyle="--", linewidth=1.5, alpha=0.6)
        means = subset.groupby("Method")["Correlation"].mean()
        # Boxplot positions are compacted separately for Alpha and Gamma when
        # a method has no value for one parameter (Dynamo has no Gamma here).
        # Derive annotation positions from the categories actually drawn so a
        # missing parameter cannot push another method's label outside the axes.
        drawn_methods = [
            f"{ORGANOID_PLOT_METHOD_NAMES.get(method, DISPLAY_METHOD_NAMES.get(method, method))}_branch{branch_id}"
            for method in method_paths
            if f"{ORGANOID_PLOT_METHOD_NAMES.get(method, DISPLAY_METHOD_NAMES.get(method, method))}_branch{branch_id}"
            in set(subset["Method"])
        ]
        method_positions = {method: position for position, method in enumerate(drawn_methods)}
        for method in method_paths:
            display_name = ORGANOID_PLOT_METHOD_NAMES.get(
                method, DISPLAY_METHOD_NAMES.get(method, method)
            )
            display_method = f"{display_name}_branch{branch_id}"
            if display_method in means:
                axis.text(
                    method_positions[display_method],
                    means[display_method] + 0.08,
                    f"{means[display_method]:.4f}",
                    ha="center",
                    color="black",
                    fontweight="bold",
                    fontsize=13,
                    bbox=dict(
                        facecolor="white", alpha=0.7, edgecolor="none", pad=1
                    ),
                )

    plt.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
    return combined


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-venn", action="store_true",
                        help="Redraw display panels without the unchanged Venn diagrams")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).with_name("figure6_config.json"),
    )
    args = parser.parse_args()

    artifact_root = args.artifact_root.expanduser().resolve()
    output = args.output.expanduser().resolve()
    ensure_output(output)
    config = json.loads(args.config.resolve().read_text(encoding="utf-8"))
    plot_data = output / "plot_data"
    plot_data.mkdir(parents=True, exist_ok=True)

    rpe1 = config["inputs"]["rpe1"]
    for parameter in ("alpha", "gamma"):
        trajectory_data = plot_trajectory(
            parameter,
            resolve(artifact_root, rpe1["reference"][parameter]),
            {
                method: resolve(artifact_root, path)
                for method, path in rpe1["trajectories"][parameter].items()
            },
            config["trajectory_genes"],
            output / f"rpe1_{parameter}_trajectories.png",
        )
        trajectory_data.to_csv(plot_data / f"rpe1_{parameter}_trajectories.csv", index=False)

    alpha_matrix, gamma_matrix = plot_heatmap(
        {
            method: resolve(artifact_root, path)
            for method, path in rpe1["heatmap_metrics"]["Alpha"].items()
        },
        {
            method: resolve(artifact_root, path)
            for method, path in rpe1["heatmap_metrics"]["Gamma"].items()
        },
        config["heatmap_genes"],
        output / "rpe1_alpha_gamma_heatmap.png",
    )
    alpha_matrix.to_csv(plot_data / "rpe1_alpha_heatmap.csv")
    gamma_matrix.to_csv(plot_data / "rpe1_gamma_heatmap.csv")

    if not args.skip_venn:
        venn = rpe1["deepvelo_venn"]
        venn_rows = plot_venn(
            resolve(artifact_root, venn["reference_alpha"]),
            resolve(artifact_root, venn["reference_gamma"]),
            resolve(artifact_root, venn["prediction_alpha"]),
            resolve(artifact_root, venn["prediction_gamma"]),
            output / "rpe1_deepvelo_kinetic_classes.png",
            "DeepVelo",
            plot_data / "rpe1_deepvelo_kinetic_class_membership.csv",
        )
        venn_rows.to_csv(plot_data / "rpe1_deepvelo_kinetic_classes.csv", index=False)
        cell_dancer_venn = rpe1["cellDancer_venn"]
        cell_dancer_rows = plot_venn(
            resolve(artifact_root, cell_dancer_venn["reference_alpha"]),
            resolve(artifact_root, cell_dancer_venn["reference_gamma"]),
            resolve(artifact_root, cell_dancer_venn["prediction_alpha"]),
            resolve(artifact_root, cell_dancer_venn["prediction_gamma"]),
            output / "rpe1_celldancer_kinetic_classes.png",
            "cellDancer",
            plot_data / "rpe1_celldancer_kinetic_class_membership.csv",
        )
        cell_dancer_rows.to_csv(
            plot_data / "rpe1_celldancer_kinetic_classes.csv", index=False
        )

    for branch_id, branch_config in config["inputs"]["organoid"]["branches"].items():
        method_paths = {
            method: {
                "color": canonical_method_color(method),
                **{
                    parameter: resolve(artifact_root, path)
                    for parameter, path in values.items()
                    if parameter != "color"
                },
            }
            for method, values in branch_config["methods"].items()
        }
        branch_data = plot_branch(
            int(branch_id),
            method_paths,
            branch_config["genes"],
            output / f"organoid_branch{branch_id}_boxplot.png",
        )
        branch_data.to_csv(
            plot_data / f"organoid_branch{branch_id}_gene_metrics.csv", index=False
        )

    print(f"Wrote Figure 6 quantitative panels from {artifact_root} to {output}")


if __name__ == "__main__":
    main()
