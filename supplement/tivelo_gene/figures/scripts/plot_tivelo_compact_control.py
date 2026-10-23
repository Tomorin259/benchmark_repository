"""Plot paired TIvelo simulation root-control scores on a fixed canvas."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from plot_style import BLUE, GREY, VARIANTS, style_axis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    score_path = args.data_root / "root_stage_intervention/root_stage_gene_scores.csv"
    summary_path = args.data_root / "root_stage_intervention/root_stage_summary.csv"
    scores = pd.read_csv(score_path, float_precision="round_trip")
    summaries = pd.read_csv(summary_path, float_precision="round_trip")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.7,
                         "axes.labelsize": 7.8, "xtick.labelsize": 6.8,
                         "ytick.labelsize": 6.8, "legend.fontsize": 6.7,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(6.82, 1.82), facecolor="white")
    ax = fig.add_axes([0.49 / 6.82, 0.39 / 1.82, 6.20 / 6.82, 1.10 / 1.82])
    groups = [("branching", "cbdir_reference30", "Branching\nCBDir-gene"),
              ("branching", "cossim", "Branching\nCosSim"),
              ("all", "cbdir_reference30", "Hybrid\nCBDir-gene"),
              ("all", "cossim", "Hybrid\nCosSim")]
    centers = [0, 1, 2.35, 3.35]
    rng = np.random.default_rng(20260926)
    plotted = []
    for center, (scenario, metric, label) in zip(centers, groups):
        frame = scores[(scores.scenario == scenario) & scores.variant.isin(VARIANTS)]
        values = frame.pivot(index="gene", columns="variant", values=metric).sort_index()
        expected_n = 20 if scenario == "branching" else 40
        if values.shape != (expected_n, 2) or not np.isfinite(values.to_numpy()).all():
            raise ValueError(f"Incomplete pairs: {scenario}/{metric}")
        jitter = rng.uniform(-0.046, 0.046, expected_n)
        for j, gene in enumerate(values.index):
            ax.plot([center - 0.14 + jitter[j], center + 0.14 + jitter[j]],
                    values.loc[gene, list(VARIANTS)], color="#A5ACB1",
                    alpha=0.20, lw=0.5, zorder=1)
            for variant in VARIANTS:
                plotted.append({"scenario": scenario, "metric": metric, "gene": gene,
                                "variant": variant, "value": float(values.loc[gene, variant])})
        for index, (variant, color, marker) in enumerate(zip(VARIANTS, (GREY, BLUE), ("o", "D"))):
            x = center + (-0.14 if index == 0 else 0.14)
            data = values[variant].to_numpy()
            ax.scatter(x + jitter, data, s=7, color=color, alpha=0.42,
                       edgecolor="none", marker=marker, zorder=2)
            mean = float(data.mean())
            expected = summaries[(summaries.scenario == scenario) & (summaries.variant == variant)][
                "mean_cbdir_reference30" if metric == "cbdir_reference30" else "mean_cossim"].item()
            np.testing.assert_allclose(mean, expected, rtol=0, atol=1e-10)
            ax.scatter(x, mean, s=40, color=color, marker=marker,
                       edgecolor="white", linewidth=0.7, zorder=4)
            ax.annotate(f"{mean:.2f}", (x, mean),
                        xytext=(-3 if index == 0 else 3, -13 if index == 0 else 9),
                        textcoords="offset points", ha="right" if index == 0 else "left",
                        va="center", color=color, fontsize=7)
    ax.axhline(0, color="#737B81", lw=0.7, zorder=0)
    ax.set(xlim=(-0.48, 3.82), ylim=(-0.8, 0.98), ylabel="Score")
    ax.set_xticks(centers, [group[2] for group in groups])
    ax.set_yticks([-0.5, 0, 0.5])
    ax.grid(axis="y", color="#EEF0F2", lw=0.5)
    handles = [Line2D([], [], ls="none", marker=marker, color=color, markersize=4.5, label=label)
               for color, marker, label in ((GREY, "o", "Original root"), (BLUE, "D", "Known early root"))]
    ax.legend(handles=handles, frameon=False, loc="lower left", bbox_to_anchor=(-0.02, 1.04),
              ncol=2, handletextpad=0.45, columnspacing=1.6)
    style_axis(ax)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(args.output_dir / f"tivelo_root_control_gene_metrics.{suffix}", dpi=300)
    plt.close(fig)
    pd.DataFrame(plotted).to_csv(args.output_dir / "simulation_plotted_pairs.csv", index=False)
    manifest = {"generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "sources": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (score_path, summary_path)},
                "canvas_inches": [6.82, 1.82], "numerical_changes": False,
                "new_titles_captions_panel_letters": False}
    (args.output_dir / "simulation_source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
