"""Axis styling shared by the retained simulation and bone marrow panels."""
import matplotlib.pyplot as plt

GREY = "#777D82"
BLUE = "#0072B2"
VARIANTS = ("correct_path_root", "oracle_early_root")


def export(fig, stem):
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(stem.with_suffix(f".{suffix}"), dpi=300,
                    bbox_inches="tight", pad_inches=0.035, facecolor="white")
    plt.close(fig)


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#9EA3A7")
        ax.spines[spine].set_linewidth(0.65)
    ax.tick_params(length=2.5, width=0.65, pad=3)
    ax.set_axisbelow(True)
