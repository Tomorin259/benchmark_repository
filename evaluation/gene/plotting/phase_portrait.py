"""Reproducible gene-level spliced--unspliced phase portraits.

The sampling idea is adapted from the cellDancer cdplt.scatter_gene API and
sampling_neighbors (BSD 3-Clause; https://github.com/GuangyuWangLab2021/cellDancer).
This implementation is maintained in the benchmark repository, does not
import cellDancer, and changes only how arrows are displayed. The source
AnnData object and all metric inputs remain untouched.

The main figure contract is intentionally conservative:

* phase coordinates stay linear and use the method-native layers;
* real cells are sampled deterministically from a normalized phase grid;
* velocity direction is preserved exactly;
* arrow length uses panel-wide p90 normalization and log1p compression;
* no neighborhood velocity smoothing is applied in the main figure.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors


PHASE_PLOT_VERSION = "gene-phase-log1p-005-010-v1"


def _dense_gene_vector(adata: Any, gene: str, layer: str) -> np.ndarray:
    """Read one gene without changing the source AnnData object."""
    if layer == "X":
        values = adata[:, [gene]].X
    elif layer in adata.layers:
        values = adata[:, [gene]].layers[layer]
    else:
        raise KeyError(f"Layer {layer!r} is not present in the AnnData object.")

    if sparse.issparse(values):
        values = values.toarray()
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    if values.size != adata.n_obs:
        raise ValueError(
            f"Layer {layer!r} for {gene!r} has {values.size} values; "
            f"expected {adata.n_obs}."
        )
    return values


def _axis_limits(
    values: np.ndarray, pad_fraction: float = 0.025
) -> Tuple[Tuple[float, float], float]:
    finite = np.isfinite(values)
    if not finite.any():
        raise ValueError("A phase coordinate contains no finite values.")

    lower = float(np.min(values[finite]))
    upper = float(np.max(values[finite]))
    span = upper - lower
    if not np.isfinite(span):
        raise ValueError("A phase coordinate has a non-finite range.")
    if span == 0:
        span = max(abs(lower) * 0.05, 1.0)
        lower -= 0.5 * span
        upper += 0.5 * span
    else:
        pad = pad_fraction * span
        lower -= pad
        upper += pad
    return (lower, upper), upper - lower


def _category_colors(
    adata: Any,
    color_key: Optional[str],
    valid: np.ndarray,
) -> Tuple[Any, Optional[str], Sequence[Line2D]]:
    """Return scatter colors and optional legend handles."""
    if color_key is None or color_key not in adata.obs:
        return "#95D9EF", None, []

    groups = adata.obs[color_key]
    if pd.api.types.is_numeric_dtype(groups):
        return np.asarray(groups)[valid], "viridis", []

    if hasattr(groups, "cat"):
        categories = list(groups.cat.categories)
    else:
        categories = list(pd.unique(np.asarray(groups)))

    palette = getattr(adata, "uns", {}).get(f"{color_key}_colors")
    if palette is None or len(palette) < len(categories):
        cmap = plt.get_cmap("tab20", max(1, len(categories)))
        palette = [cmap(i) for i in range(len(categories))]
    else:
        palette = list(palette)

    color_map = dict(zip(categories, palette))
    values = np.asarray(
        [color_map.get(value, "#808080") for value in groups], dtype=object
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=str(category),
            markerfacecolor=color_map[category],
            markersize=7,
        )
        for category in categories
    ]
    return values[valid], None, handles


def _sample_phase_cells(
    coordinates: np.ndarray,
    cell_ids: Sequence[str],
    *,
    grid: Tuple[int, int] = (25, 25),
    density_percentile: float = 15.0,
) -> np.ndarray:
    """Select representative real cells without random perturbation.

    The public CellDancer API uses a jittered grid and then selects nearest
    cells. We retain that scientific intent but make the selection stable:
    grid and density calculations happen in normalized phase coordinates,
    and cell IDs provide a stable tie-breaking order.
    """
    coordinates = np.asarray(coordinates, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("coordinates must have shape (n_cells, 2).")
    if len(cell_ids) != len(coordinates):
        raise ValueError("cell_ids and coordinates must have the same length.")
    if len(coordinates) == 0:
        return np.empty(0, dtype=np.int64)
    if not 0 <= density_percentile <= 100:
        raise ValueError("density_percentile must be between 0 and 100.")

    nx, ny = (int(grid[0]), int(grid[1]))
    if nx < 1 or ny < 1:
        raise ValueError("grid dimensions must be positive integers.")

    ids = np.asarray([str(value) for value in cell_ids], dtype=str)
    order = np.argsort(ids, kind="mergesort")
    sorted_coordinates = coordinates[order]
    finite = np.isfinite(sorted_coordinates).all(axis=1)
    if not finite.any():
        return np.empty(0, dtype=np.int64)

    sorted_coordinates = sorted_coordinates[finite]
    sorted_order = order[finite]
    lower = sorted_coordinates.min(axis=0)
    span = sorted_coordinates.max(axis=0) - lower
    safe_span = np.where(span > 0, span, 1.0)
    normalized = (sorted_coordinates - lower) / safe_span

    grid_x = np.linspace(0.0, 1.0, nx)
    grid_y = np.linspace(0.0, 1.0, ny)
    mesh_x, mesh_y = np.meshgrid(grid_x, grid_y, indexing="xy")
    grid_points = np.column_stack((mesh_x.ravel(), mesh_y.ravel()))

    nearest = NearestNeighbors(n_neighbors=1).fit(normalized)
    _, nearest_indices = nearest.kneighbors(grid_points)
    candidates = sorted_order[nearest_indices[:, 0]]

    candidate_ids = ids[candidates]
    _, first_positions = np.unique(candidate_ids, return_index=True)
    candidates = candidates[np.sort(first_positions)]
    if len(candidates) == 0 or len(normalized) < 2:
        return candidates.astype(np.int64)

    n_neighbors = min(20, len(normalized) - 1)
    density_neighbors = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(normalized)
    distances, _ = density_neighbors.kneighbors(normalized)
    distances = distances[:, 1:]
    grid_spacing = np.mean(
        [
            1.0 / max(nx - 1, 1),
            1.0 / max(ny - 1, 1),
        ]
    )
    sigma = max(3.0 * grid_spacing, np.finfo(float).eps)
    density = np.exp(-0.5 * (distances / sigma) ** 2).sum(axis=1)

    sorted_position = {
        int(original): position
        for position, original in enumerate(sorted_order)
    }
    candidate_density = np.asarray(
        [density[sorted_position[int(candidate)]] for candidate in candidates],
        dtype=float,
    )
    threshold = np.percentile(candidate_density, density_percentile)
    keep = candidate_density >= threshold
    if not keep.any():
        keep[:] = True
    return candidates[keep].astype(np.int64)


def _scale_phase_vectors(
    velocities: np.ndarray,
    axis_ranges: Tuple[float, float],
    *,
    mode: str = "log1p",
    reference_quantile: float = 0.90,
    reference_length: float = 0.05,
    max_length: float = 0.10,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Compress vector lengths while preserving every vector direction."""
    vectors = np.asarray(velocities, dtype=np.float64)
    if vectors.ndim != 2 or vectors.shape[1] != 2:
        raise ValueError("velocities must have shape (n_cells, 2).")
    ranges = np.asarray(axis_ranges, dtype=np.float64)
    if ranges.shape != (2,) or not np.isfinite(ranges).all() or (ranges <= 0).any():
        raise ValueError("axis_ranges must contain two positive finite values.")
    if mode not in {"linear", "log1p"}:
        raise ValueError("mode must be 'linear' or 'log1p'.")
    if not 0 < reference_quantile <= 1:
        raise ValueError("reference_quantile must be in (0, 1].")
    if not 0 < reference_length <= max_length:
        raise ValueError("reference_length must be positive and <= max_length.")
    if max_length <= 0:
        raise ValueError("max_length must be positive.")

    display = np.zeros_like(vectors)
    finite = np.isfinite(vectors).all(axis=1)
    normalized = np.zeros_like(vectors)
    normalized[finite] = vectors[finite] / ranges
    magnitudes = np.linalg.norm(normalized, axis=1)
    positive = finite & (magnitudes > 0)
    if not positive.any():
        return display, {
            "reference_magnitude": None,
            "input_magnitude_max": 0.0,
            "display_magnitude_max": 0.0,
            "n_capped": 0,
            "capped_fraction": 0.0,
        }

    reference = float(np.quantile(magnitudes[positive], reference_quantile))
    reference = max(reference, np.finfo(float).eps)
    ratio = magnitudes[positive] / reference
    if mode == "log1p":
        lengths = reference_length * np.log1p(ratio) / np.log(2.0)
    else:
        lengths = reference_length * ratio
    capped = lengths > max_length
    lengths = np.minimum(lengths, max_length)
    direction = normalized[positive] / magnitudes[positive, None]
    display[positive] = direction * lengths[:, None] * ranges

    return display, {
        "reference_magnitude": reference,
        "input_magnitude_max": float(np.max(magnitudes[positive])),
        "display_magnitude_max": float(np.max(lengths)),
        "n_capped": int(np.count_nonzero(capped)),
        "capped_fraction": float(np.mean(capped)),
    }


def plot_phase_portrait(
    adata: Any,
    gene: str,
    ax: Optional[Any] = None,
    layer_s: str = "Ms",
    layer_u: str = "Mu",
    layer_v_s: str = "velocity",
    layer_v_u: str = "velocity_u",
    color_key: Optional[str] = "clusters",
    s: float = 20,
    alpha: float = 0.6,
    arrow_grid: Tuple[int, int] = (25, 25),
    legend_loc: Optional[str] = "upper right",
    *,
    arrow_length_mode: str = "log1p",
    arrow_reference_quantile: float = 0.90,
    arrow_reference_length: float = 0.05,
    arrow_max_length: float = 0.10,
    arrow_density_percentile: float = 15.0,
    return_diagnostics: bool = False,
) -> Any:
    """Plot one gene's phase portrait with a stable, bounded arrow field.

    The displayed arrows are a visualization transform only. The source
    velocity layers and all values used by gene-level metrics are unchanged.
    Both velocity components are required because a missing component would
    otherwise silently create a misleading horizontal or vertical field.
    """
    if gene not in adata.var_names:
        raise ValueError(f"Gene {gene!r} is not present in the AnnData object.")
    if layer_v_u is None:
        raise ValueError("layer_v_u must name an unspliced velocity layer.")

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))

    spliced = _dense_gene_vector(adata, gene, layer_s)
    unspliced = _dense_gene_vector(adata, gene, layer_u)
    velocity_s = _dense_gene_vector(adata, gene, layer_v_s)
    velocity_u = _dense_gene_vector(adata, gene, layer_v_u)
    coordinates = np.column_stack((spliced, unspliced))
    velocities = np.column_stack((velocity_s, velocity_u))
    coordinate_valid = np.isfinite(coordinates).all(axis=1)
    velocity_valid = np.isfinite(velocities).all(axis=1)
    if not coordinate_valid.any():
        raise ValueError(f"Gene {gene!r} has no finite phase coordinates.")

    xlim, x_range = _axis_limits(spliced[coordinate_valid])
    ylim, y_range = _axis_limits(unspliced[coordinate_valid])
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_box_aspect(1)

    colors, cmap, legend_handles = _category_colors(adata, color_key, coordinate_valid)
    ax.scatter(
        spliced[coordinate_valid],
        unspliced[coordinate_valid],
        c=colors,
        cmap=cmap,
        s=s,
        alpha=alpha,
        edgecolors="none",
    )
    if legend_loc and legend_handles:
        ax.legend(
            handles=legend_handles,
            loc=legend_loc,
            bbox_to_anchor=(1.05, 1),
            frameon=False,
        )

    obs_names = np.asarray([str(value) for value in adata.obs_names], dtype=str)
    sample_indices = _sample_phase_cells(
        coordinates,
        obs_names,
        grid=arrow_grid,
        density_percentile=arrow_density_percentile,
    )
    display_velocities, scale_diagnostics = _scale_phase_vectors(
        velocities,
        (x_range, y_range),
        mode=arrow_length_mode,
        reference_quantile=arrow_reference_quantile,
        reference_length=arrow_reference_length,
        max_length=arrow_max_length,
    )
    draw_indices = sample_indices[velocity_valid[sample_indices]]
    if len(draw_indices):
        ax.scatter(
            spliced[draw_indices],
            unspliced[draw_indices],
            color="none",
            edgecolor="k",
            s=s,
            linewidth=0.8,
            zorder=19,
        )
        ax.quiver(
            spliced[draw_indices],
            unspliced[draw_indices],
            display_velocities[draw_indices, 0],
            display_velocities[draw_indices, 1],
            angles="xy",
            scale_units="xy",
            scale=1,
            units="width",
            width=0.002,
            headwidth=4.5,
            headlength=5,
            headaxislength=4.5,
            color="k",
            alpha=1.0,
            zorder=20,
        )

    ax.set_xlabel(f"Spliced ({layer_s})")
    ax.set_ylabel(f"Unspliced ({layer_u})")
    ax.set_title(f"Gene: {gene}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    diagnostics = {
        "plot_version": PHASE_PLOT_VERSION,
        "gene": gene,
        "layer_s": layer_s,
        "layer_u": layer_u,
        "layer_v_s": layer_v_s,
        "layer_v_u": layer_v_u,
        "xlim": [float(xlim[0]), float(xlim[1])],
        "ylim": [float(ylim[0]), float(ylim[1])],
        "n_obs": int(adata.n_obs),
        "n_coordinate_valid": int(np.count_nonzero(coordinate_valid)),
        "n_velocity_valid": int(np.count_nonzero(coordinate_valid & velocity_valid)),
        "n_sampled": int(len(sample_indices)),
        "n_drawn": int(len(draw_indices)),
        "sampled_cell_ids": obs_names[sample_indices].tolist(),
        "drawn_cell_ids": obs_names[draw_indices].tolist(),
        "arrow_grid": [int(arrow_grid[0]), int(arrow_grid[1])],
        "arrow_density_percentile": float(arrow_density_percentile),
        "arrow_length_mode": arrow_length_mode,
        "arrow_reference_quantile": float(arrow_reference_quantile),
        "arrow_reference_length": float(arrow_reference_length),
        "arrow_max_length": float(arrow_max_length),
        **scale_diagnostics,
    }
    return (ax, diagnostics) if return_diagnostics else ax
