"""BATC metric computation utilities.

This module packages the ad-hoc BATC (Branch-Aware Tangent Consistency)
computation that was previously performed in notebookss./script The metric
quantifies how well per-cell velocity vectors align with principal curves that
connect successive cell clusters along a lineage graph.

The implementation follows these high-level steps:

1. Extract a 2D embedding and velocity representation from an :class:`AnnData`
   object (from ``.obsm`` or ``.layers``).
2. For each directed edge (``source`` → ``target`` cluster) fit a smooth
   parametric curve using a robust PCHIP interpolation over the projected
   manifold points belonging to the two clusters.
3. Project every cell along the edge onto the fitted curve and compute the
   cosine similarity between its velocity vector and the local curve tangent.
4. Aggregate cosine scores per edge, per cell (taking the best outgoing branch
   when multiple exist), and overall.

The public entry-point :func:`compute_batc` exposes configuration knobs while
keeping the original behaviour intact. Results can optionally be written back
into ``adata.obs``/``adata.uns`` for downstream inspection.
"""

from __future__ import annotations

from typing import Dict, Sequence, Tuple, Union

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.spatial import cKDTree
try:  # Optional progress feedback – silently skip if tqdm is absent.
    from tqdm import tqdm
except Exception:  # pragma: no cover - dependency not guaranteed
    def tqdm(iterable, **kwargs):  # type: ignore
        return iterable

# from ..utils import log_verbose

__all__ = ["compute_batc"]

ArrayLike = np.ndarray


class ParametricCurve2D:
    """Piecewise-cubic parametric curve ``(x(s), y(s))`` defined on ``s ∈ [0, 1]``."""

    def __init__(self, parameters: ArrayLike, points: ArrayLike) -> None:
        s = np.asarray(parameters, dtype=float)
        pts = np.asarray(points, dtype=float)
        if s.ndim != 1 or pts.ndim != 2 or pts.shape[1] != 2:
            raise ValueError("ParametricCurve2D expects 1D parameters and Nx2 points.")
        self.fx = PchipInterpolator(s, pts[:, 0])
        self.fy = PchipInterpolator(s, pts[:, 1])

    def eval(self, u: ArrayLike) -> ArrayLike:
        u_arr = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
        return np.column_stack([self.fx(u_arr), self.fy(u_arr)])

    def deriv(self, u: ArrayLike, *, h: float = 1e-3) -> ArrayLike:
        """Return numerical first derivative using centred differences."""

        u_arr = np.asarray(u, dtype=float)
        u_minus = np.clip(u_arr - h, 0.0, 1.0)
        u_plus = np.clip(u_arr + h, 0.0, 1.0)
        xy_minus = self.eval(u_minus)
        xy_plus = self.eval(u_plus)
        denom = (u_plus - u_minus)[:, None]
        denom[denom == 0.0] = 1e-12
        return (xy_plus - xy_minus) / denom


def _curve_eval(curve: ParametricCurve2D, u: float | ArrayLike, *, derivative: bool) -> ArrayLike:
    if derivative:
        values = curve.deriv(np.atleast_1d(u))
    else:
        values = curve.eval(np.atleast_1d(u))
    if np.ndim(u) == 0:
        return values[0]
    return values


def _extract_embedding_and_velocity(
    adata,
    embedding_key: str,
    velocity_key: str,
) -> Tuple[ArrayLike, ArrayLike]:
    """Return ``(X, V)`` matrices from an :class:`AnnData` object.

    The logic mirrors the original script: first attempt direct keys, then fall
    back to ``X_<key>`` / ``velocity_<key>`` style lookups, and finally try the
    ``layers`` container.
    """

    if hasattr(adata, "obsm"):
        if embedding_key in adata.obsm and velocity_key in adata.obsm:
            return np.asarray(adata.obsm[embedding_key]), np.asarray(adata.obsm[velocity_key])

        if embedding_key.startswith("X_"):
            X_key = embedding_key
            base = embedding_key[1:]
        else:
            X_key = f"X_{embedding_key}"
            base = embedding_key

        V_key = velocity_key
        if base and f"{velocity_key}_{base}" in adata.obsm:
            V_key = f"{velocity_key}_{base}"
        elif base.startswith("_") and f"{velocity_key}{base}" in adata.obsm:
            V_key = f"{velocity_key}{base}"

        if X_key in adata.obsm and V_key in adata.obsm:
            return np.asarray(adata.obsm[X_key]), np.asarray(adata.obsm[V_key])

    if hasattr(adata, "layers"):
        layers = getattr(adata, "layers", {})
        if embedding_key in layers and velocity_key in layers:
            return np.asarray(layers[embedding_key]), np.asarray(layers[velocity_key])
        if embedding_key.startswith("X_"):
            base = embedding_key[1:]
        else:
            base = embedding_key
        alt_X = base
        alt_V = f"{velocity_key}_{base}" if base else velocity_key
        if alt_X in layers and alt_V in layers:
            return np.asarray(layers[alt_X]), np.asarray(layers[alt_V])

    raise KeyError(
        f"Could not find a matching embedding/velocity pair for keys '{embedding_key}' and '{velocity_key}'."
    )


def _fit_curve_pchip(
    points: ArrayLike,
    src_center: ArrayLike,
    tgt_center: ArrayLike,
    *,
    use_bins: bool,
    n_bins: int,
    min_per_bin: int,
    n_samples: int,
) -> Tuple[ParametricCurve2D, ArrayLike, ArrayLike, float]:
    """Fit a parametric curve through points belonging to source & target clusters."""

    P = np.asarray(points, dtype=float)
    direction = tgt_center - src_center
    direction_norm = np.linalg.norm(direction)
    unit_dir = direction / (direction_norm + 1e-12)

    projection = (P - src_center) @ unit_dir
    order = np.argsort(projection)
    P_sorted = P[order]
    proj_sorted = projection[order]

    if use_bins and len(P_sorted) >= n_bins:
        bin_edges = np.linspace(proj_sorted.min(), proj_sorted.max(), n_bins + 1)
        anchors = []
        for idx in range(n_bins):
            if idx < n_bins - 1:
                mask = (proj_sorted >= bin_edges[idx]) & (proj_sorted < bin_edges[idx + 1])
            else:
                mask = (proj_sorted >= bin_edges[idx]) & (proj_sorted <= bin_edges[idx + 1])
            if np.count_nonzero(mask) >= min_per_bin:
                anchors.append(P_sorted[mask].mean(axis=0))
        skeleton = np.vstack(anchors) if len(anchors) >= 3 else P_sorted
    else:
        skeleton = P_sorted

    if len(skeleton) < 2:
        skeleton = np.vstack([P_sorted[0], P_sorted[-1]])

    segment_lengths = np.linalg.norm(np.diff(skeleton, axis=0), axis=1)
    s = np.r_[0.0, np.cumsum(segment_lengths)]
    if s[-1] == 0.0:
        s[-1] = 1.0
    s /= s[-1]

    curve = ParametricCurve2D(s, skeleton)

    u_grid = np.linspace(0.0, 1.0, int(n_samples))
    samples = curve.eval(u_grid)
    tree = cKDTree(samples)
    _, nearest = tree.query(P, k=1)
    u_all = u_grid[nearest]

    start = curve.eval(0.0).ravel()
    end = curve.eval(1.0).ravel()
    tangent = end - start
    orient = 1.0
    if direction_norm > 0 and np.dot(tangent, direction) < 0:
        orient = -1.0

    return curve, u_all[order], order, orient


def _prepare_cluster_data(
    X: ArrayLike,
    labels: ArrayLike,
) -> Tuple[Dict[str, ArrayLike], Dict[str, ArrayLike]]:
    unique_labels = np.unique(labels)
    idx_by_cluster = {lbl: np.where(labels == lbl)[0] for lbl in unique_labels}
    centroids = {lbl: X[idx].mean(axis=0) if len(idx) else np.zeros(2) for lbl, idx in idx_by_cluster.items()}
    return idx_by_cluster, centroids


def compute_batc(
    adata,

    cluster_edges: Sequence[Tuple[str, str]],
    *,
    cluster_key: str = "clusters",
    embedding_key: str = "umap",
    velocity_key: str = "velocity",
    use_bins: bool = True,
    n_bins: int = 60,
    min_per_bin: int = 5,
    n_samples: int = 800,
    store_in_adata: bool = True,
    progress: bool = False,
    dataset: str="",
    return_details: bool = False,
) -> Union[float, Tuple[float, dict]]:


    """Compute the Branching-aware Trajectory Consistency (BATC) score.

    BATC evaluates RNA velocity fields on predefined lineage graphs by fitting
    smooth principal curves between successive clusters and comparing their
    tangents with per-cell velocity vectors.  For each directed edge
    :math:`A \rightarrow B`, a curve :math:`\gamma_{A\to B}(u)` is fitted on the
    2D embedding of cells belonging to :math:`A \cup B`.  Every cell :math:`c`
    on the edge is projected onto the curve and the cosine similarity between the
    local tangent :math:`\tau_c` and its velocity :math:`v_c` is recorded:

    . math::

        \mathrm{BATC}_{A\to B} = \frac{1}{|I_{A\cup B}|}
        \sum_{c \in I_{A\cup B}}
        \frac{v_c \cdot \tau_c}{\lVert v_c \rVert\,\lVert \tau_c \rVert}.

    To handle branching, for each source cluster :math:`A` with outgoing targets
    :math:`B_1, \ldots, B_m` we compute these cosine scores for every outgoing
    edge and, for each cell :math:`c \in I_A`, keep the best matching branch

    . math::

        b_c^{(A)} = \max_{B \in \mathrm{Out}(A)}
        \frac{v_c \cdot \tau_c^{(A\to B)}}{\lVert v_c \rVert\,\lVert \tau_c^{(A\to B)} \rVert}.

    The final dataset-level BATC is the cell-weighted mean of these branch-aware
    scores over all source clusters:

    . math::

        \mathrm{BATC}_{\mathrm{overall}} =
        \frac{1}{\sum_A |I_A|} \sum_A \sum_{c \in I_A} b_c^{(A)}.

    Zero-norm vectors are mapped to ``nan`` scores, and all averages use
    ``numpy.nanmean`` for robustness.

    Parameters
    ----------
    adata:
        Annotated data matrix containing the embedding and velocity arrays.
    cluster_edges:
        Directed edges describing permitted transitions between clusters. Nodes
        are cast to strings to match ``adata.obs[cluster_key]``.
    cluster_key:
        Column in ``adata.obs`` holding cluster labels. Defaults to ``'clusters'``.
    embedding_key:
        Key identifying the 2D embedding. Defaults to ``'umap'`` which is
        resolved to ``adata.obsm['X_umap']`` when present. Any other key is
        resolved in the same fashion (``X_<key>`` or direct entry).
    velocity_key:
        Key identifying the velocity representation (must align with the
        embedding). With the default ``'velocity'`` the function looks for
        ``adata.obsm['velocity_umap']`` when ``embedding_key='umap'``.
    use_bins, n_bins, min_per_bin:
        Controls for the skeletonisation step when fitting the principal curves.
    n_samples:
        Number of evaluation points used to project cells onto each curve.
    store_in_adata:
        If ``True`` (default) write per-cell and aggregate scores back to the
        ``adata.obs``/``adata.uns`` containers.
    progress:
        Whether to display a progress bar over edges (requires ``tqdm``).
    return_details:
        If ``True``, return the unchanged overall score together with the
        existing edge cosine and branch-max diagnostics used for audit output.

    Returns
    -------
    float or tuple
        The BATC overall score (best-per-cell cosine mean), or ``(score,
        details)`` when ``return_details`` is true. The diagnostics report
        pre-max source-cell candidate means and branch winner counts; they do
        not define a second BATC aggregate.
    """
    cluster_series = adata.obs[cluster_key].astype(str)
    X, V = _extract_embedding_and_velocity(adata, embedding_key, velocity_key)
    if X.shape != V.shape:
        raise ValueError(f"Embedding shape {X.shape} does not match velocity shape {V.shape}.")
    if X.shape[1] != 2:
        raise ValueError("BATC expects a 2D embedding.")

    clusters = cluster_series.to_numpy()
    idx_by_cluster, centroids_data = _prepare_cluster_data(X, clusters)
    edges = [(str(src), str(tgt)) for src, tgt in cluster_edges]

    edge_fits: Dict[Tuple[str, str], Dict[str, object]] = {}
    iterable = edges
    if progress:
        iterable = tqdm(iterable, desc="Fitting BATC curves")
    for src, tgt in iterable:
        src_idx = idx_by_cluster.get(src, np.empty(0, dtype=int))
        tgt_idx = idx_by_cluster.get(tgt, np.empty(0, dtype=int))
        if src_idx.size == 0 and tgt_idx.size == 0:
            continue
        points_idx = np.concatenate([src_idx, tgt_idx])
        curve, u_sorted, order, orient = _fit_curve_pchip(
            X[points_idx],
            centroids_data.get(src, X[points_idx].mean(axis=0)),
            centroids_data.get(tgt, X[points_idx].mean(axis=0)),
            use_bins=use_bins,
            n_bins=n_bins,
            min_per_bin=min_per_bin,
            n_samples=n_samples,
        )
        edge_fits[(src, tgt)] = {
            "curve": curve,
            "u": u_sorted,
            "sorted_idx": points_idx[order],
            "orient": orient,
        }
    edge_cosines: Dict[Tuple[str, str], Dict[int, float]] = {}
    for (src, tgt), fit in edge_fits.items():
        curve: ParametricCurve2D = fit["curve"]  # type: ignore[assignment]
        u_vals = fit["u"]
        sorted_idx = fit["sorted_idx"]
        orient = float(fit["orient"])
        cos_scores: Dict[int, float] = {}
        for pos, cell_idx in enumerate(sorted_idx):
            velocity = V[cell_idx]
            tangent = orient * _curve_eval(curve, u_vals[pos], derivative=True)
            if np.linalg.norm(tangent) == 0 or np.linalg.norm(velocity) == 0:
                cos_scores[int(cell_idx)] = np.nan
                continue
            cos = float(np.dot(tangent, velocity) / (np.linalg.norm(tangent) * np.linalg.norm(velocity)))
            cos_scores[int(cell_idx)] = cos
        edge_cosines[(src, tgt)] = cos_scores

    next_map: Dict[str, list[str]] = {}
    for src, tgt in edges:
        next_map.setdefault(src, []).append(tgt)

    best_cosine_per_cell = np.full(X.shape[0], np.nan, dtype=float)
    branch_winner_source_cells: Dict[Tuple[str, str], int] = {
        edge: 0 for edge in edges
    }
    candidate_edge_details = []

    for src, targets in next_map.items():
        src_indices = idx_by_cluster.get(src, np.empty(0, dtype=int))
        if src_indices.size == 0:
            continue
        stacks = []
        for tgt in targets:
            mapping = edge_cosines.get((src, tgt), {})
            scores = np.array([mapping.get(int(cell_idx), np.nan) for cell_idx in src_indices], dtype=float)
            stacks.append(scores)
            candidate_values = np.nan_to_num(scores, nan=0.0)
            candidate_edge_details.append(
                {
                    "source": src,
                    "target": tgt,
                    "candidate_mean_source_cells": float(np.mean(candidate_values)),
                    "source_cells": int(src_indices.size),
                    "candidate_nonfinite_source_cells": int(np.count_nonzero(~np.isfinite(scores))),
                }
            )
        if not stacks:
            continue
        matrix = np.vstack(stacks)
        matrix = np.nan_to_num(matrix, nan=0.0)
        best_idx = np.nanargmax(matrix, axis=0)
        row_indices = np.arange(matrix.shape[1])
        best_values = matrix[best_idx, row_indices]
        best_cosine_per_cell[src_indices] = best_values
        for target_index, tgt in enumerate(targets):
            branch_winner_source_cells[(src, tgt)] = int(
                np.count_nonzero(best_idx == target_index)
            )

    candidate_by_edge = {
        (item["source"], item["target"]): item for item in candidate_edge_details
    }
    edge_details = []
    for src, tgt in edges:
        item = dict(candidate_by_edge.get(
            (src, tgt),
            {
                "source": src,
                "target": tgt,
                "candidate_mean_source_cells": None,
                "source_cells": 0,
                "candidate_nonfinite_source_cells": 0,
            },
        ))
        item["branch_winner_source_cells"] = branch_winner_source_cells.get(
            (src, tgt), 0
        )
        edge_details.append(item)

    valid_mask = np.isfinite(best_cosine_per_cell)
    batc_overall = float(np.nanmean(best_cosine_per_cell[valid_mask])) if valid_mask.any() else np.nan

    if store_in_adata:
        adata.uns["batc_overall"] = batc_overall












    if return_details:
        return batc_overall, {"overall": batc_overall, "edges": edge_details}
    return batc_overall
