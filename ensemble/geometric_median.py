"""Row-wise Euclidean geometric median for four signed transition operators.

The four input methods are treated symmetrically.  For every source cell, this
module finds the geometric median of the four complete signed off-diagonal
transition rows.  The result remains in their convex hull and is not
renormalized after aggregation.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.optimize import minimize

N_INPUT_ROWS = 4
RELATIVE_TOLERANCE = 1e-8
ABSOLUTE_TOLERANCE = 1e-14
MAX_ITERATIONS = 1024


def _distribution_summary(values):
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return {
            "min": 0.0,
            "p10": 0.0,
            "median": 0.0,
            "mean": 0.0,
            "p90": 0.0,
            "max": 0.0,
        }
    return {
        "min": float(values.min()),
        "p10": float(np.quantile(values, 0.10)),
        "median": float(np.median(values)),
        "mean": float(values.mean()),
        "p90": float(np.quantile(values, 0.90)),
        "max": float(values.max()),
    }


def _row_gram_matrices(template, method_values):
    """Compute one 4 x 4 Gram matrix for every sparse source-cell row."""
    template = sparse.csr_matrix(template)
    values = np.asarray(method_values, dtype=np.float64)
    if template.shape[0] != template.shape[1]:
        raise ValueError("row aggregation requires a square template")
    if values.ndim != 2 or values.shape != (N_INPUT_ROWS, template.nnz):
        raise ValueError(
            "row geometric median requires a 4 x template.nnz value array"
        )
    if not np.isfinite(values).all():
        raise ValueError("method_values contains NaN or infinity")

    n_obs = template.shape[0]
    row_ids = np.repeat(np.arange(n_obs), np.diff(template.indptr))
    if np.any(row_ids == template.indices):
        raise ValueError("row geometric median expects off-diagonal support")

    gram = np.zeros((n_obs, N_INPUT_ROWS, N_INPUT_ROWS), dtype=np.float64)
    for left in range(N_INPUT_ROWS):
        for right in range(left, N_INPUT_ROWS):
            products = values[left] * values[right]
            row_dot = np.bincount(
                row_ids,
                weights=products,
                minlength=n_obs,
            )
            gram[:, left, right] = row_dot
            gram[:, right, left] = row_dot
    return gram, row_ids


def _distances_from_row_coefficients(gram, coefficients):
    """Return distances from a convex row combination to all four rows."""
    gram_times_coefficients = np.einsum("nij,nj->ni", gram, coefficients)
    center_squared_norm = np.einsum(
        "ni,ni->n", coefficients, gram_times_coefficients
    )
    squared = (
        np.diagonal(gram, axis1=1, axis2=2)
        - 2.0 * gram_times_coefficients
        + center_squared_norm[:, None]
    )
    return np.sqrt(np.maximum(squared, 0.0))


def _slsqp_geometric_median_coefficient_fallback(row_gram, initial):
    """Finish one slow four-point row on the simplex with a convex solver."""
    row_gram = np.asarray(row_gram, dtype=np.float64)
    diagonal = np.diag(row_gram)

    def distances(coefficients):
        gram_times = row_gram.dot(coefficients)
        center_squared = float(coefficients.dot(gram_times))
        squared = diagonal - 2.0 * gram_times + center_squared
        return np.sqrt(np.maximum(squared, 0.0)), gram_times

    def objective(coefficients):
        return float(distances(coefficients)[0].sum())

    def gradient(coefficients):
        row_distances, gram_times = distances(coefficients)
        inverse = 1.0 / np.maximum(row_distances, 1e-15)
        return inverse.sum() * gram_times - row_gram.T.dot(inverse)

    result = minimize(
        objective,
        np.asarray(initial, dtype=np.float64),
        jac=gradient,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * N_INPUT_ROWS,
        constraints={
            "type": "eq",
            "fun": lambda coefficients: float(coefficients.sum() - 1.0),
            "jac": lambda coefficients: np.ones(N_INPUT_ROWS, dtype=np.float64),
        },
        options={"ftol": 1e-14, "maxiter": 1000, "disp": False},
    )
    if not result.success or not np.isfinite(result.x).all():
        raise RuntimeError(
            "SLSQP geometric-median fallback failed: {}".format(result.message)
        )

    coefficients = np.clip(np.asarray(result.x, dtype=np.float64), 0.0, 1.0)
    coefficients /= coefficients.sum()
    if objective(coefficients) > objective(initial) + 1e-10:
        raise RuntimeError("SLSQP geometric-median fallback increased objective")
    return coefficients, int(result.nit)


def _row_geometric_median_coefficients(
    gram,
    *,
    relative_tolerance=RELATIVE_TOLERANCE,
    absolute_tolerance=ABSOLUTE_TOLERANCE,
    max_iterations=MAX_ITERATIONS,
):
    """Solve all four-point row medians with modified Weiszfeld updates.

    Every iterate is represented as four convex coefficients.  Distances and
    updates therefore use only each row's 4 x 4 Gram matrix, rather than
    repeatedly traversing its potentially large neighbor support.
    """
    gram = np.asarray(gram, dtype=np.float64)
    if gram.ndim != 3 or gram.shape[1:] != (N_INPUT_ROWS, N_INPUT_ROWS):
        raise ValueError("gram must have shape (n_obs, 4, 4)")

    n_obs = gram.shape[0]
    coefficients = np.full(
        (n_obs, N_INPUT_ROWS), 1.0 / N_INPUT_ROWS, dtype=np.float64
    )
    diagonal = np.diagonal(gram, axis1=1, axis2=2)
    pairwise_squared = (
        diagonal[:, :, None] + diagonal[:, None, :] - 2.0 * gram
    )
    pairwise_distances = np.sqrt(np.maximum(pairwise_squared, 0.0))
    row_scale = pairwise_distances.max(axis=(1, 2))
    identical = row_scale <= absolute_tolerance
    singular_tolerance = absolute_tolerance + relative_tolerance * row_scale

    # A geometric median can coincide exactly with an input row.  Detect that
    # case from the convex subgradient condition before Weiszfeld iteration;
    # otherwise convergence toward such a point can be arbitrarily slow.
    input_row_is_median = np.zeros((n_obs, N_INPUT_ROWS), dtype=bool)
    for candidate in range(N_INPUT_ROWS):
        distances = pairwise_distances[:, candidate, :]
        coincident = distances <= singular_tolerance[:, None]
        noncoincident = ~coincident
        inverse = np.divide(
            1.0,
            distances,
            out=np.zeros_like(distances),
            where=noncoincident,
        )
        residual_coefficients = inverse.copy()
        residual_coefficients[:, candidate] -= inverse.sum(axis=1)
        residual_norm = np.sqrt(
            np.maximum(
                np.einsum(
                    "ni,nij,nj->n",
                    residual_coefficients,
                    gram,
                    residual_coefficients,
                ),
                0.0,
            )
        )
        multiplicity = coincident.sum(axis=1)
        input_row_is_median[:, candidate] = (
            residual_norm <= multiplicity + relative_tolerance
        )

    has_input_row_median = input_row_is_median.any(axis=1)
    if np.any(has_input_row_median):
        selected = input_row_is_median[has_input_row_median].astype(np.float64)
        coefficients[has_input_row_median] = (
            selected / selected.sum(axis=1, keepdims=True)
        )
    converged = identical | has_input_row_median
    iterations = np.zeros(n_obs, dtype=np.int32)

    for iteration in range(1, max_iterations + 1):
        active = np.flatnonzero(~converged)
        if active.size == 0:
            break

        active_gram = gram[active]
        current = coefficients[active]
        distances = _distances_from_row_coefficients(active_gram, current)
        singular_tolerance = (
            absolute_tolerance + relative_tolerance * row_scale[active]
        )
        near = distances <= singular_tolerance[:, None]
        has_near = near.any(axis=1)
        updated = current.copy()
        stationary = np.zeros(active.size, dtype=bool)

        regular = ~has_near
        if np.any(regular):
            inverse = 1.0 / distances[regular]
            updated[regular] = inverse / inverse.sum(axis=1, keepdims=True)

        singular = np.flatnonzero(has_near)
        if singular.size:
            singular_near = near[singular]
            non_near = ~singular_near
            inverse = np.divide(
                1.0,
                distances[singular],
                out=np.zeros_like(distances[singular]),
                where=non_near,
            )
            inverse_sum = inverse.sum(axis=1)
            multiplicity = singular_near.sum(axis=1).astype(np.float64)

            residual_coefficients = (
                inverse - current[singular] * inverse_sum[:, None]
            )
            residual_norm = np.sqrt(
                np.maximum(
                    np.einsum(
                        "ni,nij,nj->n",
                        residual_coefficients,
                        active_gram[singular],
                        residual_coefficients,
                    ),
                    0.0,
                )
            )
            singular_stationary = (inverse_sum == 0.0) | (
                residual_norm <= multiplicity + relative_tolerance
            )
            stationary[singular] = singular_stationary

            moving = ~singular_stationary
            if np.any(moving):
                beta = inverse[moving] / inverse_sum[moving, None]
                eta = multiplicity[moving] / residual_norm[moving]
                updated[singular[moving]] = (
                    eta[:, None] * current[singular[moving]]
                    + (1.0 - eta)[:, None] * beta
                )

        delta = updated - current
        movement = np.sqrt(
            np.maximum(
                np.einsum("ni,nij,nj->n", delta, active_gram, delta),
                0.0,
            )
        )
        movement_tolerance = (
            absolute_tolerance + relative_tolerance * row_scale[active]
        )
        newly_converged = stationary | (movement <= movement_tolerance)
        coefficients[active] = updated
        iterations[active] = iteration
        converged[active[newly_converged]] = True

    fallback_iterations = np.zeros(n_obs, dtype=np.int32)
    failed = np.flatnonzero(~converged)
    for obs_id in failed:
        coefficients[obs_id], fallback_iterations[obs_id] = (
            _slsqp_geometric_median_coefficient_fallback(
                gram[obs_id], coefficients[obs_id]
            )
        )
        converged[obs_id] = True

    # Guard against harmless roundoff before using the coefficients to prove
    # convex-hull and signed-L1 invariants downstream.
    coefficients[np.abs(coefficients) < absolute_tolerance] = 0.0
    if np.any(coefficients < -absolute_tolerance):
        raise AssertionError("geometric-median coefficients left the convex hull")
    coefficients = np.maximum(coefficients, 0.0)
    coefficients /= coefficients.sum(axis=1, keepdims=True)
    return (
        coefficients,
        iterations,
        identical,
        has_input_row_median,
        fallback_iterations,
        row_scale,
    )


def row_geometric_median4(
    template,
    method_values,
    *,
    return_cell_diagnostics=False,
):
    """Aggregate each source cell's four complete signed-T rows robustly.

    The output is the Euclidean geometric median of the four native
    off-diagonal rows.  It is not L1-renormalized: loss of directional mass is
    retained as an interpretable consequence of row-level disagreement and is
    completed to the diagonal by ``operator_from_offdiagonal``.
    """
    values = np.asarray(method_values, dtype=np.float64)
    gram, row_ids = _row_gram_matrices(template, values)
    (
        coefficients,
        iterations,
        identical,
        input_row_median,
        fallback_iterations,
        row_scale,
    ) = _row_geometric_median_coefficients(gram)

    aggregate = np.zeros(values.shape[1], dtype=np.float64)
    for method_id in range(N_INPUT_ROWS):
        aggregate += values[method_id] * coefficients[row_ids, method_id]

    mean4 = values.mean(axis=0)
    difference = aggregate - mean4
    difference_l2 = np.sqrt(
        np.bincount(
            row_ids,
            weights=difference * difference,
            minlength=template.shape[0],
        )
    )
    final_distances = _distances_from_row_coefficients(gram, coefficients)
    mean_coefficients = np.full_like(coefficients, 1.0 / N_INPUT_ROWS)
    mean_distances = _distances_from_row_coefficients(gram, mean_coefficients)
    final_objective = final_distances.sum(axis=1)
    mean_objective = mean_distances.sum(axis=1)
    objective_increase = final_objective - mean_objective
    objective_tolerance = RELATIVE_TOLERANCE * np.maximum(1.0, mean_objective)
    if np.any(objective_increase > objective_tolerance):
        raise AssertionError(
            "row geometric median has a larger objective than Mean-4"
        )

    diagnostics = {
        "solver": "modified Weiszfeld in the four-row Gram geometry",
        "relative_tolerance": RELATIVE_TOLERANCE,
        "absolute_tolerance": ABSOLUTE_TOLERANCE,
        "max_iterations": MAX_ITERATIONS,
        "n_rows": int(template.shape[0]),
        "n_identical_input_rows": int(identical.sum()),
        "n_medians_at_input_rows": int(input_row_median.sum()),
        "n_slsqp_fallback_rows": int(np.count_nonzero(fallback_iterations)),
        "max_slsqp_fallback_iterations": int(
            fallback_iterations.max(initial=0)
        ),
        "n_converged_rows": int(template.shape[0]),
        "max_iterations_used": int(iterations.max(initial=0)),
        "mean_iterations_used": float(iterations.mean()),
        "max_objective_increase_vs_mean4": float(
            objective_increase.max(initial=0.0)
        ),
        "mean_objective_reduction_vs_mean4": float(
            np.mean(mean_objective - final_objective)
        ),
        "row_l2_difference_vs_mean4": _distribution_summary(difference_l2),
        "row_input_diameter": _distribution_summary(row_scale),
        "mean_implicit_convex_coefficient_by_method": [
            float(value) for value in coefficients.mean(axis=0)
        ],
        "max_coefficient_sum_error": float(
            np.max(np.abs(coefficients.sum(axis=1) - 1.0), initial=0.0)
        ),
        "min_implicit_convex_coefficient": float(
            coefficients.min() if coefficients.size else 0.0
        ),
    }
    if return_cell_diagnostics:
        cell_diagnostics = {
            "gram": gram,
            "row_ids": row_ids,
            "coefficients": coefficients,
            "row_input_diameter": row_scale,
        }
        return aggregate, diagnostics, cell_diagnostics
    return aggregate, diagnostics


def validate_row_geometric_median4():
    """Check method symmetry, sign symmetry, and basic robust-center cases."""
    rows = np.repeat(np.arange(N_INPUT_ROWS), N_INPUT_ROWS - 1)
    columns = np.asarray([1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2])
    template = sparse.csr_matrix(
        (np.ones(rows.size), (rows, columns)),
        shape=(N_INPUT_ROWS, N_INPUT_ROWS),
    )
    values = np.asarray(
        [
            [0.2, -0.1, 0.3, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, -0.2, 0.1, 0.0],
            [0.2, -0.1, 0.3, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, -0.2, 0.1, 0.0],
            [0.2, -0.1, 0.3, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, -0.2, 0.1, 0.0],
            [0.2, -0.1, 0.3, 0.9, 0.0, 0.0, -0.1, 0.0, 0.0, -0.2, 0.1, 0.0],
        ],
        dtype=np.float64,
    )
    base, base_diagnostics = row_geometric_median4(template, values)
    permuted, _ = row_geometric_median4(template, values[[2, 0, 3, 1]])
    reversed_sign, _ = row_geometric_median4(template, -values)

    checks = {
        "permutation_invariant": bool(
            np.allclose(base, permuted, atol=1e-10)
        ),
        "sign_equivariant": bool(
            np.allclose(reversed_sign, -base, atol=1e-10)
        ),
        "identical_rows_are_preserved": bool(
            np.allclose(base[:3], values[0, :3], atol=1e-10)
        ),
        "three_zero_votes_reject_one_outlier": bool(
            np.allclose(base[3:6], 0.0, atol=1e-8)
        ),
        "geometric_objective_not_worse_than_mean": bool(
            base_diagnostics["max_objective_increase_vs_mean4"] <= 1e-8
        ),
    }
    if not all(checks.values()):
        raise AssertionError(
            "RowGeometricMedian-4 invariant failed: {}".format(checks)
        )
    return checks
