#!/usr/bin/env python3
"""Add the benchmark's four Gene layers to a fitted UniTVelo mode-1 output.

This is the retained get_velocity_unspliced.py calculation, with explicit
input/output paths. The saved spliced velocity and fitted parameters are kept.
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse


def add_velocity_u(adata):
    """Differentiate the fitted Gaussian curve using the retained conventions."""
    required = {"fit_a", "fit_h", "fit_t", "fit_offset", "fit_gamma",
                "fit_beta", "fit_intercept", "velocity_genes"}
    missing = required - set(adata.var)
    if missing or "fit_t" not in adata.layers:
        raise ValueError(f"UniTVelo fitted parameters required: var={sorted(missing)}, "
                         f"layers/fit_t missing={'fit_t' not in adata.layers}")
    if "velocity_u" in adata.layers:
        raise ValueError("velocity_u already exists; select the original mode-1 output")

    values = adata.layers["fit_t"]
    t = np.asarray(values.toarray() if sparse.issparse(values) else values,
                   dtype=np.float64)
    # Mode 1 shares cell time. The historical export masks unfitted columns;
    # its reconstruction fills those entries before applying the gene mask.
    finite_columns = np.flatnonzero(np.isfinite(t).all(axis=0))
    if not np.isfinite(t).all() and finite_columns.size:
        reference_t = t[:, finite_columns[0]][:, None]
        missing_time = ~np.isfinite(t)
        t = t.copy()
        t[missing_time] = np.broadcast_to(reference_t, t.shape)[missing_time]

    a, h, t0, offset, gamma, beta, intercept = (
        np.asarray(adata.var[key], dtype=np.float64)[None, :]
        for key in ("fit_a", "fit_h", "fit_t", "fit_offset", "fit_gamma",
                    "fit_beta", "fit_intercept")
    )
    beta = beta.copy()
    beta[np.isclose(beta, 0.0)] = np.nan
    dt = t - t0
    fit_s = h * np.exp(-a * dt ** 2) + offset
    ds_dt = (fit_s - offset) * (-2.0 * a * dt)
    d2s_dt2 = (fit_s - offset) * (4.0 * a ** 2 * dt ** 2 - 2.0 * a)
    fit_u = (ds_dt + gamma * fit_s) / beta + intercept
    du_dt = (d2s_dt2 + gamma * ds_dt) / beta

    selected = np.asarray(adata.var["velocity_genes"], dtype=bool)
    for key, values in (("velocity_u", du_dt), ("velocity_s_reconstructed", ds_dt),
                        ("fit_s_reconstructed", fit_s), ("fit_u_reconstructed", fit_u)):
        values[:, ~selected] = np.nan
        adata.layers[key] = values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output exists: {args.output}; choose a new file")
    adata = ad.read_h5ad(args.input)
    add_velocity_u(adata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Wrote {args.output}: added four Gene layers to {adata.shape}")


if __name__ == "__main__":
    main()
