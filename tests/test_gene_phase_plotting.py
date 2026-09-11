#!/usr/bin/env python3
"""Synthetic tests for the public gene-phase plotting contract."""

from __future__ import annotations

import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from evaluation.gene.plotting.phase_portrait import (
    _sample_phase_cells,
    _scale_phase_vectors,
    plot_phase_portrait,
)


class GenePhasePlottingTests(unittest.TestCase):
    def test_log_scaling_preserves_direction_and_caps_length(self):
        velocities = np.vstack(
            [np.tile([[1.0, 0.0]], (99, 1)), [[1000.0, 1000.0]]]
        )
        display, diagnostics = _scale_phase_vectors(
            velocities,
            (10.0, 10.0),
            mode="log1p",
            reference_quantile=0.90,
            reference_length=0.05,
            max_length=0.10,
        )

        for source, rendered in zip(velocities, display):
            if np.linalg.norm(source) == 0:
                continue
            self.assertAlmostEqual(
                source[0] * rendered[1] - source[1] * rendered[0],
                0.0,
                places=12,
            )

        relative_lengths = np.linalg.norm(display / np.asarray([10.0, 10.0]), axis=1)
        self.assertLessEqual(float(relative_lengths.max()), 0.10 + 1e-12)
        self.assertEqual(diagnostics["n_capped"], 1)

    def test_phase_sampling_is_stable_under_row_reordering(self):
        coordinates = np.asarray(
            [
                [0.0, 0.0],
                [0.2, 0.8],
                [0.4, 0.1],
                [0.6, 0.9],
                [0.8, 0.2],
                [1.0, 1.0],
                [0.1, 0.5],
                [0.9, 0.6],
            ]
        )
        cell_ids = np.asarray([f"cell-{i}" for i in range(len(coordinates))])
        first = _sample_phase_cells(coordinates, cell_ids, grid=(5, 5))
        permutation = np.asarray([5, 2, 7, 0, 4, 1, 6, 3])
        second = _sample_phase_cells(
            coordinates[permutation],
            cell_ids[permutation],
            grid=(5, 5),
        )

        self.assertEqual(
            set(cell_ids[first].tolist()),
            set(cell_ids[permutation][second].tolist()),
        )

    def test_plot_does_not_modify_source_layers(self):
        n_obs = 16
        adata = AnnData(
            X=np.zeros((n_obs, 1)),
            obs=pd.DataFrame(
                {"clusters": pd.Categorical(["a", "b"] * (n_obs // 2))}
            ),
            var=pd.DataFrame(index=["g1"]),
        )
        adata.obs_names = [f"cell-{i:02d}" for i in range(n_obs)]
        s = np.linspace(0.0, 4.0, n_obs)[:, None]
        u = np.linspace(0.1, 1.1, n_obs)[:, None]
        vs = np.linspace(-2.0, 2.0, n_obs)[:, None]
        vu = np.linspace(2.0, -2.0, n_obs)[:, None]
        adata.layers["Ms"] = s
        adata.layers["Mu"] = u
        adata.layers["velocity"] = vs
        adata.layers["velocity_u"] = vu
        before = {
            key: value.copy()
            for key, value in adata.layers.items()
        }

        fig, ax = plt.subplots()
        _, diagnostics = plot_phase_portrait(
            adata,
            "g1",
            ax=ax,
            color_key="clusters",
            arrow_grid=(4, 4),
            legend_loc=None,
            return_diagnostics=True,
        )
        self.assertEqual(diagnostics["plot_version"], "gene-phase-log1p-005-010-v1")
        self.assertLessEqual(
            diagnostics["display_magnitude_max"],
            diagnostics["arrow_max_length"] + 1e-12,
        )
        for key, value in before.items():
            np.testing.assert_array_equal(adata.layers[key], value)
        plt.close(fig)

    def test_missing_unspliced_velocity_is_rejected(self):
        adata = AnnData(
            X=np.zeros((2, 1)),
            var=pd.DataFrame(index=["g1"]),
        )
        adata.layers["Ms"] = np.asarray([[0.0], [1.0]])
        adata.layers["Mu"] = np.asarray([[0.0], [1.0]])
        adata.layers["velocity"] = np.asarray([[0.0], [1.0]])
        with self.assertRaises(ValueError):
            plot_phase_portrait(adata, "g1", layer_v_u=None)


if __name__ == "__main__":
    unittest.main()

