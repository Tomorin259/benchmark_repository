"""Protect own-moment recovery, the affine derivative, and CosSim masking."""

import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.preprocessing import MinMaxScaler

from evaluation.gene.coordinate_units import restore_velovi_coordinates, velocity_cosine


class CoordinateUnitTests(unittest.TestCase):
    def setUp(self):
        graph = sparse.csr_matrix(np.array([[0, 2, 0, 0], [1, 0, 3, 0],
                                           [0, 2, 0, 0], [0, 0, 1, 0]], dtype=np.float32))
        binary = (graph > 0).astype(np.float32).toarray()
        np.fill_diagonal(binary, 1)
        normalized = sparse.csr_matrix(binary / binary.sum(1)[:, None])
        s = np.array([[0, 2], [2, 1], [6, 8], [4, 5]], dtype=np.float32)
        u = np.array([[3, 8], [0, 3], [5, 2], [9, 0]], dtype=np.float32)
        self.moments = [normalized.dot(sparse.csr_matrix(x)).astype(np.float32).toarray() for x in (s, u)]
        self.adata = SimpleNamespace(
            obs_names=pd.Index(["a", "b", "c", "d"]), var_names=pd.Index(["g1", "g2"]),
            obsp={"connectivities": graph},
            layers={"spliced": s, "unspliced": u,
                    "Ms": MinMaxScaler().fit_transform(self.moments[0]),
                    "Mu": MinMaxScaler().fit_transform(self.moments[1])})

    def test_own_graph_forward_exact_and_velocity_derivative(self):
        before = {key: value.copy() for key, value in self.adata.layers.items()}
        ms, mu, inverse_s, inverse_u, evidence = restore_velovi_coordinates(self.adata, [0, 1])
        for restored, expected, inverse in zip((ms, mu), self.moments, (inverse_s, inverse_u)):
            np.testing.assert_array_equal(restored, expected)
            scaler = MinMaxScaler().fit(expected)
            np.testing.assert_array_equal(inverse, 1.0 / scaler.scale_.astype(np.float64))
            # Translation enters state coordinates only; the derivative uses 1/scale.
            np.testing.assert_array_equal(np.zeros(2) * inverse, np.zeros(2))
        self.assertTrue(all(row["forward_exact"] for row in evidence))
        for key, original in before.items():
            np.testing.assert_array_equal(self.adata.layers[key], original)

    def test_nonexact_forward_and_constant_axes_fail(self):
        self.adata.layers["Ms"][0, 0] = np.nextafter(self.adata.layers["Ms"][0, 0], np.float32(1))
        with self.assertRaisesRegex(ValueError, "exactly reproduce"):
            restore_velovi_coordinates(self.adata, [0, 1])
        self.setUp()
        self.adata.layers["spliced"][:, 0] = 1
        self.adata.layers["Ms"][:, 0] = 0
        with self.assertRaisesRegex(ValueError, "exactly reproduce"):
            restore_velovi_coordinates(self.adata, [0, 1])

    def test_actual_preprocessed_input_aligns_and_checks_abundance(self):
        s, u = self.adata.layers["spliced"], self.adata.layers["unspliced"]
        moments = [s * 2.5 + 7.25, u * 0.5 - 3.5]
        for axis, moment in zip(("s", "u"), moments):
            self.adata.layers["M" + axis] = MinMaxScaler().fit_transform(moment)
        rows, columns = [2, 0, 3, 1], [1, 0]
        source = SimpleNamespace(obs_names=self.adata.obs_names[rows],
                                 var_names=self.adata.var_names[columns], layers={})
        for key, value in zip(("spliced", "unspliced", "Ms", "Mu"), (s, u, *moments)):
            source.layers[key] = value[rows][:, columns].copy()
        ms, mu, _, _, evidence = restore_velovi_coordinates(self.adata, [0, 1], source)
        np.testing.assert_array_equal(ms, moments[0])
        np.testing.assert_array_equal(mu, moments[1])
        self.assertTrue(all("actual_preprocessed_input" in row["moment_source_kind"] for row in evidence))
        self.assertTrue(any(row["saved_graph_forward_error_before_source_resolution"] > 0 for row in evidence))
        source.layers["spliced"][0, 0] += 1
        with self.assertRaisesRegex(ValueError, "abundance differs"):
            restore_velovi_coordinates(self.adata, [0, 1], source)

    def test_cossim_retains_historical_norm_mask(self):
        velocity = np.array([[1., 0.], [-1., 0.], [1e-7, 0.], [0., 0.], [np.nan, 1.]])
        truth = np.tile([1., 0.], (5, 1))
        value, count = velocity_cosine(velocity, truth)
        self.assertEqual((value, count), (0., 2))
        self.assertEqual(velocity_cosine(np.zeros((2, 2)), np.ones((2, 2))), (0., 0))


if __name__ == "__main__":
    unittest.main()
