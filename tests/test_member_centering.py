"""Scientific checks for neighbor-weight centering."""

import unittest

import numpy as np
from scipy import sparse

from ensemble.benchmark_member_centering import center_neighbor_weights
from ensemble.reproject import project_graphvelo


class MemberCenteringTests(unittest.TestCase):
    def test_uses_each_members_nonzero_offdiagonal_support(self):
        original = sparse.csr_matrix([
            [.7, .5, -.25, 0.], [0., .6, .2, 0.],
            [0., 0., 1., 0.], [.1, .1, 0., .8],
        ])
        before = original.toarray().copy()
        actual = center_neighbor_weights(original)
        np.testing.assert_allclose(actual.toarray(), [
            [0., .375, -.375, 0.], [0., 0., 0., 0.],
            [0., 0., 0., 0.], [0., 0., 0., 0.],
        ], atol=1e-15)
        np.testing.assert_array_equal(original.toarray(), before)

    def test_projection_matches_explicit_mean_geometry_subtraction(self):
        original = sparse.csr_matrix([
            [.2, .3, -.1], [.4, .2, .1], [.1, .2, .7],
        ])
        coordinates = np.asarray([[0., 0.], [2., 1.], [-1., 3.]])
        actual = project_graphvelo(center_neighbor_weights(original), coordinates)
        expected = np.zeros_like(coordinates)
        for cell in range(3):
            neighbors = np.flatnonzero(original.toarray()[cell])
            neighbors = neighbors[neighbors != cell]
            weights = original.toarray()[cell, neighbors]
            displacement = coordinates[neighbors] - coordinates[cell]
            expected[cell] = weights @ displacement - weights.mean() * displacement.sum(axis=0)
        np.testing.assert_allclose(actual, expected, atol=1e-15)


if __name__ == "__main__":
    unittest.main()
