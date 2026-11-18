"""Checks for aligning signed neighbor weights before ensemble aggregation."""
import unittest
import numpy as np
from scipy import sparse
from ensemble.benchmark_normalization import normalize_offdiagonal, aggregate_rows
from ensemble.reproject import project_graphvelo

class MemberNormalizationTests(unittest.TestCase):
    def test_signed_mass_diagonal_and_zero_rows(self):
        matrix = sparse.csr_matrix([[8., 2., -1.], [0., 9., 0.], [1e-14, 0., 1.]])
        before = matrix.toarray().copy()
        result = normalize_offdiagonal(matrix)
        np.testing.assert_allclose(result.toarray(), [[0., 2/3, -1/3], [0., 0., 0.], [0., 0., 0.]], atol=1e-15)
        np.testing.assert_array_equal(matrix.toarray(), before)

    def test_positive_rescaling_and_self_weights_do_not_change_normalized_operator(self):
        base = sparse.csr_matrix([[0., .3, -.1], [-.2, 0., .1], [.1, .4, 0.]])
        scaled = base.multiply(np.array([.01, 10., 4.])[:, None]).tocsr()
        scaled.setdiag([10., -3., 7.])
        np.testing.assert_allclose(normalize_offdiagonal(base).toarray(),
                                   normalize_offdiagonal(scaled).toarray(), atol=1e-15)
        for count in (4, 5):
            result, _, weights = aggregate_rows([normalize_offdiagonal(scaled)] * count)
            np.testing.assert_allclose(result.toarray(), normalize_offdiagonal(base).toarray(), atol=1e-14)
            self.assertEqual(weights.shape, (3, count))
        coords = np.array([[0., 0.], [1., 2.], [4., 1.]])
        expected = project_graphvelo(base, coords) / np.asarray(abs(base).sum(1))
        np.testing.assert_allclose(project_graphvelo(normalize_offdiagonal(base), coords),
                                   expected, atol=1e-14)

if __name__ == "__main__":
    unittest.main()
