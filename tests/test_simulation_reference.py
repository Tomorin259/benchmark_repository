"""Reference identity contracts without native H5ADs."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from evaluation.gene.simulation.reference import load_velocity_reference


class SimulationReferenceTests(unittest.TestCase):
    def test_named_alignment_and_subsets(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            np.savez(path / "case_truth.npz", cell_ids=["b", "a"], gene_names=["y", "x"],
                     velocity_s_true=[[1., 2.], [3., 4.]], velocity_u_true=[[5., 6.], [7., 8.]])
            s, u = load_velocity_reference("case", ["a", "b"], ["x"], path)
            np.testing.assert_array_equal(s, [[4.], [2.]])
            np.testing.assert_array_equal(u, [[8.], [6.]])
            with self.assertRaisesRegex(ValueError, "Nonunique"):
                load_velocity_reference("case", ["a", "a"], ["x"], path)
            with self.assertRaisesRegex(ValueError, "does not cover"):
                load_velocity_reference("case", ["a"], ["missing"], path)

    def test_invalid_reference_fails_before_scoring(self):
        for cells, genes, s, u in [
            (["a", "a"], ["x"], [[1.], [2.]], [[1.], [2.]]),
            (["a"], ["x"], [[np.nan]], [[1.]]),
            (["a"], ["x"], [[1., 2.]], [[1.]]),
        ]:
            with self.subTest(cells=cells, s=s), tempfile.TemporaryDirectory() as directory:
                path = Path(directory)
                np.savez(path / "case_truth.npz", cell_ids=cells, gene_names=genes,
                         velocity_s_true=s, velocity_u_true=u)
                with self.assertRaises(ValueError):
                    load_velocity_reference("case", ["a"], ["x"], path)


if __name__ == "__main__":
    unittest.main()
