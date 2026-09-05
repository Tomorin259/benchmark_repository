#!/usr/bin/env python3
"""Small behavior tests for the cell-level evaluation boundary."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from evaluation.cell.metrics import consistency_scores, resolve_artifact


class CellEvaluationTests(unittest.TestCase):
    def test_consistency_uses_neighbor_cosines_and_skips_self(self):
        velocity = np.asarray(
            [[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float
        )
        neighbors = [[0, 1, 2], [1, 0, 2], [2, 0, 1]]
        scores = consistency_scores(velocity, neighbors, n_neighbors=2)
        self.assertTrue(np.allclose(scores, [0.5, 0.5, 0.0]))

    def test_resolve_artifact_rejects_absolute_and_parent_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "result.h5ad").write_bytes(b"test")
            self.assertEqual(
                resolve_artifact(root, "result.h5ad"), root / "result.h5ad"
            )
            with self.assertRaises(ValueError):
                resolve_artifact(root, "/tmp/result.h5ad")
            with self.assertRaises(ValueError):
                resolve_artifact(root, "../result.h5ad")


if __name__ == "__main__":
    unittest.main()
