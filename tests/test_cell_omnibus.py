"""Protect Table S10's panel-wise correction and the BATC sensitivity calculation."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from evaluation.cell.omnibus import CELL, CELL_METRICS, MATCHED, MATCHED_METRICS, batc_ranks, main, omnibus, validate_identities


class CellOmnibusTests(unittest.TestCase):
    def test_committed_panels(self):
        out = CELL / "results/statistics"
        for path, metrics, shape, name in (
            (CELL / "results/cell_metrics.csv", CELL_METRICS, (14, 20), "s10a"),
            (MATCHED, MATCHED_METRICS, (7, 9), "s10b"),
        ):
            data = pd.read_csv(path, float_precision="round_trip")
            expected = pd.read_csv(out / f"{name}.csv", float_precision="round_trip")
            actual = omnibus(data, metrics, shape)
            self.assertEqual(list(actual.metric), list(expected.metric))
            np.testing.assert_allclose(actual.iloc[:, 1:], expected.iloc[:, 1:], atol=0, rtol=1e-11)
            with self.assertRaises(ValueError):
                omnibus(data.iloc[:-1], metrics, shape)

    def test_tied_ranks_and_equal_dataset_means(self):
        data = pd.DataFrame({"dataset": list("aabbcc"), "method": ["x", "y"] * 3,
                             "batc": [0.1, 0.1, 0.9, 0.2, 0.2, 0.4]})
        ranks, omissions = batc_ranks(data)
        self.assertEqual(list(ranks.loc[ranks.dataset == "a", "rank"]), [1.5, 1.5])
        after_b = omissions.loc[omissions.omitted_dataset == "b"]
        np.testing.assert_allclose(after_b.mean_batc, [0.15, 0.25])
        self.assertEqual(list(after_b["rank"]), [2, 1])


    def test_degenerate_rank_patterns(self):
        values = pd.DataFrame({"dataset": list("aaabbb"), "method": list("xyzxyz"),
                               "flat": [1, 1, 1, 2, 2, 2], "same": [1, 1, 2, 3, 3, 4]})
        result = omnibus(values, ["flat", "same"], (2, 3))
        np.testing.assert_equal(result.friedman_q.to_numpy(), [0, 4])
        np.testing.assert_equal(result.kendall_w.to_numpy(), [0, 1])
        np.testing.assert_equal(result.iman_davenport_f.to_numpy(), [0, np.inf])
        np.testing.assert_equal(result.p_holm.to_numpy(), [1, 0])

    def test_input_identities(self):
        cell = pd.read_csv(CELL / "results/cell_metrics.csv")
        matched = pd.read_csv(MATCHED)
        for table, field, value in [(cell, "dataset", None), (cell, "method", "unknown"),
                                    (matched, "config_id", "C99"), (matched, "cell_score_route", "wrong")]:
            changed = table.copy()
            changed.loc[0, field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_identities(changed if table is cell else cell, changed if table is matched else matched)
        with self.assertRaises(ValueError):
            validate_identities(pd.concat([cell, cell.iloc[:1]]), matched)

    def test_publish_only_complete_new_output(self):
        with tempfile.TemporaryDirectory() as root:
            output = Path(root) / "result"
            with patch.object(sys, "argv", ["omnibus", "--output", str(output)]):
                with patch.object(pd.DataFrame, "to_csv", side_effect=OSError("disk full")), self.assertRaises(OSError):
                    main()
                self.assertEqual(list(Path(root).iterdir()), [])
                main()
                original = {p.name: p.read_bytes() for p in output.iterdir()}
                self.assertEqual(len(original), 4)
                with self.assertRaises(SystemExit):
                    main()
                self.assertEqual({p.name: p.read_bytes() for p in output.iterdir()}, original)
