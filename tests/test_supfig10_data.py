"""Protect the S10 score table contract and the historical renderer's output boundary."""

from pathlib import Path
import sys

import pandas as pd
import pytest

from ensemble import plot_supfig10
from ensemble.core import load_config
from ensemble.supfig10_data import write_panel_data

ROOT = Path(__file__).resolve().parents[1]


def test_supfig10_source_tables_use_full5_and_sort_each_metric(tmp_path):
    datasets = list(load_config()["datasets"])
    selected = tmp_path / "Table_S5_B.csv"
    pd.DataFrame({
        "dataset": ["macro_average", *reversed(datasets)],
        "display_name": ["Macro-average", *reversed(datasets)],
        "cbdir": 0.9,
        "batc": -0.9,
    }).to_csv(selected, index=False)
    output = tmp_path / "panels"
    write_panel_data(selected, output)

    bootstrap = pd.read_csv(output / "panel_a_data.csv")
    dots = pd.read_csv(output / "plot_data.csv")
    macro = pd.read_csv(output / "macro_means.csv")
    assert (len(bootstrap), len(dots), len(macro)) == (20, 588, 21)
    assert macro.loc[macro.method.eq("GeoMedian-5"), ["cbdir", "batc"]].values.tolist() == [[0.9, -0.9]]
    assert bootstrap.loc[bootstrap.original_mean_front, "method"].tolist() == ["veloVI", "TIvelo"]
    assert set(dots.loc[dots.constituent, "method_label"]) == {
        "DeepVelo", "GraphVelo_scv", "TIvelo", "UniTVelo", "veloVI",
    }
    ensemble_dots = dots.loc[dots.method.eq("GeoMedian-5")]
    assert ensemble_dots.plot_row.tolist() == [0] * 14 + [20] * 14
    assert ensemble_dots.dataset.tolist() == datasets * 2


def test_legacy_cli_requires_explicit_input_and_refuses_existing_or_assembly_output(tmp_path, monkeypatch):
    with pytest.raises(SystemExit, match="2"):
        plot_supfig10.build_parser().parse_args(["--output-dir", str(tmp_path / "explicit")])
    protected = ROOT / "figures/assembly/suppfig10/unused-test-output"
    for destination in (tmp_path, protected):
        monkeypatch.setattr(sys, "argv", ["plot_supfig10", "--ensemble-csv", "unused.csv",
                                        "--output-dir", str(destination)])
        with pytest.raises(SystemExit, match="2"):
            plot_supfig10.main()
