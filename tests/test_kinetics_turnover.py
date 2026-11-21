"""Check the checked-in S11 results against both supported Figure 6D input layouts."""

import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


MODULE = Path(__file__).resolve().parents[1] / "evaluation/kinetics_recovery"
RESULTS = MODULE / "results/turnover"
PLOT_DATA = MODULE / "figures/fig6/plot_data"


@pytest.mark.parametrize("released_filenames", [False, True])
def test_turnover_reproduces_s11(tmp_path, released_filenames):
    pytest.importorskip("numpy")
    command = [sys.executable, str(MODULE / "turnover.py")]
    if released_filenames:
        input_dir = tmp_path / "Figure_6"
        input_dir.mkdir()
        for source in PLOT_DATA.glob("rpe1_*_kinetic_class*.csv"):
            shutil.copyfile(source, input_dir / f"D_{source.name}")
        script = tmp_path / "turnover.py"
        shutil.copyfile(command[1], script)
        command = [sys.executable, str(script), "--input-dir", str(input_dir)]
    output = tmp_path / "replay"
    subprocess.run([*command, "--output", str(output)], check=True, cwd=tmp_path)
    canonical = RESULTS / "class_f1.csv"
    assert (output / canonical.name).read_bytes() == canonical.read_bytes()
    with canonical.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [(row["class"], int(row["reference_genes"])) for row in rows] == [
        ("Cooperative", 94), ("Neutral", 142), ("Destabilizing", 14), ("Macro", 250)
    ]
    assert rows[-1]["permutation_p"] == rows[-1]["holm_p"] == ""
    manifest = json.loads((output / "manifest.json").read_text())
    expected = json.loads((RESULTS / "manifest.json").read_text())
    assert manifest.pop("numpy")
    del expected["numpy"]
    assert manifest == expected
    assert manifest["output_sha256"][canonical.name] == hashlib.sha256(canonical.read_bytes()).hexdigest()
