import csv
import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts.tsv"
FIELDS = ["artifact_id", "category", "path", "sha256", "status", "uri", "notes"]


def read_artifacts():
    with ARTIFACTS.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == FIELDS
        return list(reader)


def test_artifact_registry_schema_and_paths():
    rows = read_artifacts()
    assert rows
    ids = [row["artifact_id"] for row in rows]
    assert len(ids) == len(set(ids))
    for row in rows:
        assert re.fullmatch(r"[a-z0-9][a-z0-9._-]*", row["artifact_id"])
        assert row["status"] in {"available", "external", "pending"}
        assert not row["path"].startswith(("/", "mnt/", "Users/"))
        assert not re.search(r"/(?:mnt|Users)/|(?:^|\s)(?:/mnt|/Users)(?:/|\s|$)", row["path"])
        if row["status"] == "available":
            path = ROOT / row["path"]
            assert path.is_file(), row["path"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
        else:
            assert row["path"] == ""
        if row["status"] == "external":
            assert "Zenodo URI" in row["notes"]


def test_registry_has_no_deprecated_variant_and_formal_ensemble_is_available():
    text = ARTIFACTS.read_text()
    assert "".join(("center", "weighted")) not in text.lower()
    rows = read_artifacts()
    pending = [row for row in rows if row["status"] == "pending"]
    assert not pending
    formal = {
        row["artifact_id"]: row
        for row in rows
        if row["artifact_id"] in {
            "ensemble-mean-4-formal-all14",
            "ensemble-row-geometric-median-4-formal-all14",
        }
    }
    assert set(formal) == {
        "ensemble-mean-4-formal-all14",
        "ensemble-row-geometric-median-4-formal-all14",
    }
    assert all(row["status"] == "available" for row in formal.values())
    assert len({row["path"] for row in formal.values()}) == 1
    assert all("smoke" not in row["notes"].lower() for row in rows)


def test_assembly_manifest_materialized_panels_are_verified():
    manifest = ROOT / "figures/assembly/assembly_manifest.tsv"
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 124
    materialized = [row for row in rows if row["status"] == "materialized"]
    assert len(materialized) == 117
    panels = [row for row in rows if row["figure"] == "suppfig10"]
    assert len(panels) == 2
    assert all(row["status"] == "materialized" for row in panels)
    assert {row["panel"] for row in panels} == {"A", "B"}
    assert not any("workflow schematic" in row["role"] for row in panels)
    expected = {
        "suppfig10/A_row_geometric_median4_cbdir.svg",
        "suppfig10/B_row_geometric_median4_batc.svg",
    }
    assert {row["materialized_path"] for row in panels} == expected


def test_online_supplementary_figure_gaps_are_explicit():
    with (ROOT / "figures/assembly/assembly_manifest.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for figure in ("suppfig09", "suppfig11"):
        gaps = [row for row in rows if row["figure"] == figure]
        missing = [row for row in gaps if row["status"] == "missing_source"]
        assert len(missing) == 1
        if figure == "suppfig11":
            panels = [row for row in rows if row["figure"] == figure and row["status"] == "materialized"]
            assert {row["panel"] for row in panels} == {"A", "B"}
            for row in panels:
                source = ROOT / row["source_path"]
                target = ROOT / "figures/assembly" / row["materialized_path"]
                assert source.read_bytes() == target.read_bytes()
                assert hashlib.sha256(target.read_bytes()).hexdigest() == row["sha256"]
