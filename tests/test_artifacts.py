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


def test_registry_has_no_rowgeom_and_formal_ensemble_is_available():
    text = ARTIFACTS.read_text()
    assert "rowgeom" not in text.lower()
    rows = read_artifacts()
    pending = [row for row in rows if row["status"] == "pending"]
    assert not pending
    formal = {
        row["artifact_id"]: row
        for row in rows
        if row["artifact_id"] in {
            "ensemble-mean-4-formal-all14",
            "ensemble-centerweighted-4-formal-all14",
        }
    }
    assert set(formal) == {
        "ensemble-mean-4-formal-all14",
        "ensemble-centerweighted-4-formal-all14",
    }
    assert all(row["status"] == "available" for row in formal.values())
    assert len({row["path"] for row in formal.values()}) == 1
    assert all("smoke" not in row["notes"].lower() for row in rows)


def test_assembly_manifest_materialized_panels_are_verified():
    manifest = ROOT / "figures/assembly/assembly_manifest.tsv"
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 112
    materialized = [row for row in rows if row["status"] == "materialized"]
