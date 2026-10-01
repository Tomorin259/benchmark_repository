import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from evaluation.gene.refresh_panel_manifest import FIELDS as GENE_FIELDS, refresh_manifest


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
    for row in materialized:
        source = (ROOT / row["source_path"]).read_bytes()
        target = (ROOT / "figures/assembly" / row["materialized_path"]).read_bytes()
        assert source == target, row["materialized_path"]
        assert hashlib.sha256(source).hexdigest() == row["sha256"], row["source_path"]
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


def test_gene_panel_manifest_hashes_and_supporting_scope():
    manifest = ROOT / "evaluation/gene/figures/panels/assembly_manifest.tsv"
    original = manifest.read_bytes()
    assert refresh_manifest(manifest) == 0
    assert manifest.read_bytes() == original
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    summaries = [row for row in rows if row["panel"] == "Figure4_D"]
    assert {row["dataset"] for row in summaries} == {"murk", "branching", "all"}
    assert all(row["source_kind"] == "supporting_gene_summary" for row in summaries)
    provenance = json.loads((manifest.parent / "provenance.json").read_text())
    assert provenance["record_scope"] == "historical_generation"
    assert "panel_count" not in provenance


def test_manifest_refresh_preserves_references_and_requires_explicit_write(tmp_path):
    panel = tmp_path / "panel.png"
    panel.write_bytes(b"unchanged panel")
    provenance = tmp_path / "provenance.json"
    provenance.write_text('{"historical": true}\n')
    rows = [{"panel": figure, "dataset": "murk", "asset": "gene",
             "path": "panel.png", "score_source": "recorded/source.csv",
             "sha256": "0" * 64, "source_kind": "score_annotated_panel"}
            for figure in ("Figure4_B", "SuppFig1")]
    manifest = tmp_path / "assembly_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GENE_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    original = manifest.read_bytes()
    with pytest.raises(ValueError, match="Stale panel hashes"):
        refresh_manifest(manifest)
    assert manifest.read_bytes() == original
    script = ROOT / "evaluation/gene/refresh_panel_manifest.py"
    subprocess.run([sys.executable, "-S", str(script), "--manifest", str(manifest), "--write"],
                   check=True, capture_output=True, text=True)
    with manifest.open(newline="") as handle:
        refreshed = list(csv.DictReader(handle, delimiter="\t"))
    expected = hashlib.sha256(panel.read_bytes()).hexdigest()
    for previous, current in zip(rows, refreshed):
        assert current == {**previous, "sha256": expected}
    updated = manifest.read_bytes()
    assert refresh_manifest(manifest, write=True) == 0
    assert manifest.read_bytes() == updated
    assert panel.read_bytes() == b"unchanged panel"
    assert provenance.read_text() == '{"historical": true}\n'


@pytest.mark.parametrize("path", ["missing.png", "../outside.png", "/outside.png"])
def test_manifest_refresh_rejects_invalid_paths_without_writing(tmp_path, path):
    manifest = tmp_path / "assembly_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GENE_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow({"panel": "Figure4_B", "dataset": "murk", "asset": "gene",
                         "path": path, "score_source": "source.csv", "sha256": "0" * 64,
                         "source_kind": "score_annotated_panel"})
    original = manifest.read_bytes()
    with pytest.raises(ValueError):
        refresh_manifest(manifest, write=True)
    assert manifest.read_bytes() == original


def test_figure4d_registration_and_composite_match_current_source():
    with (ROOT / "figures/assembly/assembly_manifest.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    formal = [row for row in rows if row["figure"] == "fig4" and row["panel"] == "D"
              and row["role"] == "paired cell and gene directionality comparison"]
    assert len(formal) == 1
    source = "evaluation/cell/simulation/figures/simulation_cross_level_cbdir.png"
    assert formal[0]["source_path"] == source
    assert formal[0]["materialized_path"] == "fig4/D_summary/simulation_cross_level_cbdir.png"
    assert not (ROOT / "figures/assembly/fig4/D_summary/summary_dot_interval.png").exists()
    supporting = [row for row in rows if row["figure"] == "fig4" and row["panel"] == "D"
                  and row["role"] == "supporting gene-only scenario summary"]
    assert len(supporting) == 3
    composite = json.loads((ROOT / "figures/fig4_composite_manifest.json").read_text())
    assert composite["paper_figure"] == "fig4"
    assert composite["composite_sha256"] == hashlib.sha256(
        (ROOT / composite["composite_path"]).read_bytes()).hexdigest()
    assert composite["panels"][0]["source_path"] == source
    assert composite["panels"][0]["source_sha256"] == formal[0]["sha256"]
    for item in composite["metric_sources"]:
        assert hashlib.sha256((ROOT / item["path"]).read_bytes()).hexdigest() == item["sha256"]
