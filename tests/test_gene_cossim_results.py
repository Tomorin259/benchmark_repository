import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv"
PROVENANCE = ROOT / "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.provenance.json"
EXPECTED_SHA256 = "11fa92fec20a73ebc2105ded3825d0de3cdd5f8f6651b876a89f6dcaa99d23c9"
HISTORICAL_SOURCE_SHA256 = "25f216ca4573556ae484d41c275ddd5bc2dfd2c19b5e987f272b8aa315a97662"


def test_simulation_cossim_table_provenance_and_hash():
    digest = hashlib.sha256(TABLE.read_bytes()).hexdigest()
    assert digest == EXPECTED_SHA256
    provenance = json.loads(PROVENANCE.read_text())
    assert provenance["direct_source_locator"] == "artifact://velocity/plot/to_use_sim/details/results_new/cbdir_cossim_scores_long.csv"
    assert provenance["direct_source_sha256"] == HISTORICAL_SOURCE_SHA256
    assert provenance["direct_source_sha256"] != digest
    assert "historical" in provenance["direct_source_role"]
    assert provenance["table_sha256"] == digest
    restoration = provenance["unit_restoration"]
    assert (restoration["updated_method"], restoration["updated_field"], restoration["updated_rows"]) == ("veloVI", "CosSim", 70)
    path = ROOT / restoration["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == restoration["sha256"]
    evidence = json.loads(path.read_text())
    assert evidence["historical_baseline"]["table_sha256"][str(TABLE.relative_to(ROOT))] == HISTORICAL_SOURCE_SHA256
    for name, expected in evidence["evidence"].items():
        assert hashlib.sha256((path.parent / name).read_bytes()).hexdigest() == expected


def test_simulation_cossim_table_schema_and_complete_coverage():
    with TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == [
            "Dataset", "Gene", "Method_raw", "Method", "CBDir", "CosSim"
        ]
        rows = list(reader)

    assert len(rows) == 630
    assert Counter(row["Dataset"] for row in rows) == {
        "branching": 180, "murk": 90, "all": 360
    }
    assert Counter(row["Method"] for row in rows) == {
        method: 70 for method in {
            "cellDancer", "DeepVelo", "DeepVelo_scv", "Graphvelo_scv",
            "MultiVelo", "scVelo_dyn", "TIvelo", "UniTVelo", "veloVI"
        }
    }
    keys = {(row["Dataset"], row["Method"], row["Gene"]) for row in rows}
    assert len(keys) == len(rows)
    assert all(math.isfinite(float(row["CBDir"])) for row in rows)
    assert all(math.isfinite(float(row["CosSim"])) for row in rows)


def test_velovi_cossim_matches_paired_restoration_without_mask_changes():
    unit = ROOT / "evaluation/gene/results/unit_restoration"
    with (unit / "velovi_score_pairs.csv").open(newline="") as handle:
        paired = {(row["dataset"], row["gene"]): row for row in csv.DictReader(handle)}
    with TABLE.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["Method"] == "veloVI"]
    for row in rows:
        source = paired[row["Dataset"], row["Gene"]]
        assert float(row["CosSim"]) == float(source["restored_cossim"])
        assert source["native_cossim_cells"] == source["restored_cossim_cells"]
        assert int(source["cossim_mask_changes"]) == 0
