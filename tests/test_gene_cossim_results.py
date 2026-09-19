import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv"
PROVENANCE = ROOT / "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.provenance.json"
EXPECTED_SHA256 = "25f216ca4573556ae484d41c275ddd5bc2dfd2c19b5e987f272b8aa315a97662"


def test_simulation_cossim_table_provenance_and_hash():
    digest = hashlib.sha256(TABLE.read_bytes()).hexdigest()
    assert digest == EXPECTED_SHA256
    provenance = json.loads(PROVENANCE.read_text())
    assert provenance["direct_source_locator"] == "artifact://velocity/plot/to_use_sim/details/results_new/cbdir_cossim_scores_long.csv"
    assert provenance["direct_source_sha256"] == EXPECTED_SHA256
    assert provenance["table_sha256"] == digest


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
