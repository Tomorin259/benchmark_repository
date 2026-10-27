import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv"
PROVENANCE = ROOT / "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.provenance.json"
EXPECTED_SHA256 = "4e2aaa1c6809a15aba648c0fd2501a5e967d5aacbe4f8fbe7657c8f5af578697"


def test_simulation_cossim_table_provenance_and_hash():
    digest = hashlib.sha256(TABLE.read_bytes()).hexdigest()
    assert digest == EXPECTED_SHA256
    provenance = json.loads(PROVENANCE.read_text())
    assert provenance["table_sha256"] == digest
    truth = ROOT / provenance["velocity_reference"]["path"]
    assert hashlib.sha256(truth.read_bytes()).hexdigest() == provenance["velocity_reference"]["sha256"]
    truth_manifest = json.loads(truth.read_text())
    assert hashlib.sha256((ROOT / "evaluation/gene/simulation/make_truth_reference.py").read_bytes()).hexdigest() == truth_manifest["generator_sha256"]
    assert hashlib.sha256((ROOT / provenance["reproduction"]["script"]).read_bytes()).hexdigest() == provenance["reproduction"]["script_sha256"]
    for record in truth_manifest["scenario_records"]:
        asset = truth.parent / record["output"]
        assert hashlib.sha256(asset.read_bytes()).hexdigest() == record["output_sha256"]
    units = provenance["coordinate_units"]
    path = ROOT / units["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == units["sha256"]
    evidence = json.loads(path.read_text())
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


def test_velovi_unit_restoration_preserves_historical_cell_masks():
    unit = ROOT / "evaluation/gene/results/unit_restoration"
    with (unit / "velovi_score_pairs.csv").open(newline="") as handle:
        paired = {(row["dataset"], row["gene"]): row for row in csv.DictReader(handle)}
    with TABLE.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["Method"] == "veloVI"]
    assert len(rows) == 70
    assert any(float(row["CosSim"]) != float(paired[row["Dataset"], row["Gene"]]["restored_cossim"])
               for row in rows)
    for row in rows:
        source = paired[row["Dataset"], row["Gene"]]
        assert source["native_cossim_cells"] == source["restored_cossim_cells"]
        assert int(source["cossim_mask_changes"]) == 0


def test_simulation_replay_scores_other_methods_with_raw_velocity_masks(tmp_path, monkeypatch):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from evaluation.gene import run_reference30_gene_cbdir as replay

    data = ad.AnnData(np.ones((4, 1)),
                     obs=pd.DataFrame({"bin": ["A", "A", "B", "B"]},
                                      index=["a", "b", "c", "d"]),
                     var=pd.DataFrame(index=["g"]))
    for key, values in {"Ms": [0., 1., 2., 3.], "Mu": [1., 2., 3., 4.],
                        "velocity_s_true": [1.] * 4, "velocity_u_true": [0.] * 4}.items():
        data.layers[key] = np.array(values)[:, None]
    truth = tmp_path / "truth.npz"
    # Truth rows deliberately differ from method order. Embedded H5AD truth
    # points the opposite way, so the assertion also checks source selection.
    data.layers["velocity_s_true"] = -np.ones((4, 1))
    np.savez(truth, cell_ids=np.asarray(["d", "b", "a", "c"]),
             gene_names=np.asarray(["g"]), velocity_s_true=np.array([-1., 1., 1., 1.])[:, None],
             velocity_u_true=np.zeros((4, 1)))
    reference = tmp_path / "neighbors.npz"
    np.savez(reference, obs_names=np.asarray(data.obs_names).astype(str),
             indices=np.tile(np.arange(4), (4, 1)))
    monkeypatch.setattr(replay, "ARTIFACT_ROOT", tmp_path)
    methods = [method for method in replay.METHOD_ORDER if method != "veloVI"]
    paths = {}
    for method in methods:
        s_key, u_key = replay.VELOCITY_KEYS[method]
        data.layers[s_key] = np.array([1., 0., np.nan, -1.])[:, None]
        data.layers[u_key] = np.zeros((4, 1))
        paths[method] = f"{method}.h5ad"
        data.write_h5ad(tmp_path / paths[method])
    *_, cosines = replay.evaluate_dataset(
        "synthetic", {"genes": ["g"], "edges": [("A", "B")]}, paths, {},
        tmp_path, reference, simulation=True, truth_path=truth, methods=methods)
    assert len(cosines) == 8
    assert {row["method"] for row in cosines} == set(methods)
    for row in cosines:
        assert row["cossim"] == 1.
        assert row["evaluable_cells"] == 2
