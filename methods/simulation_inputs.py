"""Simulation inputs shared by the nine method launchers.

Paths come from the public input-file manifest and are relative to the
velocity artifact root. Dataset IDs come from data/datasets.tsv.
"""

import csv
from pathlib import Path


def simulation_inputs():
    manifest = Path(__file__).resolve().parents[1] / "data/input_file_manifest.tsv"
    with manifest.open() as handle:
        return {
            row["dataset_id"]: {
                "input": row["path"], "label_key": "bin",
                "embedding_key": None, "preprocessed": False, "simulation": True,
            }
            for row in csv.DictReader(handle, delimiter="\t")
            if row["role"] == "fixed_simulation"
        }
