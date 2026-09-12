# Fixed-neighbor ensemble provenance

This record distinguishes the three commits that are otherwise easy to
confuse in the fixed-neighbor release:

| Field | Value | Meaning |
| --- | --- | --- |
| `lineage/current_upstream` | `transition_ensemble/main@4f2a621f5c1814635f33ec70f75c9bba2e7fd23f` | Current upstream source lineage used for the migration. |
| `release_index_commit` | `d9fae0bb8e67fd49b719746914bc06645c9dd738` | Commit recorded when the upstream all-14 release index was generated. |
| `numerical_reproduction_commit` | `baa79856522f2528ec798a5a0641574280a90a58` | Clean source commit used for the numerical reproduction. |

The formal release is the path-free record in
`transition-ensemble:fixed-neighbors-all14:schema6:2026-09-12`. Its scores are
the two signed-transition variants over all 14 datasets:

```text
Mean-4:               CBDir 0.349882266127113, BATC 0.48541164067602616
RowGeometricMedian-4: CBDir 0.37668835002093143, BATC 0.5144182874429892
```

The Pancreas and Forebrain migration checks compared the benchmark runner with
the formal generated outputs. `ensemble_transition_signed` and
`ensemble_velocity_umap` were equal entry by entry for both variants. The
derived CBDir/BATC values differed only at approximately 1e-8 floating-point
precision.

The 20-method comparison continues to use the repository's fixed-reference-30
table, `evaluation/cell/results/cell_metrics.csv`. The public fixed-neighbor
NPZ files are retained unchanged: their numerical arrays are elementwise equal
to the audited upstream arrays, while their metadata is portable and therefore
has benchmark-specific SHA256 values.

No server absolute path is part of this provenance contract. Large H5AD files
remain external artifacts; Git contains the compact score tables, selected
stream panels, manifests, and this provenance record.
