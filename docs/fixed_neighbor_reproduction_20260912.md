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

## BATC NaN-policy reproduction

The repository BATC scorer follows GRAVITY's metric construction but applies
`np.nan_to_num(..., nan=0.0)` to the matrix of candidate branch scores before
selecting the per-source-cell maximum. The upstream GRAVITY selector instead
uses `np.nanargmax` directly. The repository rule is therefore an explicit
benchmark robustness adaptation: it defines an all-NaN source-cell case as
zero, while the upstream implementation raises `ValueError: All-NaN slice
encountered`. For a mixture of NaN and finite negative branches, the two rules
can also differ in principle because the zero fallback may become the maximum.

An isolated shadow rerun at benchmark commit `72d6754` changed only this branch
selection rule and retained the formal input manifest, fixed-reference-30
assets, dataset edges, BATC bin counts, embeddings, and velocity projections.
The Cell result was:

| Outcome | Count | Evidence |
| --- | ---: | --- |
| Upstream policy completed; BATC exactly matched the repository | 278/280 | Zero differing pair scores; maximum absolute difference 0 |
| Upstream policy raised on an all-NaN source-cell column | 2/280 | `bonemarrow/TSvelo`; `rpe1_gen/TFvelo` |

The repository fallback-defined BATC values for the two upstream-undefined
pairs are `0.31053362099524107` for `bonemarrow/TSvelo` and
`-0.12657291592050282` for `rpe1_gen/TFvelo`. These values are benchmark
results, not values emitted by the unmodified upstream selector.

The same comparison was run directly on the 28 formal fixed-neighbor Ensemble
H5AD results (14 datasets × Mean-4 and RowGeometricMedian-4). All 28 upstream
policy evaluations completed; every value matched the current scorer on the
same H5AD exactly, with maximum absolute difference 0. The approximately 1e-8
archived-H5AD reproduction differences noted above are independent floating-
point rerun effects, not effects of the BATC NaN policy.

Decision: retain the repository fallback and the formal Cell and Ensemble
numbers. No score table, rank, panel, or manuscript value changes as a result
of this policy audit. Public descriptions must identify the zero fallback as a
benchmark adaptation rather than claiming line-for-line identity with
GRAVITY's branch selector.

The 20-method comparison continues to use the repository's fixed-reference-30
table, `evaluation/cell/results/cell_metrics.csv`. The public fixed-neighbor
NPZ files are retained unchanged: their numerical arrays are elementwise equal
to the audited upstream arrays, while their metadata is portable and therefore
has benchmark-specific SHA256 values.

No server absolute path is part of this provenance contract. Large H5AD files
remain external artifacts; Git contains the compact score tables, selected
stream panels, manifests, and this provenance record.
