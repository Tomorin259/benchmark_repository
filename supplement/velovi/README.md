# veloVI coordinate sensitivity

Supplementary Table S6 is a post-hoc coordinate diagnostic on 14 Cell datasets.
It compares three coordinate conditions using stored Cell-role velocities.

## Inspect and reproduce

| File | Contents |
| --- | --- |
| [results/full_dataset_three_condition.csv](results/full_dataset_three_condition.csv) | Full-precision 14 × 3 conditions |
| [results/table_s6.csv](results/table_s6.csv) | Wide table, C − B and unweighted macro-average |
| [results/manifest.json](results/manifest.json) | Frozen labels/edges, input hashes, settings and output identities |

Regenerate the table without scientific dependencies:

```bash
python3 supplement/velovi/run_coordinate_sensitivity.py \
  --table-only --output-root /path/to/table-export
```

Recalculate the diagnostic with the recorded H5ADs and reference moments:

```bash
python supplement/velovi/run_coordinate_sensitivity.py \
  --artifact-root /path/to/velocity-artifacts \
  --reference-moments-root /path/to/reference-moments \
  --output-root /path/to/fresh-coordinate-results
```

Materialize the H5ADs and reference moments identified in the manifest.
The script verifies input and TIvelo metric-source hashes before scoring.
Use Python 3.9.23, scVelo 0.3.1, TIvelo 0.1.4, NumPy 1.26.4, SciPy 1.13.1,
AnnData 0.10.8 and pandas 2.3.1 for the checked reproduction environment.
Only scVelo was recorded in the original log; the other versions were inspected
in the later reproduction environment, not independently archived originally.

## Conditions

Cells, common genes, fixed-neighbor candidates, labels and transitions stay
fixed:

| Condition | State | Velocity |
| --- | --- | --- |
| A | Native stored `Ms` | Stored velocity |
| B | Reference `Ms` | Same stored velocity |
| C | Reference `Ms` | Gene-wise transformed `V'_ig = a_g V_ig` |

OLS fits `M_ref,ig = a_g M_native,ig + b_g` per gene/dataset without transition
labels or scores. The intercept is not applied to velocity.
The fitted map describes gene-wise coordinate scaling and can also reflect
preprocessing and moment smoothing.
Original latent-time scaling stays in all conditions.
Constant native genes use slope 1; constant reference genes retain the archived
R-squared convention of 1. Match genes case-insensitively and cells by ID.

All conditions recompute velocity graphs with scVelo 0.3.1, 30 neighbors,
one recursive expansion and `mode_neighbors="distances"`.
TIvelo's positive-graph TransCosine averages finite source-cell values per
annotated edge, then across evaluable edges; it is not a normalized transition
probability.

The retained macro-average is A=0.1147, B=0.0778, C=0.2193, C−B=0.1415.
TransCosine increases in all 14 datasets after gene-wise velocity transformation
relative to condition B. CBDir-HD values are reported in the full table.

```bash
python3 -m unittest discover -s tests -p test_velovi_coordinate_sensitivity.py -v
```

Checks cover the 42 rows, source hashes, B/primary agreement (within 2e-8),
table regeneration and shuffled-cell affine transformation.
