# veloVI coordinate-sensitivity diagnostic

This directory publishes the post-hoc coordinate diagnostic behind
Supplementary Table S6. It does not retrain veloVI and it does not replace the
primary cell benchmark.

## Conditions

For each of the 14 cell datasets, the diagnostic keeps the scored cells,
common genes, reference-neighbor candidates, labels, and transition edges
fixed:

* A: native stored Ms with stored velocity;
* B: reference Ms with the unchanged stored velocity;
* C: reference Ms with gene-wise transformed velocity
  V'_ig = a_g V_ig.

The coefficients are fitted independently for each dataset and gene, without
using transition labels or benchmark scores, by
ordinary least squares,

    M_ref,ig = a_g M_native,ig + b_g,

and the intercept is not applied to velocity. Thus a_g is a diagnostic
coordinate map, not an official veloVI inverse-normalization factor. The
mapping can reflect preprocessing and moment smoothing in addition to
gene-specific scale differences.

The original stored gene-wise latent-time scaling is retained in all three
conditions. For a constant native gene, the archived calculation uses a slope
of 1; constant reference genes receive the archived R-squared convention of 1.
Common genes are matched case-insensitively. Cells are aligned by their IDs.
TransCosine is calculated using TIvelo's positive scVelo graph metric: finite
source-cell values are averaged per annotated edge, then across evaluable
edges. This is not a row-normalized transition probability.

All three conditions recompute the velocity graph with the same scVelo 0.3.1 settings:
30 neighbors, one recursive neighbor expansion, and mode_neighbors="distances".
The stored cell-role layers["velocity"] output is used; the gene-role
velocity_u output is not part of Table S6.

## Reproduce

Large H5AD files and reference moments are external artifacts. The repository
contains the fixed reference-neighbor assets and records artifact hashes in
artifacts.tsv. The frozen manifest additionally records all 14 input paths,
labels, transition edges, and checksums, so later benchmark-configuration
changes cannot silently redefine this diagnostic. Supply the artifact root and the external reference-moment
directory explicitly:

    python supplement/velovi/run_coordinate_sensitivity.py \
        --artifact-root /path/to/velocity-artifacts \
        --reference-moments-root /path/to/reference_moments_union_scvelo025_20260907 \
        --output-root supplement/velovi/results

The script writes:

* full_dataset_three_condition.csv: full-precision rows for all 14 × 3
  conditions;
* table_s6.csv: a wide, publication-facing view with C - B and the
  unweighted macro-average;
* manifest.json: software versions, graph settings, input hashes, script
  hash, and output hashes.

Recomputation verifies the frozen input checksums and the TIvelo metric source
hash before scoring. The verified environment uses Python 3.9.23, scVelo
0.3.1, TIvelo 0.1.4, NumPy 1.26.4, SciPy 1.13.1, AnnData 0.10.8, and pandas
2.3.1. The original log recorded only the scVelo version; the other versions
were inspected in the reproduction environment on 2026-10-02 and must not be
read as independently archived original-run metadata.

The result CSV and Table S6 can be inspected without the large input files.
Public download URIs for those inputs remain unpublished in the repository's
artifact index; a clean checkout alone cannot run the full recomputation.
The input SHA-256 values identify the required external files.

The committed CSV is the archived run used for Table S6. A clean rerun should
reproduce its rounded values and the following macro-average:

    A = 0.1147, B = 0.0778, C = 0.2193, C - B = 0.1415

The CSV also retains CBDir-HD for transparency; Table S6's 14/14 improvement
claim concerns TransCosine, not CBDir-HD. This post-hoc result does not establish
better biological inference or isolate the causal effect of training-time
normalization.

Regenerate the full-precision wide table without scientific dependencies:

    python3 supplement/velovi/run_coordinate_sensitivity.py \
        --table-only --output-root /path/to/table-export

Run the lightweight regression checks:

    python3 -m unittest discover -s tests -p test_velovi_coordinate_sensitivity.py -v

These check the 42-row contract, hashes, agreement of B with the primary
benchmark (within 2e-8), table regeneration, and a synthetic shuffled-cell
affine transformation. Full recomputation requires the external artifacts.

The uniform-scaling and time-scaling controls, and the separate
reference-moment retraining experiments, are intentionally kept outside this
module because they answer different questions.
