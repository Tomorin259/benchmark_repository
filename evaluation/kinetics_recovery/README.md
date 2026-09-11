# Kinetics-recovery

This module evaluates method-estimated kinetic parameters against the scEU-seq
reference tables for `rpe1` and `organoid` systems. It consumes saved method-result
H5AD files; it does not prepare datasets or train methods.

The configuration keys and result filenames retain the internal system IDs
`rpe1` and `organoid`. Compact result tables use the display names
`rpe1` and `organoid`; this is a label change only and does not
alter the selected cells, genes, transformations, or scores.

## Files

~~~text
config.json   Dataset rules, reference paths, method-result paths, layer keys,
              branch definitions, and fixed display gene lists.
metrics.py    H5AD extraction, trajectory processing, per-gene scores, and
              component/system/method summaries.
plotting/     Figure 5 panel generators, stream plotting, and plot configuration.
README.md     This input, processing, and output contract.
~~~

The paths in config.json are relative to the --artifact-root argument. The
artifact root is an external archive containing the reference tables and the
canonical method-result H5AD files. H5AD files are not tracked by this
repository. The reviewed historical archive uses `references/kinetics/` and
`method-results/kinetics_legacy/`. The latter contains the H5AD inputs selected
for the historical CSV baseline; the separate `method-results/kinetics/`
aliases are current-result inputs and are not used by this reproduction run.
The promoted `results/final_seed123/manifest.json` records the selected
artifact-relative paths and SHA-256 values.

## Method result keys

| Method | System | Extracted values | Input-specific rule |
| --- | --- | --- | --- |
| cellDancer | rpe1, organoid | celldancer_alpha, celldancer_gamma | Direct layers |
| DeepVelo | rpe1, organoid | cell_specific_alpha, cell_specific_gamma | Uses the unspliced-result family and checks velocity_unspliced |
| DeepKINET | rpe1, organoid | degradation_rate as gamma | Direct layer; formal result uses seed 123 |
| Dynamo(label data) | rpe1 | cell_wise_alpha | Read from sparse CSR storage, bin by cycle position, and remove all-zero genes |
| TSvelo | rpe1 | alpha | Read from two lineage H5AD files and average overlapping cells |

## Processing and scoring

For `rpe1`, the evaluator follows the historical sequence:

1. Extract each configured alpha/gamma layer and group cells by
   `Cell_cycle_possition`, taking the mean and reindexing positions 0--300.
2. Intersect the genes across the method parameters and their references, in
   reference-table order.
3. Write the raw grouped matrix to CSV and reload it.
4. For the configured DeepVelo RPE1 result only, clip negative grouped values
   at zero. CellDancer and DeepKINET are not clipped.
5. Apply a centered 20-point moving mean with `min_periods=1`, divide each
   gene by its own smoothed-series median, and calculate
   `log2(ratio + 1e-6)`.
6. Restore the original missing-position mask. The processed trajectory is
   written to CSV and reloaded before scoring; aligned references use the same
   CSV round-trip in memory.

For `organoid`, each branch requires an exact configured Monocle branch label
and one of its configured SOM clusters. Branch 1 uses time `3.5 - pseudotime`
and branch 2 uses `pseudotime - 6.5`; negative aligned times are excluded.
Cells are placed into 200 intervals, averaged within each interval, smoothed
with a centered 5-point moving mean, and divided by the mean of the first 50
bins plus `1e-6`. The processed trajectory is written to CSV and reloaded
before scoring; the aligned reference uses the same CSV round-trip in memory.
No log2 transform is applied to this path.

The evaluator aligns genes and timepoints before calculating per-gene Spearman
correlation, cosine similarity, and the number of valid samples. Component
means are aggregated to system means and then to equal-system method means.
The display-only z-score used by the `rpe1` trajectory plots is not used for
benchmark scores.

The current DeepKINET compact result bundle is generated from the seed-123
H5AD artifacts and contains gamma only; the final H5AD files do not expose a
cell-specific alpha layer.

The kinetics Dynamo input is the labelled-data RPE1 result, displayed as
`Dynamo(label data)`; this is distinct from the global RNA-seq configuration
`Dynamo(RNA_seq)`. The Dynamo all-zero-gene removal and TSvelo overlapping-cell
mean are fixed method-specific rules, not configurable alternatives.

The evaluator applies the same CSV round-trip to aligned reference tables in
memory for scoring. It does not write method-specific `*_reference.csv` files
to the public result bundle; the shared `reference_*.csv` files loaded from
the artifact root remain the complete reference tables used by plotting.

## Run

Validate all configured inputs:

~~~bash
python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --check-only
~~~

Write a fresh metrics output directory:

~~~bash
python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --output /path/to/metrics
~~~

Render Figure 5 quantitative panels and plot-specific tables:

~~~bash
python evaluation/kinetics_recovery/plotting/figure5.py \
  --artifact-root /path/to/kinetics-artifacts \
  --output /path/to/figures
~~~

Render the `organoid` stream panels (the CLI key remains the internal ID):

~~~bash
python evaluation/kinetics_recovery/plotting/streams.py --artifact-root /path/to/velocity --panel organoid
~~~

The metrics output contains component_means.csv, gene_metrics.csv,
system_means.csv, method_summary.csv, extracted matrices under matrices/,
and processed trajectories under trajectories/; method-specific aligned
reference copies are not written. The plotting output contains
descriptively named figures and a plot_data/ directory. Figure 5 quantitative
panels use the legacy input mapping in plotting/figure5_config.json: the
historical method CSVs remain fixed, while the DeepKINET RPE1 and organoid
gamma inputs point to the final seed-123 result. Both kinetic-class panels are
computed from their mapped alpha/gamma trajectories; the CellDancer panel is no
longer copied from a static PNG. The organoid stream panels
used by Figure 5E are generated by plotting/streams.py. The current DeepKINET
RPE1 stream in Figure 5A is generated by plotting/streams.py; the other four
RPE1 stream assets remain explicitly recorded manual assets.
