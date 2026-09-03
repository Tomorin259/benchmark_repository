# Kinetic-parameter recovery

This module evaluates method-estimated kinetic parameters against the scEU-seq
reference tables for RPE1 and organoid systems. It consumes saved method-result
H5AD files; it does not prepare datasets or train methods.

## Files

~~~text
config.json   Dataset rules, reference paths, method-result paths, layer keys,
              branch definitions, and fixed display gene lists.
metrics.py    H5AD extraction, trajectory processing, per-gene scores, and
              component/system/method summaries.
plots.py      Figures and plot-specific tables from a metrics output directory.
README.md     This input, processing, and output contract.
~~~

The paths in config.json are relative to the --artifact-root argument. The
artifact root is an external archive containing the reference tables and the
canonical method-result H5AD files. H5AD files are not tracked by this
repository.

## Method result keys

| Method | System | Extracted values | Input-specific rule |
| --- | --- | --- | --- |
| cellDancer | RPE1, organoid | celldancer_alpha, celldancer_gamma | Direct layers |
| DeepVelo | RPE1, organoid | cell_specific_alpha, cell_specific_gamma | Uses the unspliced-result family and checks velocity_unspliced |
| DeepKINET | RPE1, organoid | degradation_rate as gamma | Direct layer |
| Dynamo | RPE1 | cell_wise_alpha | Read from sparse CSR storage, bin by cycle position, and remove all-zero genes |
| TSvelo | RPE1 | alpha | Read from two lineage H5AD files and average overlapping cells |

## Processing and scoring

For RPE1, cells are grouped by Cell_cycle_possition into positions 0--300.
The configured methods are then smoothed with a centered window of 20,
median-normalized per gene, and transformed with log2(value + epsilon).
The configured DeepVelo RPE1 result is clipped at zero before this transform.

For organoid, each branch is selected using the configured branch and cluster
labels. Branch 1 uses reverse time from 3.5 and branch 2 uses forward time from
6.5. Cells are placed into 200 pseudotime bins, smoothed with a centered window
of 5, and normalized by the mean of the first 50 bins.

The evaluator aligns genes and timepoints before calculating per-gene Spearman
correlation, cosine similarity, and the number of valid samples. Component
means are aggregated to system means and then to equal-system method means.
The display-only z-score used by the RPE1 trajectory plots is not used for
benchmark scores.

The Dynamo all-zero-gene removal and TSvelo overlapping-cell mean are fixed
method-specific rules, not configurable alternatives.

## Run

Validate all configured inputs:

~~~bash
python evaluation/parameter_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --check-only
~~~

Write a fresh metrics output directory:

~~~bash
python evaluation/parameter_recovery/metrics.py \
  --artifact-root /path/to/artifacts \
  --output /path/to/metrics
~~~

Render figures and plot-specific tables:

~~~bash
python evaluation/parameter_recovery/plots.py \
  --metrics /path/to/metrics \
  --output /path/to/figures
~~~

The metrics output contains component_means.csv, gene_metrics.csv,
system_means.csv, method_summary.csv, extracted matrices under matrices/,
and processed trajectories under trajectories/. The plotting output contains
descriptively named figures and a plot_data/ directory. Cell-embedding
velocity-stream panels belong to the cell-level visualization workflow and are
not duplicated here.
