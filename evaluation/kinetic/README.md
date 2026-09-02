# Case 3 kinetic-parameter evaluation

This directory is the complete Case 3 analysis layer. It consumes existing
method result objects; it does not prepare datasets or train methods.

## Files

~~~text
data/case3.json  Direct paths to references, method H5ADs, parameter layers,
                  and the fixed Figure 5 gene lists.
metrics.py        H5AD extraction, RPE1/organoid processing, per-gene
                  Spearman/cosine calculation, and Figure 6 summaries.
plots.py          Plot-ready tables plus compact Figure 5/6 quantitative
                  review panels from an existing metrics output directory.
README.md         This interface description.
~~~

The data map is deliberately small. cellDancer, DeepVelo, and DeepKINET are
read directly from their existing H5AD outputs. DeepVelo is mapped to its
unspliced-result family.

Dynamo and TSvelo have no generic, one-file parameter-layer contract, but
they are still read from H5AD rather than from historical metrics:

- Dynamo reads sparse cell_wise_alpha from the labelled RPE1 H5AD, bins
  Cell_cycle_possition with the original pd.cut(range(0, 302)) rule, and
  removes all-zero genes.
- TSvelo reads alpha from l0_TSvelo.h5ad and l1_TSvelo.h5ad, averages
  duplicate cells across lineages, keeps the first cell-cycle label, then
  bins positions 0--300.

sceu_alpha.csv remains only as the historical 250-gene selection. The
reference values are always reread from Table S1. The old Dynamo/TSvelo alpha
CSV and Spearman files are regression baselines, never runtime inputs.

After a full metrics run, verify the two special H5AD reconstructions with:

~~~bash
python evaluation/kinetic/validate_historical_extraction.py +  --working-root /mnt/second19T/zhaozelin/velocity +  --metrics /path/to/new_case3_metrics
~~~

This test checks the archived raw matrix, processed matrix, and per-gene score.
The TSvelo score tolerance is intentionally 1e-5 to accommodate its historical
library-version-level floating-point difference.

## Run metrics

~~~bash
python evaluation/kinetic/metrics.py \
  --working-root /mnt/second19T/zhaozelin/velocity \
  --output /path/to/new_case3_metrics
~~~

Use --methods to evaluate a subset, or --check-only to verify that all direct
inputs, labels, and parameter layers can be read without writing any files.

The output directory must be new or empty. It contains:

~~~text
component_means.csv  One row per method, system component, and parameter.
gene_metrics.csv     Per-gene Spearman/cosine results for plotting.
system_means.csv     RPE1 and organoid system summaries.
method_summary.csv   Figure 6 style equal-system summaries.
matrices/            Reference and extracted, binned parameter matrices.
trajectories/        Processed parameter trajectories for Figure 5 panels.
manifest.json        Sources, adapters, and every generated matrix path.
~~~

## Figure boundary

~~~bash
python evaluation/kinetic/plots.py \
  --metrics /path/to/new_case3_metrics \
  --output /path/to/new_case3_figures
~~~

plots.py writes the plot tables and compact quantitative review panels. It does
not recompute metrics or reopen H5AD files. Cell-embedding velocity-stream
panels remain in the cell-level visualization workflow and are deliberately
not duplicated here; this directory owns the parameter-recovery trajectories,
boxplots, and summaries.
