# Gene-matched cell-level CBDir/BATC rerun

This supplementary rerun uses the repository's existing `evaluation/cell/run_cell_metrics.py`
entry point with `evaluation/cell/gene_matched/inputs.tsv`. It contains the gene-analysis H5AD
outputs for DeepVelo, UniTVelo, and veloVI. This directory is separate from the canonical Cell
rerun so those gene-matched outputs are never mixed with canonical DeepVelo or veloVI files.
The six gene-kinetic datasets are included; `capo_pre` is evaluated as `capofibro` according to
the Cell configuration.

`results/cell_metrics_edge_details_extended.csv` contains native two-dimensional CBDir and BATC
diagnostics. `batc_candidate_source_mean` is the pre-branch-max source-cell candidate mean after
the existing NaN-to-zero policy; it is not the overall BATC. The winner-count field records the
source cells selecting each outgoing target after the existing branch maximization. Overall BATC
is the scalar in `cell_metrics.csv`.

Validation completed for 18/18 pairs with status `ok`. UniTVelo scalar values were checked against
the canonical Cell rows and matched within 1e-12. DeepVelo and veloVI are expected to differ where
their gene-matched H5AD outputs differ from canonical Cell outputs.
