# Cell-level CBDir/BATC rerun

This supplementary rerun uses the repository's existing `evaluation/cell/run_cell_metrics.py`
entry point with the fixed reference-30 neighbors, the configured reference moments, and the
canonical `inputs.tsv`/`config.json`. It covers six datasets and the nine methods used by the
gene-kinetic phase portraits. The Cell configuration maps `capo_pre` to the evaluation dataset
`capofibro`.

The original `evaluation/cell/results/cell_metrics_edge_details.csv` semantics are retained in
`results/cell_metrics_edge_details.csv` for high-dimensional CBDir and TransCosine. The new
`results/cell_metrics_edge_details_extended.csv` contains native two-dimensional CBDir and BATC
diagnostics. Its `batc_candidate_source_mean` rows are the mean candidate cosine over source cells
for one directed edge after the existing NaN-to-zero policy, before branch maximization. They are
not the overall BATC score. `branch_winner_source_cells` records how many source cells selected
that target after the existing per-source branch maximization. The overall BATC remains the scalar
in `cell_metrics.csv` and is computed by the unchanged aggregate path.

Validation completed for 54/54 dataset-method pairs: every pair has status `ok`, all extended
scores and source-cell counts are finite/positive, branch winner counts sum to the source-cell
count for every source cluster, and the CBDir/BATC scalar absolute differences from the canonical
Cell table are both 0 (tolerance 1e-12).
