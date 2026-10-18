# Cell-level results for gene–cell comparisons

Cell metrics for six datasets and the nine methods used by the gene-kinetic
phase portraits (54 dataset–method pairs). Evaluation uses
`evaluation/cell/run_cell_metrics.py`, fixed reference-30 neighbors, configured
reference moments and the canonical `inputs.tsv`/`config.json`. The Cell
configuration maps `capo_pre` to `capofibro`.

`results/cell_metrics_edge_details.csv` uses the same high-dimensional CBDir and TransCosine
definitions as `evaluation/cell/results/cell_metrics_edge_details.csv`.
`results/cell_metrics_edge_details_extended.csv` contains native two-dimensional CBDir and BATC
diagnostics. Its `batc_candidate_source_mean` rows are the mean candidate cosine over source cells
for one directed edge after the existing NaN-to-zero policy, before branch maximization. They are
not the overall BATC score. `branch_winner_source_cells` records how many source cells selected
that target after the existing per-source branch maximization. The overall BATC remains the scalar
in `cell_metrics.csv` and is computed by the unchanged aggregate path.

All 54 pairs have status `ok`. Extended scores are finite and source-cell counts
are positive; branch winner counts sum to the source-cell count for each source
cluster. CBDir/BATC scalars match the canonical Cell table exactly (comparison
tolerance 1e-12).
