# Gene-output-matched cell-level results

Cell metrics for six gene-kinetic datasets and three methods: DeepVelo, UniTVelo
and veloVI (18 dataset–method pairs). Evaluation uses their Gene-analysis H5ADs
with `evaluation/cell/run_cell_metrics.py` and
`evaluation/cell/gene_matched/inputs.tsv`. These source files are recorded
separately from the canonical Cell inputs. The Cell configuration evaluates
`capo_pre` as `capofibro`.

`results/cell_metrics_edge_details_extended.csv` contains native two-dimensional CBDir and BATC
diagnostics. `batc_candidate_source_mean` is the pre-branch-max source-cell candidate mean after
the existing NaN-to-zero policy; it is not the overall BATC. The winner-count field records the
source cells selecting each outgoing target after the existing branch maximization. Overall BATC
is the scalar in `cell_metrics.csv`.

All 18 pairs have status `ok`. UniTVelo scalar scores match the canonical Cell
rows within 1e-12. DeepVelo and veloVI use distinct Gene-analysis outputs, so
their scores can differ from the canonical Cell results.
