# Simulated-dataset Cell evaluation

Compare the nine saved Gene-evaluation method outputs across three historical
simulation realizations. This 27-pair comparison is separate from the
14-real-dataset, 20-configuration Cell benchmark.

## Inspect results and check the table

Run commands from the repository root.

| Material | Contents |
| --- | --- |
| [results/cell_metrics.csv](results/cell_metrics.csv) | 27 scenario–method rows; five Cell metrics and 135 full-precision scores |
| [results/](results/) | Supplementary metrics and source-hash audits (27 rows each), native edge tables (810 rows each), and 54 TiVelo parity rows |
| [results/cross_level_summary.csv](results/cross_level_summary.csv) | 36 paired rows: nine methods across four Gene subsets |
| [figures/](figures/) | Cell CBDir/BATC overview and the Figure 4D cross-level source panel |
| [results/simulation_cell_scores.tex](results/simulation_cell_scores.tex) | Supplementary Table S7; four-decimal raw scores and [table manifest](results/table_manifest.json) |

Check the saved Table S7 exports without writing files:

```bash
python evaluation/cell/simulation/build_score_table.py --check
```

The check compares all 135 formatted scores and the complete LaTeX fragment
against the full-precision source. The saved primary metric quintets are finite,
and the retained TiVelo parity checks pass.

## Comparison and reference state

- Scenarios are `murk` (normal + MURK), `branching` (normal + branching), and
  `all` (hybrid). Historical reference H5ADs and directed adjacent-bin edges
  follow `SIM_CONFIGS` in [the Gene evaluator](../../gene/run_reference30_gene_cbdir.py):
  nine Branch-B edges for `murk`, nine per branch for the other scenarios.
- Method H5ADs and spliced velocity layers match the
  [Gene score sources](../../gene/results/gene_cbdir_reference30_all.csv).
  The comparison unit is one scenario by one method, using the registered
  historical inputs.
- Every method within a scenario shares the fixed 30-neighbor reference graph,
  expression state and historical UMAP. The Cell reference adds distances,
  connectivities and `Ms`/`Mu` to the Gene reference's cell names and indices;
  the latter must agree exactly. [reference_manifest.json](reference_manifest.json)
  records source hashes; full reference moments remain external artifacts.
- [config.json](config.json) identifies the historical reference UMAP and fixes
  20 BATC bins. DeepVelo's missing `velocity_umap` uses
  `reference_graph_projection` on that shared graph; the other eight methods
  retain their native embedding-velocity path.

Consistency, ICCoh, CBDir, TransCosine and BATC use the
[Cell definitions](../README.md); CBDir-HD and edge/evaluability tables are
additional diagnostics. Cross-level rows pair Cell CBDir/BATC with finite-gene
means and sample SDs of reference-30 CBDir-gene and
[ground-truth CosSim](../../gene/results/simulation/gene_cbdir_cossim_reference.csv).
The join is one-to-one by scenario, method and gene, with identical gene sets;
the existing Gene metrics are not recalculated.

Hybrid MURK and branching subsets share a scenario-level Cell score but retain
separate Gene means and SDs. SD describes variation across genes, with one
historical realization per scenario. The simulation comparison has its own
result table and scenario-specific summaries.

## Re-evaluate saved outputs

Materialize the historical artifacts and use compatible scVelo/TiVelo
environments. Use new directories for reference construction, moments,
evaluation and temporary work:

```bash
python evaluation/cell/simulation/build_simulation_inputs.py \
  --artifact-root /path/to/artifacts \
  --reference-moments-root /path/to/new-reference-moments \
  --output-root outputs/simulation-inputs

python evaluation/cell/run_cell_metrics.py \
  --artifact-root /path/to/artifacts \
  --inputs outputs/simulation-inputs/inputs.tsv \
  --config evaluation/cell/simulation/config.json \
  --reference-neighbors-root outputs/simulation-inputs/reference_neighbors \
  --reference-moments-root /path/to/new-reference-moments \
  --output-root outputs/simulation-evaluation \
  --work-root /path/to/new-temporary-work \
  --indicators-python /path/to/indicators-python \
  --tivelo-python /path/to/tivelo-python \
  --prepare-script evaluation/cell/prepare_cell_metrics.py \
  --official-script evaluation/cell/score_tivelo_metrics.py
```

Reference construction follows Gene simulation preprocessing and verifies its
cell names and neighbor indices against the stored Gene assets. Evaluation
scores the saved method outputs and writes under
`outputs/simulation-evaluation/results/`.

## Export tables and render retained panels

The table exporter reads this module's saved results and input map, using
registry method order. Write a new copy with:

```bash
python evaluation/cell/simulation/build_score_table.py \
  --output-root outputs/simulation-tables
```

The plotter reads the saved module tables, not the fresh evaluation directory
above. It requires `--output-dir` and writes both plots and a new
`cross_level_summary.csv` there. Use an empty directory outside the repository
or under `outputs/`; existing result tables and paper assets stay unchanged:

```bash
python evaluation/cell/simulation/plot_simulation_cell_comparison.py \
  --output-dir outputs/simulation-panels
```

`--cross-level-only` skips the Cell overview but still writes the summary and
Figure 4D source into that output directory. Both stems export PDF/SVG/PNG;
the Figure 4D PNG uses 600 dpi.
Its 9.0 × 4.125-inch canvas uses alternating method bands, offset markers and
gene-wise SD. Blue circles denote CBDir-gene, purple triangles CosSim, orange
diamonds Cell CBDir, and green squares BATC. Names and relative method order
come from [the registry](../../../methods/configurations.tsv); scores retain
full precision rather than rounded-value ranks.

The author-exported composites are [Figure 4](../../../figures/figFour.pdf)
and [Supplementary Figure 12](../../../figures/suppFigTwelve.pdf), with source
boundaries in their [Figure 4 manifest](../../../figures/fig4_composite_manifest.json)
and [S12 manifest](../../../figures/suppfig12_composite_manifest.json).
S12 contains three nine-method stream rows and a five-metric boxplot row;
its former heatmap is represented by Table S7. The renderer supplies comparison
panels; the manifests identify the remaining stream and assembly sources.
The Figure 4 manifest retains the renderer hash for its original generation;
that source is available at [revision 1e2d236](https://github.com/Tomorin259/benchmark_repository/blob/1e2d236fc163ab15cf4181e27383b7175676f9b1/evaluation/cell/simulation/plot_simulation_cell_comparison.py).
