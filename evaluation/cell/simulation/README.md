# Simulated-dataset cell-level evaluation

This directory owns a supplementary, method-output-matched cell-level
evaluation of the three historical simulated datasets. It is separate from
the canonical 14-real-dataset, 20-configuration Cell benchmark. A new
simulation realization must not silently replace the historical inputs.

## Inputs and comparison unit

- Scenarios: `murk` (normal + MURK), `branching` (normal + branching), and
  `all` (hybrid). Their historical reference H5ADs are identified by
  `SIM_CONFIGS` in `evaluation/gene/run_reference30_gene_cbdir.py`.
- Score exactly the nine method-output H5ADs used by the gene-level
  simulation evaluation. Resolve the source H5AD and spliced velocity layer
  from `evaluation/gene/results/gene_cbdir_reference30_all.csv`; use the
  method registry for operational and paper-facing names. The comparison
  unit is one scenario by one method (27 pairs), not one gene.
- Use `obs['bin']` and the directed adjacent-bin edges defined in
  `SIM_CONFIGS`: nine Branch-B edges for `murk`, and nine edges on each of
  Branch A and Branch B for `branching` and `all`.

## Shared reference state

The simulation Gene neighbor assets contain cell names and neighbor indices
only. The Cell evaluator additionally requires distances, connectivities,
and common `Ms`/`Mu` matrices. Build these from each historical reference
H5AD using the fixed 30-neighbor protocol, then verify that cell names and
neighbor indices agree exactly with the Gene assets. All methods within a
scenario must use the same reference graph, expression state, and UMAP. Do
not substitute a method-specific neighbor graph or regenerate the simulation.
`config.json` points `reference_umap_path` to the historical simulation H5AD
for each scenario. The Cell preparation code uses that common UMAP instead of
requiring an additional scVelo-deterministic output. DeepVelo's simulation
output has no stored `velocity_umap`; its `reference_graph_projection`
alignment computes the embedding velocity on the same reference graph. The
other eight methods retain their native embedding-velocity path. The compact
graph assets and their source hashes are recorded in `reference_neighbors/`
and `reference_manifest.json`. Full `Ms`/`Mu` reference moments are external
artifacts and are not committed to Git.

## Results and presentation

Use the existing Cell definitions for Consistency, ICCoh, CBDir,
TransCosine, and BATC; retain CBDir-HD and edge-level/evaluability audits as
supplementary diagnostics. Keep full precision in numerical result tables.
The principal cross-level comparison pairs cell-level CBDir and BATC with the
mean CBDir-gene and ground-truth CosSim of the corresponding simulated gene
set. CosSim comes from
`evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv` and is
joined one-to-one to the reference-30 CBDir-gene results by scenario, gene,
and method. The gene sets must match exactly; neither metric is recalculated.
In the hybrid scenario,
the MURK and branching gene subsets remain separate; they share the same
scenario-level cell score. Do not merge these rows into the real-data Cell
table or the 20-method overall ranking.

The completed bundle has 27 method-by-scenario rows in each primary,
supplementary, and source-hash audit table, 810 native edge-detail rows in
each edge table, and 54 TiVelo parity rows. All 27 primary metric quintets
are finite; all parity checks pass. `cross_level_summary.csv` has 36 rows:
nine methods for each of the normal + MURK, normal + branching, hybrid MURK,
and hybrid branching gene subsets. The hybrid subsets use the same cell-level
score but distinct gene-level means and sample standard deviations for both
CBDir-gene and CosSim. These
standard deviations describe variation across the evaluated genes, not
independent simulation runs; there is one historical realization per scenario.

`figures/simulation_cell_cbdir_batc` compares CBDir and BATC across the three
scenarios. `figures/simulation_cross_level_cbdir` compares CBDir and BATC with
the mean CBDir-gene and CosSim, each with gene-wise SD, for the four subsets.
Each method occupies one alternating light-gray/white band spanning all four
panels and their gutters, with cell-level scores above and gene-level scores
below. The panel canvas is 9.0 by 4.125 inches.
Small offsets within each subrow avoid marker occlusion. Blue circles encode
CBDir-gene, purple triangles encode CosSim, orange diamonds encode CBDir, and
green squares encode BATC; marker shapes, colors, and sizes are unchanged.
Each stem is exported as PDF, SVG, and PNG. The Figure 4D cross-level PNG uses
600 dpi; its graphic contains the legend and plot labels, with explanatory
prose reserved for the manuscript caption. The panels use the relative method
order and names in `methods/configurations.tsv`; the CSV files retain full
precision, with no ranking based on rounded display values. BATC is shown as
a separate cell-level overview and is not treated as a gene-level score.

The final Figure 4D source is registered in the root paper-assembly manifest.
`figures/figFour.pdf` is the approved, user-exported and cropped composite;
`figures/fig4_composite_manifest.json` records its checksum and D-panel source.
The renderer reproduces the D panel and its summary, not the manually
assembled whole-figure PDF.

## Reproduce

### Supplementary Figure 12 and Table 7

`results/cell_metrics.csv` is the full-precision source for the 27 rows and
135 scores in Supplementary Table 7. The table CSV
`results/simulation_cell_scores_4dp.csv` and the manuscript-ready fragment
`results/simulation_cell_scores.tex` report raw scores, not ranks, to four
decimal places. `results/table_manifest.json` records the source checksum,
method order, formatting, and LaTeX dependencies. No second copy of the
full-precision results is maintained.

The table exporter reuses the final material bundle's original table script;
its input paths and method ordering now resolve through this module's inputs
and the repository method registry. Regenerate or check the frozen exports:

```bash
python evaluation/cell/simulation/build_score_table.py
python evaluation/cell/simulation/build_score_table.py --check
```

These commands do not rerun velocity methods, change metric definitions, or
regenerate the simulation. The check compares all 135 formatted scores and
the complete LaTeX fragment against the full-precision source.

The approved, user-exported and cropped Supplementary Figure 12 composite is
stored once at `figures/suppFigTwelve.pdf` at the repository root.
`figures/suppfig12_composite_manifest.json` records its checksum and source
boundary. It contains the three nine-method stream rows and the five-metric
boxplot row; the heatmap was replaced by Supplementary Table 7. The table
exporter does not regenerate this manually assembled PDF or its stream panels.

### Evaluation and comparison panels

From the repository root, using the historical artifacts and compatible
scVelo/TiVelo environments:

```bash
python evaluation/cell/simulation/build_simulation_inputs.py \
  --artifact-root /path/to/artifacts \
  --reference-moments-root /path/to/external-reference-moments

python evaluation/cell/run_cell_metrics.py \
  --artifact-root /path/to/artifacts \
  --inputs evaluation/cell/simulation/inputs.tsv \
  --config evaluation/cell/simulation/config.json \
  --reference-neighbors-root evaluation/cell/simulation/reference_neighbors \
  --reference-moments-root /path/to/external-reference-moments \
  --output-root evaluation/cell/simulation \
  --work-root /path/to/external-temporary-work \
  --indicators-python /path/to/indicators-python \
  --tivelo-python /path/to/tivelo-python \
  --prepare-script evaluation/cell/prepare_cell_metrics.py \
  --official-script evaluation/cell/score_tivelo_metrics.py

python evaluation/cell/simulation/plot_simulation_cell_comparison.py
```

Pass `--cross-level-only` to refresh the paired summary table and Figure 4D
replacement without regenerating the three-scenario Cell overview.

The reference build uses the Gene simulation preprocessing and checks its
cell names and 30-neighbor indices against the stored Gene reference assets.
`config.json` fixes 20 BATC bins for all three scenarios and scores the same
adjacent-bin directed transitions as the Gene simulation benchmark: nine on
Branch B in `murk`, and nine on each branch in `branching` and `all`. Neither
the original H5ADs nor the canonical 14×20 real-data Cell tables are changed.
