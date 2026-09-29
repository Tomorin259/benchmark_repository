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
The principal cross-level comparison pairs cell-level CBDir with the mean
CBDir-gene of the corresponding simulated gene set. In the hybrid scenario,
the MURK and branching gene subsets remain separate; they share the same
scenario-level cell score. Do not merge these rows into the real-data Cell
table or the 20-method overall ranking.

The completed bundle has 27 method-by-scenario rows in each primary,
supplementary, and source-hash audit table, 810 native edge-detail rows in
each edge table, and 54 TiVelo parity rows. All 27 primary metric quintets
are finite; all parity checks pass. `cross_level_summary.csv` has 36 rows:
nine methods for each of the normal + MURK, normal + branching, hybrid MURK,
and hybrid branching gene subsets. The hybrid subsets use the same cell-level
score but distinct gene-level means and gene-wise standard deviations. These
standard deviations describe variation across the evaluated genes, not
independent simulation runs; there is one historical realization per scenario.

`figures/simulation_cell_cbdir_batc` compares the two cell directionality
metrics across the three scenarios. `figures/simulation_cross_level_cbdir`
compares Cell CBDir with Gene CBDir mean and gene-wise SD for the four subsets.
Each stem is exported as PDF, SVG, and PNG. The panels use the relative method
order and names in `methods/configurations.tsv`; the CSV files retain full
precision, with no ranking based on rounded display values. BATC is shown as
a separate cell-level overview and is not treated as a gene-level score.

## Reproduce

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

The reference build uses the Gene simulation preprocessing and checks its
cell names and 30-neighbor indices against the stored Gene reference assets.
`config.json` fixes 20 BATC bins for all three scenarios and scores the same
adjacent-bin directed transitions as the Gene simulation benchmark: nine on
Branch B in `murk`, and nine on each branch in `branching` and `all`. Neither
the original H5ADs nor the canonical 14×20 real-data Cell tables are changed.
