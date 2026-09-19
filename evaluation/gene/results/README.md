# Canonical gene-level results

The files in this directory are the fixed-neighbor reference-30 Gene
evaluation. The two earlier runs had complementary dataset coverage and are
represented here as one date-free result set:

```text
real/<dataset>/
simulation/<scenario>/
reference_neighbors/<dataset>.npz
gene_cbdir_reference30_all.csv
gene_cbdir_edge_details_reference30_all.csv
neighbor_audit_all.csv
provenance.json
```

The neighbor NPZ files are committed so that the exact reference neighbor
sets and cell ordering can be inspected without access to the original H5AD
files. Method-result H5AD files and simulation inputs remain external
artifacts. Their public locators use `artifact://velocity/<relative-path>`;
the evaluator accepts the artifact root through `--artifact-root` or
`VELOCITY_ARTIFACT_ROOT`.

The metric values are copied from the fixed-neighbor results without
recalculation. The additional edge, neighbor, and comparison columns are
retained for supplementary inspection; manuscript summaries use the CBDir
columns.

The simulation Figure 3D source table is
simulation/gene_cbdir_cossim_reference.csv. It contains the 630 rows from
the three upstream results_new/score_rerun_results.txt files
(branching: 180, murk: 90, all: 360), with columns Dataset, Gene,
Method_raw, Method, CBDir, and CosSim. The source files remain external
artifacts and are recorded in the gene figure assembly manifest.
