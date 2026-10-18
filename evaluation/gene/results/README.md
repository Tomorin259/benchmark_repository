# Canonical gene-level results

The files in this directory are the fixed-neighbor reference-30 Gene
evaluation with veloVI's additional per-axis min-max transformations reversed
for numerical state and velocity scoring:

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

veloVI scores cover 47 real and 70 simulated genes, with 1,405 edge rows.
`old_cbdir` and `old_edge_score` remain the earlier comparison scores, so their
deltas include neighbor changes and, for veloVI, unit restoration.
[unit_restoration/provenance.json](unit_restoration/provenance.json) records
234 axis factors, paired native/restored scores, input hashes, scoring masks
and counts, and the numerical replay.
[The historical provenance](unit_restoration/historical_reference30_provenance.json)
retains the original generation identity and neighbor metadata.

`simulation/gene_cbdir_cossim_reference.csv` supplies simulation CosSim
and retains the historical CBDir values. Its 630 rows originated from the three
upstream `results_new/score_rerun_results.txt` files (branching: 180,
murk: 90, all: 360), with columns Dataset, Gene, Method_raw, Method, CBDir,
and CosSim. The 70 veloVI CosSim values use restored units; the other 560
values retain their upstream source. Historical source and current table hashes are in
[the compact-table provenance](simulation/gene_cbdir_cossim_reference.provenance.json).
The external historical sources remain in the
[historical panel provenance](../figures/panels/provenance.json). The
[current panel manifest](../figures/panels/assembly_manifest.tsv) binds CBDir
annotations to the canonical reference-30 table; current-generation records
identify the separate CosSim source.

`common_genes/` contains the finite nine-method gene intersection, dataset
means/sample SDs/ranks, and equal-dataset macro means/ranks for the seven real
datasets. Rebuild or check these three tables with
`evaluation/gene/refresh_common_gene_summary.py`.

Current Figure 4D uses `new_cbdir` from `gene_cbdir_reference30_all.csv`
for CBDir-gene and CosSim from the compact simulation table, paired with
cell-level CBDir and BATC. The final panel and its
`results/cross_level_summary.csv` are owned by
the [Cell simulation module](../../cell/simulation/README.md).
