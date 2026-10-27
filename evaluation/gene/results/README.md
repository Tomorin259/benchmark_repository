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

Use `simulation/gene_cbdir_cossim_reference.csv` for simulation CosSim. Its
`CBDir` column retains historical comparison values; current CBDir is
`new_cbdir` in `gene_cbdir_reference30_all.csv`.
Its 630 rows cover branching (180), MURK (90), and hybrid (360), with columns
Dataset, Gene, Method_raw, Method, CBDir, and CosSim. CosSim uses the paired
noise-free ODE velocities in [`../simulation/truth/`](../simulation/truth/);
all values retain full precision. The 70 veloVI values use restored units.
The [Gene evaluator](../run_reference30_gene_cbdir.py) recalculates all nine
methods into `simulation/gene_cossim.csv` in a fresh output root. Compare by
dataset/method/gene with absolute tolerance `1e-12`.
Table, velocity-reference and coordinate-unit sources are in
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
