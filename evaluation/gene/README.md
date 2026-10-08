# Gene-level evaluation

This module owns gene-level velocity evaluation on simulated and real datasets.
It preserves the method-native `Ms`, `Mu`, spliced-velocity, and
unspliced-velocity layers. The reference-neighbor calculation is the only
controlled replacement used by the current CBDir results.

## Reproducible entry point

`config.py` owns the shared input paths, velocity-layer keys, selected gene
lists, and transition edges used by evaluation and panel rendering. Display-only
summary dimensions and typography live in `plotting/layouts.json`.

`run_reference30_gene_cbdir.py` evaluates the nine gene-level methods on:

- real: `retina`, `DentateGyrus_lamanno`, `gastrulation`, `capo_pre`,
  `bonemarrow`, `mousebrain`, and `pancreas`;
- simulation: `branching`, `murk`, and `all`.

Large method-result H5AD files and simulation inputs are external artifacts.
Pass their common root explicitly, or set `VELOCITY_ARTIFACT_ROOT`:

```bash
python run_reference30_gene_cbdir.py \
  --artifact-root /path/to/velocity-artifacts \
  --neighbor-mode reference \
  --output-root evaluation/gene/results
```

The canonical neighbor mode uses the shared 30-neighbor assets. The native
neighbor mode is a controlled comparison only:

```bash
python run_reference30_gene_cbdir.py \
  --artifact-root /path/to/velocity-artifacts \
  --neighbor-mode native \
  --output-root evaluation/gene/results/native_control
```

The committed result layout is organized by scientific scope:

```text
results/
├── real/<dataset>/
├── simulation/<scenario>/
├── reference_neighbors/<dataset>.npz
└── gene_cbdir_reference30_all.csv
```

The two earlier reference-30 runs had complementary coverage; their values
are now represented in this single date-free layout. The reference-neighbor
NPZ files are retained so that the fixed neighbor sets can be inspected
without reopening the original H5AD files. The result tables retain the
method-specific audit columns, while the compact gene-level CBDir columns are
the fields used by manuscript summaries.

## Real-data gene selection provenance

[`selection/`](selection/) preserves the recovered MURK and branching candidate
rankings, upstream/project screening sources, and the frozen real-data evaluation
lists. Its reconstruction checks separate historical evidence from newly run
results and distinguish computational screening from subsequent manual selection.
Dataset-specific configurations that have not been recovered remain documented.
The checks do not modify the reference-30 evaluation results.

## Repository-owned phase portraits

`plotting/phase_portrait.py` is the public phase-portrait renderer. It changes
only vector display and never writes to the source AnnData object or to any
metric input. Phase coordinates stay linear, representative cells are selected
deterministically, and no neighborhood velocity smoothing is applied.

The single canonical arrow contract is:

```text
length mode:          log1p
reference quantile:   0.90
reference length:     0.05
maximum length:       0.10
grid:                 25 x 25
density percentile:  15
sampling:             deterministic
```

This is the enlarged-arrow contract used by the paper panels. It is the only
display contract supported by the public renderer.

The paper panels are materialized under one date-free directory:

```bash
python evaluation/gene/plot_gene_assets.py \
  --artifact-root /path/to/velocity-artifacts \
  --output-root outputs/gene-panels
```

The renderer is owned entirely by this repository. It reads reference-30
`new_cbdir` from `results/gene_cbdir_reference30_all.csv` and joins `CosSim`
from `results/simulation/gene_cbdir_cossim_reference.csv` by dataset, method,
and gene. The compact simulation table's separate CBDir column is not a
replacement for the reference-30 values used by these panels. Missing tables,
duplicate keys, mismatched inputs, and missing velocity components are errors;
there is no external score-TXT fallback or dynamic import of legacy scripts.

For CSV-only summary boxplots/heatmaps, no H5AD root, Scanpy, or scVelo import
is needed:

```bash
python evaluation/gene/plot_gene_assets.py --summaries-only \
  --output-root outputs/gene-summary
```

For a focused phase-panel check, use `--cases bonemarrow --genes LMNA`;
summaries still cover every registered gene in the selected dataset. `--dpi`
defaults to 600. Output directories must be empty and must not be the committed
`figures/panels/` directory; the default is `outputs/gene-panels/`. Independently
review a new export before promoting any assets to the paper bundle.

The phase API is unchanged. Native `Ms`/`Mu` are read without recomputing
moments or neighbor graphs. All 90 frozen result inputs have both coordinate
layers; missing layers stop the export. The historical real-panel velocity
cleaning (NaN/Inf to zero and clipping to +/-1e10) applies only to selected
plotting copies. Simulation vectors retain their native values. Real cluster
labels/colors are taken from each registered reference H5AD; the historical
equal-length label-assignment fallback is retained. H5AD files are never saved
or modified by the renderer. Summary stripplot jitter inherits the display RNG
and has no effect on metric values, means, or ordering.

Each export writes `scores_used.csv`, source/configuration hashes in
`provenance.json`, and an `assembly_manifest.tsv` for generated paper-panel
references. It generates the previously scripted nine-method comparison panels
and summaries, not the manually assembled whole PDFs or retained raw phase
plates. The committed `figures/panels/assembly_manifest.tsv` remains the
panel-to-file mapping and hash authority for current paper assets.

After a display-only update, check or refresh the existing panel hashes from
the repository root without importing scVelo, opening H5AD files, or redrawing:

```bash
python evaluation/gene/refresh_panel_manifest.py
python evaluation/gene/refresh_panel_manifest.py --write
```

The refresh preserves panel references, row order, and recorded score sources;
it changes only SHA256 fields. Shared images may have separate references for
main and supplementary panels. `figures/panels/provenance.json` retains the
historical generation record, not the current script or image checksums.
Current image hashes belong to the panel manifest. Its Figure 4D scenario
summaries are supporting gene-only analyses; the final paired Cell/Gene panel
is owned by [`../cell/simulation/`](../cell/simulation/).

The sampling design is adapted from cellDancer 1.1.7, which is distributed
under the BSD 3-Clause license. The repository implementation does not require
cellDancer to be installed; the retained notice is in
`plotting/THIRD_PARTY_NOTICES.md`.

## Figure display precision

Display precision is separate from the stored metric values. Individual
method-by-dataset, method-by-gene, or method-by-transition scores are formatted
to four decimal places, including both CBDir-gene and cell-level CBDir in a
single-gene phase portrait. Positive scores do not carry a leading plus sign;
negative signs are retained. Aggregated summaries across genes or datasets
(means, standard deviations, medians, macro-averages, branch averages,
differences between summary values, and ranking-panel score annotations) use
four decimal places. Dense heatmap annotations use two decimal places.
Formal manuscript numerical tables use four decimal places, including their
summary entries. Sorting and ranking always use unrounded values. All CSV
result tables retain their full numeric precision and are the source of
record for re-analysis.

The same rule applies to staged and final gene-level figure regeneration:
formatting changes must not alter the underlying CSV values.

The audited simulation inputs have a clean regeneration entry point under
[`simulation/`](simulation/). It uses the three frozen scenario parameter sets
without setting a random seed. Each output is therefore an independent random
realization rather than a replacement for the historical canonical artifact.
