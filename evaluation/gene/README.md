# Gene-level evaluation

This module owns gene-level velocity evaluation on simulated and real datasets.
It preserves the method-native `Ms`, `Mu`, spliced-velocity, and
unspliced-velocity layers. The reference-neighbor calculation is the only
controlled replacement used by the current CBDir results.

## Reproducible entry point

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
python plot_current_gene_assets_from_originals.py \
  --artifact-root /path/to/velocity-artifacts \
  --output-root evaluation/gene/figures/panels
```

`figures/panels/assembly_manifest.tsv` records the panel-to-file mapping and
SHA256 values. The renderer records external inputs as logical artifact
locators and never writes machine-specific paths into public metadata.

The sampling design is adapted from cellDancer 1.1.7, which is distributed
under the BSD 3-Clause license. The repository implementation does not require
cellDancer to be installed; the retained notice is in
`plotting/THIRD_PARTY_NOTICES.md`.

## Figure display precision

Display precision is separate from the stored metric values. For gene-level
phase portraits, individual gene CBDir-gene annotations are formatted to four
decimal places so that values close to zero retain their sign and magnitude.
Aggregated summaries across genes or datasets (means, standard deviations,
and ranking panels) are formatted to three decimal places. All CSV result
tables retain their full numeric precision and are the source of record for
re-analysis.

The same rule applies to staged and final gene-level figure regeneration:
formatting changes must not alter the underlying CSV values.

The audited simulation inputs have a clean regeneration entry point under
[`simulation/`](simulation/). It uses the three frozen scenario parameter sets
without setting a random seed. Each output is therefore an independent random
realization rather than a replacement for the historical canonical artifact.
