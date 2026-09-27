# Supplementary Figure 11: veloAE representation diagnostics

## Scope

The online `merged.tex` defines one two-panel supplementary figure:

- **A — PCA summary.** PC1 explained variance for the state matrix (`Ms`, left)
  and velocity (right), comparing veloAE, scVelo stochastic, and TIvelo across
  14 real datasets. Filled points are raw PCA; open points are PCA after
  per-dimension standardization. Lines pair the same dataset and black bars
  show cross-dataset means.
- **B — shared-axis projection plate.** For each of the 14 datasets, state
  coordinates and velocity vectors are shown in the same PC1–PC2 axes fitted to
  raw `Ms`. State coordinates are centered; velocity is projected without
  centering. Points use the dataset annotations, and arrows are spatial-bin
  averages with display-only scaling and length caps.

There is no manuscript panel C and no independent S10 numbering in this
source bundle. A final page-level composite is a downstream manual-layout
artifact; the repository's authoritative source panels are kept separately.

## Authoritative source and output paths

- `plot_pca_summary.py` reads the existing PCA comparison table and writes the
  A source asset.
- `plot_shared_projection.py` reads the 14 trusted projection NPZ archives,
  resolves display names through `data/datasets.tsv`, and writes the individual
  panels plus the B plate.
- `export_ppt_assets.py` exports the two native A/B handoff assets under
  `results/ppt_assets_final/`.
- `results/ppt_assets_final/manifest.json` records output hashes and the 5×3
  layout of panel B.
- `figures/assembly/assembly_manifest.tsv` maps the formal repository sources
  to `figures/assembly/suppfig11/A_pca_summary.pdf` and
  `figures/assembly/suppfig11/B_veloae_pca_velocity.pdf`.

## Scientific and display boundaries

The scripts consume saved state/velocity projections and do not refit the
benchmark models, change metric definitions, or recompute benchmark scores.
Panel B uses the same raw-`Ms` PCA axes for state and velocity within each
dataset; it does not claim a shared PCA scale across datasets. Arrow lengths and
angles are display quantities and cannot be read as cross-dataset velocity
magnitudes or high-dimensional cosine values.

The 14 panel identities and display names are resolved through the dataset
registry. Labels, category counts, numeric ticks, source hashes, and final
output hashes are recorded in the manifests. Missing labels or malformed
projection arrays must fail validation rather than fall back to a single color.

## Validation checklist

1. Panel A contains three methods × 14 paired observations for both matrices.
2. Panel B contains 14 unique registered datasets in a five-column, three-row
   plate with one shared key and one orientation marker.
3. Every projection archive has matching `X2`, `V2`, and label lengths and
   finite coordinates.
4. The source hashes in the manifest match the trusted projection archives and
   the materialized A/B assets.
5. The figure remains descriptive: it explains representation concentration
   and shared-axis projection, without changing rankings or claiming causal
   effects on TransCosine.
