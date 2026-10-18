# veloAE representation diagnostics

Source material for Supplementary Figure 11:
A compares state/velocity PC1 explained variance, and B displays each of the 14
datasets in its own state-PCA coordinate system. The final PDF is linked in the
[figure map](../../docs/figure_map.md).

## Inputs and reproduction

The compact inputs are committed:
[inputs/manifest.json](inputs/manifest.json) identifies the 84-row PCA
comparison table and 14 projection NPZs.
Use Python, NumPy, Matplotlib and Pillow, or the checked
[plotting environment](../../docs/environments.md).

```bash
python supplement/veloae/plot_pca_summary.py --output-dir outputs/veloae
python supplement/veloae/plot_shared_projection.py --output-dir outputs/veloae
python supplement/veloae/export_ppt_assets.py --input-dir outputs/veloae --output-dir outputs/veloae/ppt_assets_final
```

Run from the repository root with fresh destinations.
The first two commands render A and the B plate; the last exports letter-free
materials for manual layout. PDF companions depend on the converter.
The renderers use saved projections, without PCA refitting or metric changes.
NPZ label arrays contain pickled objects and require trusted input files.

## Interpretation and display

A compares veloAE, scVelo stochastic and TIvelo, with 14 paired datasets per
matrix and method. PCA is fitted independently to each state and velocity
matrix. Filled/open points denote raw/standardized PCA; lines pair datasets
and bars show cross-dataset means.
The raw/standardized veloAE PC1 fractions are 95.7%/70.0% for state and
84.8%/52.7% for velocity.

B reads saved `X2`, `V2`, labels and PCA axes. State coordinates and velocities
use the same PC1–PC2 axes fitted to the dataset's raw `Ms`. State coordinates
are centered; velocity is projected without centering. Cells are colored by
dataset annotations, and dataset names follow the registry.
Arrows show spatial-bin mean projected velocities, scaled and length-capped
for display. Axis ranges and arrow scaling are dataset-specific.
The 300-dpi dataset panels form a row-major 5 × 3 plate, without interpolation,
with a shared cluster key and coordinate-orientation marker.
The marker is not a velocity scale bar, and unequal axis scaling means visual
angles cannot be read as high-dimensional cosines.
Missing labels or malformed arrays stop rendering.

A is exported at 180 mm/600 dpi; the B export retains its native 300-dpi panels.
[results/ppt_assets_final/](results/ppt_assets_final/) contains the accepted
source materials and manifest.
The final two-panel PDF uses manual author assembly. Independent reconstruction
of the whole-page layout remains incomplete.
Superseded layouts are retained at
[revision 51edf35](https://github.com/Tomorin259/benchmark_repository/tree/51edf35695e2b1db1449596f2878c1be3ee466c3/supplement/veloae/archive/non_online_20260927).
