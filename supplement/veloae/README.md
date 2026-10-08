# veloAE supplementary diagnostics

This directory contains source renderers and final A/B source assets for online
Supplementary Figure 11. Plotting does not refit PCA, modify the saved
state/velocity projections, or recompute benchmark scores. The exact compact inputs are included under `inputs/`.

## Inputs and reproduction

`plot_pca_summary.py` reads the frozen `inputs/pca_comparison.csv` (three methods,
14 paired observations per matrix and method). The accepted paired-dot layout
is retained; PNG output is 180 mm wide at 600 dpi, with PDF/SVG companions.

```sh
python supplement/veloae/plot_pca_summary.py --output-dir outputs/veloae
python supplement/veloae/plot_shared_projection.py --output-dir outputs/veloae
python supplement/veloae/export_ppt_assets.py --input-dir outputs/veloae --output-dir outputs/veloae/ppt_assets_final
```

`inputs/manifest.json` records the original logical source and SHA-256 for the
84-row PCA table and all 14 projection archives. Their bytes match the frozen
input hashes in the retained output manifests. The mean raw/standardized PC1
fractions for veloAE are 95.7%/70.0% (state) and 84.8%/52.7% (velocity).

`plot_shared_projection.py` reads the 14 trusted diagnostic NPZ archives, using
their saved `X2`, `V2`, per-cell labels and PCA axes. It renders each dataset at
300 dpi, then assembles the PNGs without interpolation. Display names are joined
to `data/datasets.tsv`; no independent paper-name dictionary is maintained.
Numeric ticks show the actual projected coordinates. Panels follow the original
row-major five-column, three-row order, with a compact key in the final slot.
One shared PC1/PC2 schematic sits at the lower left of the complete plate;
it indicates positive coordinate directions, not a velocity scale bar. A compact
shared key occupies the unused bottom-right slot.

The environment requires Python, NumPy, Matplotlib and Pillow. Pickled label
arrays in these existing NPZ archives require trusted inputs. No silent
single-color fallback is permitted. Manifests record source hashes and panel
parameters. Check the manifests, all 14 label/category counts, numeric ticks,
and the finished figures before use. Do not infer high-dimensional cosines
from angles in a panel with unequal x/y display scaling.

`export_ppt_assets.py` is the handoff for final manual layout. It exports two
letter-free assets: the native side-by-side `pca_summary` and the complete
`veloae_pca_velocity` plate, each as PNG/SVG and PDF where the converter is
available. The former is 600 dpi; the latter preserves its native 300-dpi
dataset panels. It does not create or modify a PowerPoint.

The online `merged.tex` assigns the paired PCA summary to panel A and the
14-dataset shared-axis projection plate to panel B. The final manual page
composite is not part of this source bundle.
