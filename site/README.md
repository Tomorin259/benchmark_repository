# Static results website

[Open the website](https://tomorin259.github.io/benchmark_repository/).
The Vue site displays committed scores, registered images and saved trajectories.
Python export prepares the static data consumed by the Vue pages.

## Build and preview

Requires Python 3.9+ and Node.js 18+. From the repository root:

```bash
cd site
npm ci
npm run build
npm run preview -- --port 4173
```

`npm run build` runs `export_data.py`, then Vite.
`public/data/`, `dist/` and `node_modules/` are ignored; deploy only `dist/`.
The [Pages workflow](../.github/workflows/pages.yml) builds and deploys relevant
pushes to `main`; the GitHub Pages source must be **GitHub Actions**.

For a server preview, serve an existing build and forward its loopback port:

```bash
python3 -m http.server 4174 --bind 127.0.0.1 --directory site/dist
ssh -N -L 127.0.0.1:4174:127.0.0.1:4174 <server>
```

Run the SSH command from your computer, then open `http://127.0.0.1:4174/`.

## Numerical sources

| Display | Source / interpretation |
| --- | --- |
| Cell scores and previews | [Cell input map/results](../evaluation/cell/README.md); five metrics and complete 14 × 20 preview bindings |
| Gene scores | Reference-30 `new_cbdir`; CosSim joined by dataset/method/gene from its separate frozen table |
| Kinetics | Canonical scores, selected Figure 6 plot-data and [expanded-reference results](../evaluation/kinetics_recovery/results/expanded_reference/README.md) are separate scopes |
| Overall | [Committed Figure 7 summaries](../evaluation/summary/README.md), using original-scale scores |
| Ensemble | [Frozen scores/ranks](../ensemble/README.md); both variants remain outside the method registry |
| Cross-level | Committed simulation case rows and selected real-data gene–transition images, with links to the corresponding dataset-level results |
| Diagnostics | Registered veloAE, veloVI and TIvelo materials |
| Downloads | Curated tables, registered final/summary figures and [actual data access](../data/README.md) |

Filters select display rows; summary statistics come from the stored tables.
Nonfinite values become JSON `null` and display as an em dash.
Values display to four decimals and sort at source precision.
Canonical kinetics branch/system weighting, common-valid memberships and
coverage come from their owning tables.
Labelled Dynamo has a distinct identity from the RNA-seq configuration.

Cell Forebrain streamlines show PC1/PC2 of the full PCA30 projection used by
CBDir; BATC uses shared UMAP2. Simulation Cell and its 27 rows are separate
from the real-data summary.
RPE1 turnover strategy uses the stored 250-gene memberships and Macro-F1.

## Regenerate saved visuals

The committed `public/visuals/` bundle lets Pages build without large H5ADs.
Use the scientific plotting environment and materialized inputs only when
regenerating that bundle:

```bash
python site/export_visuals.py --output site/public/visuals \
  --gene-artifact-root /path/to/velocity-artifacts \
  --kinetics-artifact-root /path/to/kinetics-artifacts
```

The destination must be empty; preserve the old bundle outside `public/`
before replacing it. Cell previews use
[render_stream_previews.py](../evaluation/cell/README.md#website-previews).
The visual exporter renders the stored outputs and scores.

RPE1 display curves are independently z-scored with sample SD and keep their
own 0–1 grids; Dynamo has 290 positions. Organoid prediction ratios and Table S4
references keep separate labelled axes and gaps, without interpolation.
TSvelo organoid trajectories match the expanded table's 74 rows and undefined
reasons. Kinetics velocity-field images use the Figure 6 inputs and their
kinetic-recovery analysis roles.

## Content and checks

[content/README.md](content/README.md) identifies reader-facing text and image
sources. The TSV registries own names; the evaluation modules own scientific
definitions. Update the owning sources, then rebuild the exported JSON.

```bash
python site/export_data.py
python -m unittest discover -s tests -p test_site_export.py -v
```

The exporter validates source/display hashes, complete preview bindings and
recorded scope before copying files. `manifest.json` records relative sources,
sizes, schemas, hashes and build identity. Historical visual materialization
and current validation commits remain distinct.
Figures marked `pending_pdf_export` are excluded from final-figure downloads.
Large-file links follow published URIs and the current availability index.

Navigation supports summary/single-gene views, phase/trajectory inspection and
Cell method comparison. Share links retain selected task, scope and display
state; older route aliases remain supported by the existing implementation.
