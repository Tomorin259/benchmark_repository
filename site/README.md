# Static result website

[Open the benchmark website](https://tomorin259.github.io/benchmark_repository/).

This Vue website reads existing benchmark tables and registry metadata. It
does not train methods, run evaluators, modify result tables, or need a
database/API server. The website is a display projection, not a new result
authority.

## Build and preview

From the repository root, using Python 3.9+ and Node.js 18+:

```bash
cd site
npm ci
npm run build
npm run preview -- --port 4173
```

Open the printed local URL. The build runs `export_data.py` before Vite and
copies only the selected compact source files into `public/data/`. Generated
data, `dist/`, and `node_modules/` are ignored. Deploy only `dist/`; the
relative asset paths and hash navigation support both a root domain and a
project subdirectory. `.github/workflows/pages.yml` validates and builds the
site on relevant pushes to `main`, then publishes `dist/` to GitHub Pages.
The repository's Pages publishing source must be **GitHub Actions**.

On a server without Node.js, copy the built `site/dist/` and serve the static
files from the repository root:

```bash
python3 -m http.server 4174 --bind 127.0.0.1 --directory site/dist
```

To view that server on your computer, open a separate terminal and forward
the same loopback port through your normal SSH connection:

```bash
ssh -N -L 127.0.0.1:4174:127.0.0.1:4174 <server>
```

Then open `http://127.0.0.1:4174/`. No public inbound HTTP port is needed.

The exporter uses the Python standard library and can also be run separately:

```bash
python site/export_data.py
python -m unittest discover -s tests -p test_site_export.py -v
```

## Display contract

- Cell scores are joined to registry identities through the frozen
  `evaluation/cell/inputs.tsv`, with all primary source values preserved.
- Gene CBDir uses `new_cbdir` from the reference-30 table. Simulation CosSim
  is joined by dataset/method/gene from the separate frozen compact table;
  that table's CBDir does not replace the reference-30 score.
- Kinetics offers expanded-reference genes, canonical scores and Figure 6
  selected plot-data as separate scopes. Canonical component means and counts
  come from `component_means.csv`. Expanded scores restrict canonical rows
  to `gene_scope.csv` and add the isolated TSvelo organoid alpha table;
  means and coverage come directly from `coverage_summary.csv`. Common-valid
  comparisons use the committed summary and exact membership, without
  recomputing their intersection when visitors filter methods. Organoid
  branches remain separate; canonical equal-branch system summaries stay in
  the owning results module. Expanded scores do not replace Figure 7 aggregation.
- The expanded manifest's canonical/output hashes are verified before export.
  Curated result CSVs and registered summary figures are downloadable;
  provenance remains in the source manifest and owning module guides.
  Source figures are not labelled as a final paper composite based only on
  their filenames. Their manifest owns the current registration boundary.
- Cell results show all five metrics with per-method stream previews and a 2–4 method comparison. Gene and Kinetics provide explicit Summary / Single gene selectors. Results default to a table, with a plot toggle. Kinetic method order follows
  `PAPER_METHOD_ORDER` in the Figure 6 renderer; score sorting is an explicit
  visitor choice. No new method ranking or scientific mean is calculated.
- Labelled-data Dynamo has its own display identity; it is not the C07
  label-free configuration.
- Missing or nonfinite scores become JSON `null` and are displayed as `—`.
  Values are displayed to four decimal places and sorted at source precision.
  Method filtering changes the displayed rows, not any stored score or mean.
- `manifest.json` records source commit, relative paths, row/schema information,
  sizes and SHA-256. Original compact downloads are copied byte-for-byte.
- Large-file downloads are shown only for actual HTTPS URIs from the external
  artifact registry. The input release manifest and expanded reference-gene
  tables are published as lightweight release metadata; they do not imply
  that source H5AD files or full method outputs are publicly downloadable.
  Upstream data references are not relabelled as frozen benchmark input
  downloads.

## Single-gene visual assets

In Results, choose **Single gene**, select a gene and then **Phase portraits** (Gene-level) or
**Trajectories** (Kinetic-parameter recovery). Mean views remain score views.
Gene and Kinetics method filters apply to their score rows and visual. Cell
results use Dataset and Preview method for the single streamline preview.
Compare methods opens an optional 2–4 method view; closing it restores one method.
Every preview shows all five Cell scores. Table headers control sorting, and
Plot metric appears only beside the score plot; neither changes stream images.
Trajectories
support pointer inspection, drag or numeric range zoom, and displayed-series
method filtering. Phase portraits support enlargement. These single-gene
displays have no download controls. Downloads exposes only the curated result
tables and registered summary figures; H5AD links follow `artifacts.tsv` and
the H5AD members in `data/release_map.tsv` once public URIs are provided.

Generate these static displays once in the scientific plotting environment
(NumPy, pandas, SciPy, AnnData, Matplotlib, seaborn, scVelo and scikit-learn):

```bash
python site/export_visuals.py --output site/public/visuals \
  --gene-artifact-root /path/to/velocity-artifacts \
  --kinetics-artifact-root /path/to/kinetics-artifacts
cd site
npm run build
```

Both roots use the owning modules' relative input contracts. Use an empty
output directory; retain an old bundle outside `public/` before replacing it.
No method training is run. The exporter verifies the saved kinetic matrix
hashes, uses the current phase renderer with native `Ms`/`Mu` and velocity
layers, and records absent genes separately from undefined scores. Hybrid
portraits are generated from the existing frozen simulation outputs.

RPE1 applies the Figure 6 display z-score independently to each saved curve
using sample SD. Each curve retains its own 0–1 relative-position grid;
Dynamo retains 290 positions. Original bin IDs remain in the static data; the pointer tooltip displays relative position.
Organoid retains the processed early-50-bin prediction ratios and the saved
Table S4 reference values/direction on separate labelled y-axes. No numerical
conversion is applied to make these different scales coincide. Missing
samples form gaps; no interpolation is applied.
TSvelo organoid alpha is processed from the hash-verified lineage inputs with
the existing evaluator's functions; all 74 metric rows and undefined reasons
must match the committed expanded table before display export. These generated
trajectories do not add canonical or selected Figure 6 scores.

`public/visuals/` is a committed display bundle so Pages can build without
the large H5AD inputs. Its manifest records the scientific inputs and renderer
hashes; it is not a replacement for the source result tables. Regenerate it
explicitly when its recorded sources change. The build checks all recorded
source and display-file hashes, and the Pages workflow rejects a missing
bundle. `public/data/`, `dist/`, and `node_modules/` remain ignored. Deploy
the assembled `dist/` to serve everything statically; browsers never read
H5AD or run Python.

Scientific website copy follows the manuscript's evaluation terminology.
Do not add promotional claims or decorative result interpretations. Verify
manuscript metadata and final figure captions against the authoritative
manuscript before adding them.

## Source ownership

Method and dataset display labels come from their TSV registries. Existing
result identifiers are translated to these identities only at export. The
scientific contracts remain in each evaluation module's README, linked from
the website. Do not edit generated JSON to change scientific results.

Update an owning result table or registry through its established workflow,
then rebuild and check the source manifest. For a release, freeze the source
commit and public artifact locations before enabling website publication.

## Cell previews and reader-facing navigation

Cell previews are exported from `evaluation/cell/streams/preview_manifest.json`.
All 280 keys must match `evaluation/cell/inputs.tsv`, including source H5AD hashes.
Image, input, renderer and palette drift stops the build. `public/data/previews/`
is generated; the Cell module owns the source images. Only the currently opened
method images are fetched by the browser. See `evaluation/cell/README.md` for
materialization. The retained Figure 1 PDF is rendered as the static overview
preview with a source/display checksum record.

Cell score tables omit repeated dataset cell counts. Only partial alignment
coverage is explained in the method preview. Gene and Kinetics retain their
valid-gene or paired-position counts; CosSim does not reuse CBDir sample counts.
The share link records task, dataset, parameter, branch, scope, gene, metric,
display mode, and the selected Cell method or comparison set. Older `streams`
links remain readable for compatibility.

Overview, Results, Ensemble, Diagnostics, Datasets & Methods and Downloads are
the main navigation. Documentation is available from the footer. Simulation
inputs are grouped in the dataset directory; previous hash routes remain valid.
SHA-256 values and commit identities stay in machine records and the repository.
Normal browsing shows objects, selectors, plots, scores and concise scope labels.

The Ensemble page binds both frozen variants to their registered stream images
on all 14 Cell datasets. It loads the selected dataset images, supports
enlargement and original-image download, and preserves the independent rankings.
Neither variant enters the 20-configuration method registry.

Diagnostics exposes existing veloAE PCA summaries/projections, veloVI coordinate
conditions and TIvelo root/ordering controls. It copies scores and plotted data
without running a method or evaluator. The published diagnostic images are
pinned by `content/diagnostic_assets.json` and checked against owning manifests
and frozen input identities. This display binding does not assign paper panels.
The veloAE summary uses the verified SVG; its PNG is not used as a substitute.
TIvelo simulation refit gene counts and paired bone-marrow optimizer seeds stay
separate from primary benchmark scores. Bone-marrow maps retain archived t-SNE.
These diagnostic results never replace method-ranking inputs.

Additional analysis CSVs are copied byte-for-byte. All existing result, release,
figure and audit downloads remain available; metadata/audits form one technical
resources section. There is no claim that the unpublished H5AD artifacts are
downloadable merely because their upstream sources are linked.

Scientific context uses registry-backed names and the documented manuscript read.
Model/dataset facts and metric terms stay in `content/website_content.json`;
`content/README.md` explains their authority. The build records the historical
visual materialization commit separately from its validation commit. Source
and display hashes must continue to match; final paper figure mapping is separate.
