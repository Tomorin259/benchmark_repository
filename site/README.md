# Website source

Source for the [published benchmark results website](https://tomorin259.github.io/benchmark_repository/).
The site displays the saved scores and figure materials from the evaluation modules.

Diagnostics includes MultiVelo's RNA-only phase-assignment controls at
`#diagnostics?method=multivelo`: the retained phase portraits and paired-score
figure from Supplementary Figure 15, plus the four-condition summary in
Supplementary Table 9. Download links use the published Source Data files.

Requires Python 3.9+ and Node.js 18+. From the repository root:

```bash
cd site
npm ci
npm run build
npm run preview -- --port 4173
```

`npm run build` exports static data with `export_data.py`, then builds the Vue
site with Vite. The output is `site/dist/`.

Downloads presents figure/table Source Data, benchmark inputs and method outputs.
The exporter reads their manifests in `data/`, including URLs, sizes and hashes.
Method outputs can be filtered by dataset, method and output type, with 25
packages per page. Source Data includes main figures, supplementary figures,
supplementary tables and simulation velocity references, plus a complete ZIP.

## Annotated evaluation transitions

Single-dataset Cell results and Gene CBDir results include a reference diagram
on the saved UMAP coordinates (PC1/PC2 for Forebrain). Arrows come from the
corresponding evaluation configuration; Retina has separate Cell and Gene
versions. The diagrams share the reference cell labels and published palette.
PNG/SVG downloads and the complete transition list are available beside each
image.

Use the [plotting environment](../docs/environments.md#csv-only-plotting)
with h5py. Run from the repository root:

```bash
python -m pip install h5py
python site/render_transition_references.py \
  --artifact-root /path/to/velocity-artifacts \
  --output site/public/transitions
```

Here the artifact root is the velocity directory itself, as for Gene evaluation.
The renderer verifies registered input hashes and records the source files,
coordinate identity, task-specific edges and output hashes in
`site/public/transitions/index.json`. The normal website build checks the
committed diagrams without reading external H5ADs or recalculating scores.
