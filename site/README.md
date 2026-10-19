# Website source

Source for the [published benchmark results website](https://tomorin259.github.io/benchmark_repository/).
The site displays the saved scores and figure materials from the evaluation modules.

Requires Python 3.9+ and Node.js 18+. From the repository root:

```bash
cd site
npm ci
npm run build
npm run preview -- --port 4173
```

`npm run build` exports static data with `export_data.py`, then builds the Vue
site with Vite. The output is `site/dist/`.

Result downloads use the verified `data/result_package_manifest.tsv` and
`data/result_member_manifest.tsv`. The exporter exposes `result_packages`
and `result_members` and copies both manifests into Technical resources at
the end of the Downloads page.
All package URLs are pinned to the `v1.0.0` release.
Archive hashes and original H5AD/NPZ member hashes remain separate; no
training or metric calculation is performed by the exporter.

Downloads starts with the figure/table-indexed Source data packages, followed
by Benchmark input files and Full method outputs. The original kinetic-reference
access guide is linked from the input section. Result tables and summary figures
remain in their owning repository modules; the download page uses Source Data as
their consolidated entry. Existing individual download files and URLs remain
available, including downloads beside results and diagnostic figures.

The Full method outputs section reads `result_packages`, with dataset,
registered-method and output-type filters, 25 packages per page, explicit
TAR.GZ/NPZ download buttons and package sizes. Checksums remain in the
downloadable manifests. Native H5AD
artifact URIs are not used as archive download links.

The Source data section reads `data/source_data_manifest.tsv`; its ZIP links
are pinned to the published Hugging Face content commit.
