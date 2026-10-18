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
and `result_members` and copies both manifests into Technical resources.
Archive hashes and original H5AD/NPZ member hashes remain separate; no
training or metric calculation is performed by the exporter.
