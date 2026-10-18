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
