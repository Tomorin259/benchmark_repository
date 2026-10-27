# Maintainer tools

Run from the repository root. Scientific analysis commands are listed in the
owning method and evaluation modules.

| Entry | Purpose |
| --- | --- |
| `prepare_release.py` | Verify registered artifacts and build release archives using the Python standard library |
| `archive_multivelo_results.py` | Import a complete diagnostic bundle with verified current reference, row coverage and file hashes; material-only import retains the accepted source map |
| `compare_figure6.py` | Compare supplied archived and current kinetics tables and panels; uses NumPy, pandas and Matplotlib |

```bash
python tools/prepare_release.py --help
python tools/archive_multivelo_results.py --help
python tools/compare_figure6.py --help
```

MultiVelo result imports require `results/`, `inputs/` and `provenance.json`
in `--diagnostics-root`. The importer checks the 490 scVelo and 250 MultiVelo
rows and the shared ODE reference identity before copying any result.
Use `--materials-only --materials-dir` for accepted rendering without
changing scientific tables. No directory-name or dated-run fallback is used.
