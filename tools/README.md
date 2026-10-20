# Maintainer tools

Run from the repository root. Scientific analysis commands are listed in the
owning method and evaluation modules.

| Entry | Purpose |
| --- | --- |
| `prepare_release.py` | Verify registered artifacts and build release archives using the Python standard library |
| `archive_multivelo_results.py` | Import identified MultiVelo runs into compact result tables and plotting inputs; uses the MultiVelo analysis environment |
| `compare_figure6.py` | Compare supplied archived and current kinetics tables and panels; uses NumPy, pandas and Matplotlib |

```bash
python tools/prepare_release.py --help
python tools/archive_multivelo_results.py --help
python tools/compare_figure6.py --help
```
