# RNA velocity benchmark

Code, compact results and figure sources for a benchmark of 20 RNA velocity
method configurations across 14 real datasets and three simulated scenarios.
The study evaluates cell-level trajectories, gene-level velocity directions
and kinetic-parameter recovery against two scEU-seq reference systems.

[Browse results](https://tomorin259.github.io/benchmark_repository/) ·
[Download data](data/README.md) ·
[Find a paper figure](docs/figure_map.md)

## Start here

| Goal | Entry |
| --- | --- |
| Inspect scores without running methods | [Cell](evaluation/cell/results/cell_metrics.csv), [Gene](evaluation/gene/results/README.md), [kinetic parameters](evaluation/kinetics_recovery/README.md#results-and-scope) and [Figure 7 summaries](evaluation/summary/README.md) |
| Recalculate metrics from saved method outputs | [Cell evaluation](evaluation/cell/README.md), [Gene evaluation](evaluation/gene/README.md), [kinetics recovery](evaluation/kinetics_recovery/README.md) |
| Run a velocity method | [Method launchers](methods/README.md) and [environment requirements](docs/environments.md) |
| Reproduce plots or assemble paper figures | [Figure map](docs/figure_map.md) |
| Generate a new simulation realization | [Simulation generator](evaluation/gene/simulation/README.md) |
| Inspect supplementary controls | [veloAE](supplement/veloae/README.md), [veloVI](supplement/velovi/README.md), [TIvelo](supplement/tivelo_gene/README.md), [MultiVelo](supplement/multivelo/README.md) and [gene–transition cases](supplement/gene-kinetic/README.md) |

For example, render Gene summary plots from the committed CSVs after installing
the [plotting environment](docs/environments.md):

```bash
python evaluation/gene/plot_gene_assets.py --summaries-only \
  --output-root outputs/gene-summary
```

Run commands from the repository root and use a fresh output directory.

## Data and reproducibility

- [Canonical inputs](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-inputs)
  provide 25 original-format files at distribution `v1.0.0`.
  Three Battich reference tables are obtained separately; see
  [reference access and format requirements](data/UPSTREAM_REFERENCES.md).
- [Method outputs and evaluation assets](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results)
  are available as `v1.0.0`: 339 download objects totaling 149.58 GB,
  including the MultiVelo diagnostic outputs. Select packages from the
  [download manifest](data/result_package_manifest.tsv), map original H5AD/NPZ
  files through the [member manifest](data/result_member_manifest.tsv), and follow the
  [download and extraction guide](data/RESULT_DOWNLOADS.md).
- Committed tables support inspection and several CSV-only plotting workflows.
  Recomputing metrics requires the saved H5AD outputs and reference inputs;
  retraining additionally requires method-specific software. See the checked
  [plotting environment and method-training requirements](docs/environments.md).
- The benchmark uses the registered simulation H5ADs. The
  [generator](evaluation/gene/simulation/README.md) creates additional simulated datasets.
- Some final figures use author-edited PowerPoint and manual PDF export.
  Their sources and assembly steps are listed in the [figure map](docs/figure_map.md).

## Repository layout

| Directory | Contents |
| --- | --- |
| `data/` | Dataset identities, downloads and input layouts |
| `methods/` | Registered configurations and method launchers |
| `evaluation/` | Cell, Gene and kinetics metrics, results and plots |
| `ensemble/` | Mean-4 and GeoMedian-4 post-hoc comparisons |
| `supplement/` | Method diagnostics and additional analyses |
| `figures/` | Paper composites and editable assembly sources |
| `site/` | Static results website |
| `tests/` | Scientific and file-format regression checks |

Names come from [methods/configurations.tsv](methods/configurations.tsv) and
[data/datasets.tsv](data/datasets.tsv). [artifacts.tsv](artifacts.tsv) identifies
large files and retained assets.
See the [reproduction guide](docs/reproducibility.md) for required inputs,
software environments and evaluation commands.
