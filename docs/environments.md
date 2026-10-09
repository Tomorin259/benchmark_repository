# Environment evidence and release requirements

Supplementary Table S1 reports the package and selected dependency versions used in the benchmark. A clean-installed plotting lock is available under [`environments/`](environments/). Reviewed method-training Conda specifications and containers remain pending. This document records the requirements that future reviewed environment files must satisfy. Historical environment assignments and import checks remain in the private audit records.

## Plotting environment

`environments/plotting.requirements.txt` pins the seven direct dependencies;
`environments/plotting.lock.txt` records the installed transitive versions.
It was installed in a fresh Python 3.9 virtual environment on Linux. Validation
ran the frozen veloAE PCA summary, all 14 saved projection panels, and all ten
CSV-only Gene summaries (the latter at 60 dpi). These commands neither train
methods nor modify committed figure assets or scientific inputs.

```bash
python3.9 -m venv plotting-env
plotting-env/bin/python -m pip install -r docs/environments/plotting.lock.txt
plotting-env/bin/python supplement/veloae/plot_pca_summary.py --output-dir outputs/veloae-summary
plotting-env/bin/python supplement/veloae/plot_shared_projection.py --output-dir outputs/veloae-projections
plotting-env/bin/python evaluation/gene/plot_gene_assets.py --summaries-only --dpi 60 --output-root outputs/gene-summary
```

The lock covers these plotting workflows only. It does not contain the
AnnData/Scanpy/scVelo dependencies needed for H5AD phase panels or model training.
Cross-runtime pixel or PDF byte identity is not implied.

## Method environments

Existing environment assignments are candidate runtime evidence, not reviewed
locks for the paper's original runs. Package-import smoke tests do not establish
clean installation, launcher execution, GPU compatibility, or scientific
reproduction. Current package versions can differ from the manuscript-reported
versions. Historical assignments and detailed probe results remain private.

Do not publish raw Conda exports: unrelated packages, editable installs and
local paths obscure the method's actual dependency contract. Source-backed
TFvelo, TSvelo and DeepKINET runs also require exact source identity; a
non-Git source tree must not be assigned an invented upstream commit.

## Requirements for a reviewed specification

Each future environment file must:

- identify exactly one or an explicitly shared set of configuration IDs;
- pin the benchmarked method version or upstream commit;
- reproduce the manuscript-reported Python and key dependency versions, or document a validated deviation;
- avoid local filesystem paths and unreviewed editable installs;
- record GPU/CUDA constraints where required;
- install successfully in a clean environment;
- pass the corresponding runner smoke test;
- preserve third-party license notices without vendoring unapproved assets.

The audited Conda paths are provenance evidence only and will not appear in the final public execution interface.
