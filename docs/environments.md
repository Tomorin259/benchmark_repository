# Environments

The benchmark uses method-specific Python environments.
Supplementary Table S1 reports package versions and execution devices;
[methods/configurations.tsv](../methods/configurations.tsv) records method
versions and source references.

Create environments with your preferred package manager and names. Commands
use `python` from the environment selected for that task. When a command
accepts a Python path, supply the executable in your own installation; for
example, `/path/to/cell-evaluation/bin/python`.

## CSV-only plotting

The committed plotting lock was installed in a fresh Python 3.9 environment on
Linux and checked with the veloAE diagnostics and CSV-only Gene summaries:

```bash
python3.9 -m venv plotting-env
plotting-env/bin/python -m pip install -r docs/environments/plotting.lock.txt
plotting-env/bin/python supplement/veloae/plot_pca_summary.py --output-dir outputs/veloae-summary
plotting-env/bin/python supplement/veloae/plot_shared_projection.py --output-dir outputs/veloae-projections
plotting-env/bin/python evaluation/gene/plot_gene_assets.py --summaries-only --dpi 60 --output-root outputs/gene-summary
```

[plotting.requirements.txt](environments/plotting.requirements.txt) lists direct
dependencies; [plotting.lock.txt](environments/plotting.lock.txt) pins the tested
environment for CSV-only plotting. H5AD phase plots and model training require
the AnnData/Scanpy/scVelo stack and method-specific dependencies.

## Evaluation and method training

### Evaluation

| Task | Python and key packages | Version evidence / entry |
| --- | --- | --- |
| Cell preparation and local metrics | Python 3.7.12; scVelo 0.2.5, UniTVelo 0.2.5.2, Scanpy 1.9.3, AnnData 0.8.0, NumPy 1.21.6, SciPy 1.7.3, pandas 1.3.5 | Retained worker package metadata, inspected 2026-10-19; [Cell instructions](../evaluation/cell/README.md#recalculate-metrics) |
| Official TIvelo Cell metrics | Python 3.9.23; TIvelo 0.1.4, scVelo 0.3.1, Scanpy 1.10.3, AnnData 0.10.8, NumPy 1.26.4, SciPy 1.13.1, pandas 2.3.1 | Retained metric-worker package metadata, inspected 2026-10-19; [official scorer](../evaluation/cell/score_tivelo_metrics.py) |
| Gene coordinate restoration and score replay | Python 3.9.24; scVelo 0.3.1, Scanpy 1.10.3, NumPy 1.26.4, SciPy 1.13.1, scikit-learn 1.1.3; AnnData and pandas for H5AD/CSV input | [Recorded replay environment](../evaluation/gene/results/unit_restoration/verification.json); [Gene instructions](../evaluation/gene/README.md) |
| Kinetic-parameter evaluation | Python 3.9+; recorded packages NumPy 1.26.4, pandas 2.3.1, SciPy 1.13.1, AnnData 0.10.8; h5py for H5AD input | [Evaluation manifest](../evaluation/kinetics_recovery/results/central_results/manifest.json), `evaluator_versions`; [kinetics instructions](../evaluation/kinetics_recovery/README.md) |
| Signed-transition ensemble | scVelo 0.3.4 with AnnData, Scanpy, NumPy, pandas, SciPy, scikit-learn and UniTVelo | [Protocol](../ensemble/config.json) and [ensemble instructions](../ensemble/README.md) |

The Cell controller uses the standard library and starts its two workers
through the Python executables supplied to `--indicators-python` and
`--tivelo-python`. Use Python 3.9 or later for the controller. These options
identify worker roles, so environment names and installation paths are arbitrary.

Package metadata is a record of an installed combination. The plotting lock
above has a separate clean-install check; evaluation validation is described
in the linked run records. The original Gene reference-30 provenance records
scVelo 0.2.5 / Scanpy 1.9.3; current unit-restoration replay has its own record.

### Method execution

Install each method at the version in the method registry and follow its
upstream installation instructions, including framework and accelerator
dependencies. The public launchers import the following method interfaces:

| Launchers | Method-specific dependency |
| --- | --- |
| `scvelo.py`, `unitvelo.py`, `multivelo.py`, `tivelo.py`, `velovi.py`, `cell2fate.py`, `celldancer.py`, `sdevelo.py` | Corresponding registered package, together with its Scanpy/scVelo and model-framework dependencies |
| `dynamo.py` | Distribution `dynamo-release`, imported as `dynamo` |
| `graphvelo_scv.py`, `graphvelo_dynamo.py` | `graphvelo`; upstream scVelo or Dynamo stack according to the configuration |
| `veloae.py` | Distribution `veloae`, providing the `veloproj` import |
| `deepvelo.py` | `deepvelo`, PyTorch and scVelo 0.2.5; the launcher checks the scVelo version |
| `deepvelo2022.py` | TensorFlow/Keras reconstruction code in the launcher, plus Scanpy/scVelo; C18 is a distinct configuration |
| `tfvelo.py`, `tsvelo.py`, `deepkinet.py` | Upstream source-backed software; see [source requirements](../methods/README.md#source-backed-methods) |

Complete historical training environments are not recorded for every
configuration. The [method guide](../methods/README.md#source-backed-methods)
provides fixed TFvelo/TSvelo sources and DeepKINET installation patches;
[source fingerprints](../methods/source_manifest.json) identify the reconstructed
files. Original training records did not capture these fingerprints. C18's
full historical software identity remains as recorded in the registry.
Use the selected module's source and run record when matching a frozen output.

The retained DeepKINET installation reports Python 3.9.23, deepkinet 0.0.1,
PyTorch 2.1.0+cu118, scVelo 0.2.5, Scanpy 1.9.6, AnnData 0.10.3 and
NumPy 1.23.5. This is installed-package evidence, rather than a clean-install
lock. Install the patched package in your own environment; the public launcher
requires scVelo 0.2.5 and records per-dataset checkpoints and randomness.

### Hardware and rendering tools

The registry identifies the execution devices used for each configuration.
Choose a PyTorch/TensorFlow build and CUDA runtime compatible with your own
GPU and driver. `CUDA_VISIBLE_DEVICES` selects the visible devices for
DeepVelo, DeepVelo_scv and SDEvelo. TSvelo exposes `--cuda`; UniTVelo and
veloAE use the `GPU_INDEX` constant described in the method guide. Thread
counts are exposed by CLI options or top-level constants in each launcher.

Numerical plots use the plotting dependencies above, with H5AD dependencies
where required. Kinetics panels also use `matplotlib-venn` (recorded 1.1.2).
TIvelo A/B PDF assembly needs `pypdf`, `reportlab` and a user-supplied Arial TTF;
MultiVelo A/B assembly needs `pypdf` and Pillow. The veloAE projection PDF
export uses the external `rsvg-convert` executable from librsvg. Editable
whole-figure assemblies use the PowerPoint/WPS export steps in the
[figure map](figure_map.md).
