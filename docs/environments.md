# Environments

The benchmark uses method-specific Python environments.
Supplementary Table S1 reports package versions and execution devices;
[methods/configurations.tsv](../methods/configurations.tsv) records method
versions and source references.

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

Use the dependencies and software versions in the selected module guide and
run manifest. Reviewed training specifications remain incomplete; available
environment records identify the packages inspected in each recorded runtime.

Source-backed TFvelo, TSvelo and DeepKINET workflows also need the recorded
source tree or release. Historical source revisions are identified where recovered.
