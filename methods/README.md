# Run a method

[configurations.tsv](configurations.tsv) lists the 20 configurations, versions,
source references, scientific options and required result keys.
Use [canonical inputs](../data/README.md) and the method's
[software environment](../docs/environments.md).

Randomness handling follows the corresponding upstream workflows, with
configuration-specific settings listed in the registry's `randomness` column.

Run commands from the repository root. Launchers write native H5ADs under an
explicit output root or the ignored `outputs/` default. Most retain scientific
settings and dataset lists near the top of the script; CLI dataset selection
is available where described below. Evaluation selects its frozen result files
through separate input maps.

## Modes and example commands

| Method | Selection |
| --- | --- |
| scVelo, C01–C03 | `--mode deterministic`, `stochastic` or `dynamical`; repeat for several modes; default `MODES_TO_RUN` |
| UniTVelo, C05–C06 | `--fit-option 1` or `2`; repeat for both; default `FIT_OPTIONS_TO_RUN` |
| TFvelo / TSvelo | `--dataset <name>` or `all`; `--check-only` checks the selected source/input layout |
| veloVI | `--evaluation-role cell` or `gene` |
| DeepVelo | `--evaluation-role cell`, `gene` or `kinetics` |

```bash
python methods/scvelo.py --mode dynamical --dataset Pancreas \
  --output-root /path/to/scvelo-run
python methods/unitvelo.py --fit-option 2 --dataset Pancreas \
  --output-root /path/to/unitvelo-run
python methods/velovi.py --evaluation-role gene --dataset Pancreas \
  --output-root /path/to/empty/velovi-gene-run
python methods/deepvelo.py --evaluation-role cell \
  --output-root /path/to/empty/deepvelo-cell-run
```

scVelo and UniTVelo accept repeated `--dataset` selections, plus `--data-root`
and `--output-root`; modes use separate directories and reject existing H5ADs.
scVelo dynamical results use `scvelo_dynamic/`.
veloVI also accepts `--data-root`; DeepVelo uses its top-level `DATA_ROOT`.
Use separate empty roots for distinct evaluation roles.

## Method-specific inputs and workflows

- **cell2fate and DeepKINET:** raw-count inputs, integer-count checks, cell/gene
  alignment and protection against repeated preprocessing follow the launcher.
- **DeepKINET (C13):** the benchmark uses seed 123 with deterministic
  PyTorch/CUDA settings. Source identity is recorded in the registry.
- **TFvelo (C10):** provide its official checkout with `--source-root` or
  `TFVELO_ROOT`.
- **TSvelo (C14):** provide its official checkout with `--source-root` or
  `TSVELO_ROOT`. The launcher follows the upstream workflow.
  Both source-backed launchers stage the upstream demo's relative input layout.
- **veloVI (C12):** Cell calls package-default training and divides spliced
  velocity and rates by the per-gene `20 / max(latent_time)` factor.
  Gene trains for up to 1,000 epochs, with batch size 64 on GPU or 32 on CPU,
  saves unscaled spliced velocity and extracts `velocity_u`.
  Both use 25 posterior samples. The evaluation role selects the corresponding
  training and output workflow.
- **DeepVelo (C16):** Cell sets `pred_unspliced=false`; Gene and Kinetics set
  it to true, changing the training head and required outputs.
  `--reset-seed-per-dataset` reproduces the documented historical RNG scope.

scVelo, UniTVelo and veloVI record their selected settings in
`uns["benchmark_run"]`; DeepVelo writes `run_config.json`.
Entry scripts avoid sibling-file package shadowing.

## Use the results

[Cell](../evaluation/cell/README.md), [Gene](../evaluation/gene/README.md) and
[kinetics](../evaluation/kinetics_recovery/README.md) document the exact result
family, layer and alignment consumed by each task. Ensembles are
[post-hoc comparisons](../ensemble/README.md) outside the method registry.
See [reproduction routes](../docs/reproducibility.md) for data/version
identities and [repository design](../REPOSITORY_DESIGN.md) before changing code.
