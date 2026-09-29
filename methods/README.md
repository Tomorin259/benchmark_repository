# Method reproduction

This directory owns the direct method launchers used by the benchmark. The
authoritative list of configurations, paper-facing names, versions, scientific
options, and required result keys is
[`configurations.tsv`](configurations.tsv). It contains 20 registered
configurations; ensemble comparators are not method entries.

## Entry points

Each launcher keeps its dataset map and scientific configuration visible near
the top of the script. C01--C03 share `scvelo.py`, and C05--C06 share
`unitvelo.py`; every other registered configuration has its own script in this
directory. Edit the top-level dataset list for a focused run, or use a script's
command-line dataset option where one is provided.

The launchers follow the corresponding upstream workflow, validate the input
and result fields used by evaluation, and write method-native H5AD files under
the ignored `outputs/` directory. They do not define benchmark metrics. The
evaluation modules own the exact mappings from released method results to
paper analyses.

Direct launchers remove their own script directory from Python's import search
path before importing method packages, so sibling files such as `scvelo.py`
cannot shadow the installed package.

Large H5AD files, checkpoints, logs, upstream source trees, and Conda
installations are not stored here. Released external files are identified by
logical ID and checksum in [`../artifacts.tsv`](../artifacts.tsv).

## Method-specific requirements

- **cell2fate and DeepKINET:** these launchers select raw-count inputs where
  required, align cells and genes, verify integer counts, and prevent a second
  preprocessing pass. The dataset-specific rules remain in the scripts.
- **DeepKINET (C13):** the formal benchmark configuration uses seed 123 and
  deterministic PyTorch/CUDA settings. The audited source revision is recorded
  in `configurations.tsv`.
- **TFvelo (C10):** provide the official source checkout with `--source-root`
  or `TFVELO_ROOT`. The launcher stages the relative paths required by the
  upstream demos.
- **TSvelo (C14):** provide the official source checkout with `--source-root`
  or `TSVELO_ROOT`. The launcher stages the upstream path layout and preserves
  the source workflow's no-global-seed behavior.
- **veloVI (C12):** select `--evaluation-role cell` or `gene` and provide
  an empty `--output-root`. Cell retains the package-default training call,
  then divides spliced velocity and kinetic rates by `20 / max(latent_time)`
  separately for each gene. It also saves the tutorial's latent-time and fit
  fields. Gene trains for up to 1,000 epochs with batch size 64 on GPU or 32
  on CPU and saves unscaled `velocity` plus `velocity_u`. Both use 25
  posterior samples. These are distinct training runs; the role selects the
  entire workflow. Selected settings are saved in H5AD `uns["benchmark_run"]`.
- **DeepVelo (C16):** select `--evaluation-role cell`, `gene`, or `kinetics`
  and provide an empty `--output-root`. Cell sets `pred_unspliced=false`;
  Gene and Kinetics set it to true. The training head and required result
  fields follow that choice; selected settings are saved in `run_config.json`.

For example, a Cell-role DeepVelo run starts with:

```bash
python methods/deepvelo.py \
  --evaluation-role cell \
  --output-root /path/to/empty/deepvelo-cell-run
```

Use `--reset-seed-per-dataset` only when reproducing the historical DeepVelo
launcher RNG scope. The choice is written to `run_config.json`.

A focused Gene-role veloVI run uses the same explicit role/output boundary:

```bash
python methods/velovi.py \
  --evaluation-role gene \
  --dataset Pancreas \
  --output-root /path/to/empty/velovi-gene-run
```

Use `--evaluation-role cell` and a separate empty output directory for Cell.
veloVI also accepts `--data-root`; DeepVelo uses its top-level `DATA_ROOT`.
The method launchers do not replace the frozen files selected by evaluation.

## Contract boundaries

- Method names and scientific settings belong in `configurations.tsv`.
- Dataset-specific execution details belong in the launcher that consumes
  them.
- Published result selection belongs in the relevant evaluation module.
- Reproducibility policy shared across methods belongs in
  [`../docs/reproducibility.md`](../docs/reproducibility.md).
- Ensemble construction and results belong in [`../ensemble/`](../ensemble/).

General scientific-script design rules are defined in
[`../REPOSITORY_DESIGN.md`](../REPOSITORY_DESIGN.md).
