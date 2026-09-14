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
- **DeepVelo (C16):** select `--evaluation-role cell`, `gene`, or `kinetics`
  and provide an empty `--output-root`. The role determines whether
  `pred_unspliced` is false or true.

For example, a Cell-role DeepVelo run starts with:

```bash
python methods/deepvelo.py \
  --evaluation-role cell \
  --output-root /path/to/empty/deepvelo-cell-run
```

Use `--reset-seed-per-dataset` only when reproducing the historical DeepVelo
launcher RNG scope. The choice is written to `run_config.json`.

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
