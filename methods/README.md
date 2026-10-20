# Run a method

[configurations.tsv](configurations.tsv) lists the 20 configurations, versions,
source references, scientific options and required result keys.

`config_id` identifies a benchmark configuration, `operational_id` is used in
file paths and table joins, and `manuscript_name` gives its display name.
DeepVelo and DeepVelo_scv are distinct configurations. DeepVelo_scv follows
the official DeepVelo Figure 2 AE reconstruction workflow; GraphVelo_scv and
GraphVelo_dyn use different upstream velocity inputs. The operational ID
`velovi` is displayed as **veloVI**.

Use [canonical inputs](../data/README.md) and the method's
[software environment](../docs/environments.md).

The input root contains dataset directories such as `Pancreas/Pancreas.h5ad`;
it points to the scientific data directory described in the
[data layout](../data/README.md#local-data-roots). The repository's own `data/`
directory contains registries and download information. Launchers default to
a scientific `data/` directory beside the checkout; select another location
with `--data-root`, or edit the launcher's `DATA_ROOT` where listed below.

Randomness handling follows the corresponding upstream workflows, with
configuration-specific settings listed in the registry's `randomness` column.

Run commands from the repository root. Launchers write native H5ADs under an
explicit output root or the ignored `outputs/` default. Most retain scientific
settings and dataset lists near the top of the script; CLI dataset selection
is available where described below. Evaluation selects its frozen result files
through separate input maps.

## Paths and output files

| Script | Where to select paths / datasets | Main output relative to the selected output root |
| --- | --- | --- |
| `scvelo.py` | `--data-root`, `--output-root`, repeated `--dataset` and `--mode` | `scvelo_<mode>/scv_<mode>_<dataset>.h5ad`; dynamical uses `scvelo_dynamic/scv_dynamic_<dataset>.h5ad` |
| `unitvelo.py` | `--data-root`, `--output-root`, repeated `--dataset`, `--fit-option` | `unitvelo/` or `unitvelo_mode2/`, then `unitvelo_<dataset>.h5ad` |
| `velovi.py` | `--data-root`, `--output-root`, repeated `--dataset`, `--evaluation-role` | `<dataset>/veloVI_<dataset>_all.h5ad` |
| `deepvelo.py` | Top-level `DATA_ROOT`; `--output-root`, repeated `--dataset`, `--evaluation-role` | `<dataset>/deepvelo_<dataset>.h5ad`, plus `run_config.json` |
| `tfvelo.py`, `tsvelo.py` | `--source-root`, `--data-root`, `--output-root`, `--dataset` | `TFvelo_<official_name>_demo/TFvelo.h5ad` or `<official_name>/TSvelo.h5ad` |
| `cell2fate.py` | Top-level `DATA_ROOT`, `RESULTS_ROOT`, `OUTPUT_ROOT`, `DATASETS_TO_RUN` | `cell2fate/<dataset>/cell2fate_<dataset>.h5ad` |
| `dynamo.py`, `multivelo.py`, `tivelo.py`, `veloae.py` | Top-level `DATA_ROOT`, `OUTPUT_ROOT`, `DATASETS_TO_RUN` | Method subdirectory, then `<prefix><dataset>.h5ad`; prefixes are `dynamo_`, `multivelo_`, `tivelo_` and `veloAE_` |
| `celldancer.py`, `deepkinet.py` | Top-level `DATA_ROOT`, `OUTPUT_ROOT`, `DATASETS_TO_RUN` | `celldancer_<dataset>.h5ad` or `DeepKINET_<dataset>.h5ad` |
| `deepvelo2022.py`, `sdevelo.py` | Top-level `DATA_ROOT`, `OUTPUT_ROOT`, `DATASETS_TO_RUN` | Dataset subdirectory with the configured output stem; C18 uses `deepvelo_<stem>_recon.h5ad` |
| `graphvelo_scv.py`, `graphvelo_dynamo.py` | Top-level `DATA_ROOT`, `OUTPUT_ROOT`, `DATASETS_TO_RUN` | `<dataset>/graphvelo<dataset>.h5ad` or `<dataset>/graphvelo_dynamo_<dataset>.h5ad` |

In constant-configured scripts, `OUTPUT_ROOT` is the exact root used in the
source; its default may already include the method name. Dataset names in
commands follow each launcher's `DATASETS` keys. The dataset registry gives
their paper names; use the explicit input map when moving to an evaluator.

## Modes and example commands

| Method | Selection |
| --- | --- |
| scVelo, C01–C03 | `--mode deterministic`, `stochastic` or `dynamical`; repeat for several modes; default `MODES_TO_RUN` |
| UniTVelo, C05–C06 | `--fit-option 1` or `2`; repeat for both; default `FIT_OPTIONS_TO_RUN` |
| TFvelo / TSvelo | `--dataset <name>` or `all`; `--check-only` checks the selected source/input layout |
| veloVI | `--evaluation-role cell` or `gene` |
| DeepVelo | `--evaluation-role cell`, `gene` or `kinetics` |

```bash
python methods/scvelo.py --mode dynamical --dataset Pancreas --data-root /path/to/scientific-data \
  --output-root /path/to/scvelo-run
python methods/unitvelo.py --fit-option 2 --dataset Pancreas --data-root /path/to/scientific-data \
  --output-root /path/to/unitvelo-run
python methods/velovi.py --evaluation-role gene --dataset Pancreas --data-root /path/to/scientific-data \
  --output-root /path/to/empty/velovi-gene-run
python methods/deepvelo.py --evaluation-role cell \
  --output-root /path/to/empty/deepvelo-cell-run
```

scVelo and UniTVelo accept repeated `--dataset` selections, plus `--data-root`
and `--output-root`; modes use separate directories and reject existing H5ADs.
scVelo dynamical results use `scvelo_dynamic/`.
veloVI also accepts `--data-root`; DeepVelo uses its top-level `DATA_ROOT`.
Use separate empty roots for distinct evaluation roles.

### Hardware selection

Set the device for your machine before running. TSvelo accepts `--cuda`
(default 2); specify the intended visible-device index explicitly. UniTVelo
uses top-level `GPU_INDEX` (default 0); veloAE uses `GPU_INDEX` (default 3).
These indices select devices, rather than specify a minimum GPU count.
DeepVelo, DeepVelo_scv and SDEvelo inherit `CUDA_VISIBLE_DEVICES`.
CPU parallelism is set by `--n-jobs` where available, or constants such as
`N_JOBS`, `GRAPH_N_JOBS` and `CELL_DANCER_N_JOBS` in the selected script.

### Source-backed methods

Use the fixed upstream commits and patches below. TFvelo and TSvelo need the
benchmark's data-entry and relative-path adaptations; DeepKINET needs its
package-import and checkpoint-path fixes. Model code and training settings
follow the retained workflows.

[source_manifest.json](source_manifest.json) records the commits, patches and
SHA-256 hashes of the resulting runtime files and required TF databases.
The pinned commits and patches reproduce the retained source content;
historical checkout fingerprints were not recorded.

Run these commands from the benchmark repository root, using your own source
directories and the selected method's environment. Install the dependencies
described by the pinned upstream repository before launching TFvelo/TSvelo.

**TFvelo (C10)**

```bash
BENCHMARK_ROOT="$PWD"
TFVELO_ROOT="/path/to/TFvelo"
git clone https://github.com/xiaoyeye/TFvelo.git "$TFVELO_ROOT"
git -C "$TFVELO_ROOT" checkout --detach ec6cdb940af94f02fe32c8fcdb98494cdb4beb96
git -C "$TFVELO_ROOT" apply "$BENCHMARK_ROOT/methods/source_patches/tfvelo.patch"
unzip -q "$TFVELO_ROOT/ENCODE.zip" -d "$TFVELO_ROOT"
python methods/tfvelo.py --source-root "$TFVELO_ROOT" --dataset Pancreas \
  --data-root /path/to/scientific-data --output-root /path/to/empty/tfvelo-run
```

The patched source keeps the upstream `args.data_path` analysis input. The
launcher supplies each dataset's input and labels in a temporary relative
layout, then copies its native result to the selected output root.

**TSvelo (C14)**

```bash
BENCHMARK_ROOT="$PWD"
TSVELO_ROOT="/path/to/TSvelo"
git clone https://github.com/lijc0804/TSvelo.git "$TSVELO_ROOT"
git -C "$TSVELO_ROOT" checkout --detach 6e0ab8264e02aa511576278219fd7c1c65efec80
git -C "$TSVELO_ROOT" apply "$BENCHMARK_ROOT/methods/source_patches/tsvelo.patch"
unzip -q "$TSVELO_ROOT/ENCODE.zip" -d "$TSVELO_ROOT"
python methods/tsvelo.py --source-root "$TSVELO_ROOT" --dataset Pancreas --cuda 0 \
  --data-root /path/to/scientific-data --output-root /path/to/empty/tsvelo-run
```

The source patch adds the benchmark input branch and the relative resource
and output paths used by the launcher. `--source-root` may also be supplied
through the exported `TFVELO_ROOT` or `TSVELO_ROOT` variable. For either method,
add `--check-only` to check required source/input paths without training.

**DeepKINET (C13)**

```bash
BENCHMARK_ROOT="$PWD"
DEEPKINET_ROOT="/path/to/DeepKINET"
git clone https://github.com/3254c/DeepKINET.git "$DEEPKINET_ROOT"
git -C "$DEEPKINET_ROOT" checkout --detach fd9cb5de02e9139322fdef12d1629cf9fc2929e5
git -C "$DEEPKINET_ROOT" apply "$BENCHMARK_ROOT/methods/source_patches/deepkinet.patch"
python -m pip install "$DEEPKINET_ROOT"
```

The patch fixes relative imports in three package files and allows the
default checkpoint location to be selected with `DEEPKINET_CHECKPOINT`.
The benchmark launcher passes a separate checkpoint path for every dataset.
Set its `DATA_ROOT`, `OUTPUT_ROOT` and `DATASETS_TO_RUN`, then run
`python methods/deepkinet.py`. The launcher requires scVelo 0.2.5; retained
package-version evidence is described in the [environment guide](../docs/environments.md).

**Check reconstructed source content**

After applying a patch, compare the source files against the manifest.
For TFvelo, run the following from the benchmark root; substitute
`"$TSVELO_ROOT" C14` or `"$DEEPKINET_ROOT" C13` for the other sources.

```bash
python - "$TFVELO_ROOT" C10 <<'PY'
from pathlib import Path
import hashlib
import json
import sys

record = json.loads(Path("methods/source_manifest.json").read_text())["methods"][sys.argv[2]]
for relative, expected in record["files"].items():
    path = Path(sys.argv[1]) / relative
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise SystemExit(f"Source differs: {relative}")
print(f"Verified {len(record['files'])} source/resource files")
PY
```

## Method-specific inputs and workflows

- **cell2fate and DeepKINET:** raw-count inputs, integer-count checks, cell/gene
  alignment and protection against repeated preprocessing follow the launcher.
- **DeepKINET (C13):** the benchmark uses seed 123 with deterministic
  PyTorch/CUDA settings. Source identity is recorded in the registry.
- **TFvelo (C10):** provide the patched checkout above with `--source-root` or
  `TFVELO_ROOT`.
- **TSvelo (C14):** provide the patched checkout above with `--source-root` or
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

## Use the results

The launcher's output path and an evaluator's frozen input locator may have
different directory layouts. For a new run, use a separate evaluation input
map/configuration with the new relative path, velocity layer and actual file
identity. Cell accepts `--inputs`; Gene uses `evaluation/gene/config.py`;
kinetics accepts `--config`. Arrange the file under the
[root expected by that evaluator](../data/README.md#local-data-roots).

Successful method execution writes the H5AD shown in the table above.
scVelo, UniTVelo and veloVI additionally store `uns["benchmark_run"]`;
DeepVelo writes `run_config.json`. The evaluation guides identify the layers,
cell/gene alignment and result tables needed to check the next stage.

[Cell](../evaluation/cell/README.md), [Gene](../evaluation/gene/README.md) and
[kinetics](../evaluation/kinetics_recovery/README.md) document the exact result
family, layer and alignment consumed by each task. Ensembles are
[post-hoc comparisons](../ensemble/README.md) outside the method registry.
See the [reproduction guide](../docs/reproducibility.md) for input identities,
software environments and evaluation commands.
