# Run a method

[configurations.tsv](configurations.tsv) lists the 20 configurations, versions,
source references, scientific options and required result keys.
`release_status=public_launcher` records code access; version and source fields
retain the available historical identity.

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
| `deepvelo.py` | `--data-root`, `--output-root`, repeated `--dataset`, `--evaluation-role` | `<dataset>/deepvelo_<dataset>.h5ad`, plus `run_config.json` |
| `tfvelo.py`, `tsvelo.py` | `--source-root`, `--data-root`, `--output-root`, `--dataset` | `TFvelo_<official_name>_demo/TFvelo.h5ad` or `<official_name>/TSvelo.h5ad` |
| `cell2fate.py` | Top-level `DATA_ROOT`, `RESULTS_ROOT`, `OUTPUT_ROOT`, `DATASETS_TO_RUN` | `cell2fate/<dataset>/cell2fate_<dataset>.h5ad` |
| `multivelo.py`, `tivelo.py` | `--data-root`, `--output-root`, repeated `--dataset` | `multivelo/multivelo_<dataset>.h5ad` or `tivelo/tivelo_<dataset>.h5ad` |
| `celldancer.py` | `--data-root`, `--output-root`, repeated `--dataset` | `celldancer_<dataset>.h5ad` |
| `deepvelo2022.py`, `graphvelo_scv.py` | `--data-root`, `--output-root`, repeated `--dataset` | `<dataset>/deepvelo_<stem>_recon.h5ad` or `<dataset>/graphvelo<dataset>.h5ad` |
| `dynamo.py`, `veloae.py`, `deepkinet.py`, `sdevelo.py`, `graphvelo_dynamo.py` | Top-level `DATA_ROOT`, `OUTPUT_ROOT`, `DATASETS_TO_RUN` | Method-native output path declared in the script |

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
veloVI and DeepVelo also accept `--data-root`.
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
The pinned revisions and patches reconstruct the retained runtime files.
File hashes and source-matching details are recorded in the manifest.

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

- **cell2fate and DeepKINET:** use raw spliced/unspliced counts aligned by
  cell and gene ID.
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

### UniTVelo Gene layers

Run C05 (`--fit-option 1`), then differentiate its saved fitted curve with
[unitvelo_velocity_u.py](unitvelo_velocity_u.py):

```bash
python methods/unitvelo_velocity_u.py \
  --input /path/to/unitvelo-run/unitvelo/unitvelo_Pancreas.h5ad \
  --output /path/to/unitvelo-gene/unitvelo_Pancreas_with_velocity_u_mode1.h5ad
```

This is the retained benchmark `get_velocity_unspliced.py` calculation.
It adds `velocity_u`, `velocity_s_reconstructed`, `fit_s_reconstructed` and
`fit_u_reconstructed`, using the saved Gaussian curve, `fit_t`, fitted rates
and `velocity_genes` mask. Mode-1 missing times use the first complete shared
cell-time column; near-zero beta gives undefined reconstructed values.
The existing `velocity`, `Ms` and `Mu` are retained. This step uses AnnData,
NumPy and SciPy and performs no model fitting. Select the new output in the
Gene input configuration.

### Run methods on simulations

The nine Gene methods accept `sim_only_murk`, `sim_only_branching` and
`sim_hybrid` through `--dataset`. Their inputs are the three `fixed_simulation`
rows in [input_file_manifest.tsv](../data/input_file_manifest.tsv), selected
by [simulation_inputs.py](simulation_inputs.py). Supply the **velocity root**
as `--data-root`: these paths start with `indicators/gene/simulation/`.
Real-data commands keep the scientific-input root described above.

Select the method's environment before each command. For example, run the
normal + branching scenario with separate output directories:

```bash
SIM_INPUT_ROOT=/path/to/velocity-artifacts
SIM_OUTPUT_ROOT=/path/to/new-simulation-runs
python methods/scvelo.py --dataset sim_only_branching --mode dynamical \
  --data-root "$SIM_INPUT_ROOT" --output-root "$SIM_OUTPUT_ROOT/scvelo"
python methods/unitvelo.py --dataset sim_only_branching --fit-option 1 \
  --data-root "$SIM_INPUT_ROOT" --output-root "$SIM_OUTPUT_ROOT/unitvelo"
python methods/velovi.py --dataset sim_only_branching --evaluation-role gene \
  --data-root "$SIM_INPUT_ROOT" --output-root "$SIM_OUTPUT_ROOT/velovi"
python methods/deepvelo.py --dataset sim_only_branching --evaluation-role gene \
  --data-root "$SIM_INPUT_ROOT" --output-root "$SIM_OUTPUT_ROOT/deepvelo"
python methods/deepvelo2022.py --dataset sim_only_branching \
  --data-root "$SIM_INPUT_ROOT" --output-root "$SIM_OUTPUT_ROOT/deepvelo2022"
python methods/graphvelo_scv.py --dataset sim_only_branching \
  --data-root "$SIM_INPUT_ROOT" --output-root "$SIM_OUTPUT_ROOT/graphvelo"
python methods/multivelo.py --dataset sim_only_branching \
  --data-root "$SIM_INPUT_ROOT" --output-root "$SIM_OUTPUT_ROOT/multivelo"
python methods/tivelo.py --dataset sim_only_branching \
  --data-root "$SIM_INPUT_ROOT" --output-root "$SIM_OUTPUT_ROOT/tivelo"
python methods/celldancer.py --dataset sim_only_branching \
  --data-root "$SIM_INPUT_ROOT" --output-root "$SIM_OUTPUT_ROOT/celldancer"
```

Repeat `--dataset` for the other two scenarios. Fourteen real datasets remain
the default; simulations require explicit selection. The output filenames use
the selected dataset ID. Gene evaluation's frozen aliases `murk`, `branching`
and `all` correspond to `sim_only_murk`, `sim_only_branching` and `sim_hybrid`;
point its input configuration at the new native files before reevaluation.
UniTVelo additionally needs the four-layer reconstruction above.

The simulation workflows preserve 20-count/2,000-gene preprocessing and
30-PC/30-neighbor moments. Historical settings include:

| Method | Simulation-specific setting |
| --- | --- |
| scVelo_dyn / GraphVelo_scv | All genes marked before dynamical recovery |
| DeepVelo_scv | All genes used for the AE reconstruction target |
| UniTVelo | Mode 1, `R2_ADJUST=True`, `IROOT=None`, `AGENES_R2=1` |
| TIvelo | `group_key=bin`, `filter_genes=False`, 100 epochs, 12 workers |
| DeepVelo | Unspliced head, seed 123 reset per scenario, `mask_zero=False` |
| veloVI | Gene training and posterior extraction described above |

[source_manifest.json](source_manifest.json), `simulation_workflows`, records
the retained source hashes and settings. The public entry points use explicit
input and output paths. Retraining stochastic models produces a new run;
the published native files remain the inputs for replaying reported scores.

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
