# Method reproduction

This directory contains one direct reproduction script per method or shared
official implementation. The 20-method benchmark configurations remain
explicit in configurations.tsv, the sole 20-method registry. It has exactly
20 unique config_id values and 20 unique operational_id values; ensemble
variants are not entries. Configurations that share an implementation are
listed in visible top-level mode constants in one script.

Each script begins with visible input/output paths and the dataset/configuration
list to edit for a focused run. It keeps the upstream tutorial sequence
visible, validates the result keys consumed by evaluation, and writes a
method-native H5AD result. It does not use a universal runner or a shared
method adapter. C10 and C14 have method-specific source-backed launchers
because their official demos require relative-path compatibility trees and
per-dataset input overrides. Upstream repositories, checkpoints, Conda
installations, logs, and generated H5AD files are not vendored here.

Methods that require raw counts own that input logic locally. In particular,
DeepKINET and cell2fate document the `capo_pre`, `mousebrain`, `neuron`, and
`rpe1_gen` `_counts.h5ad` selections, cell/gene alignment, integer-count
checks, and the rule preventing an accidental second preprocessing pass. No
count-specific branch is added to methods that do not consume counts.

The formal DeepKINET result uses seed 123 with deterministic PyTorch/CUDA
settings. This is a benchmark reproducibility setting; it is not presented as
an additional default from the upstream tutorial.

TFvelo and TSvelo retain the official no-global-seed behavior. Their
method-specific launchers keep the 14-dataset input map at the top of each
script, stage the official relative paths, and write generated outputs under
the ignored outputs/ tree. Set TFVELO_ROOT or TSVELO_ROOT (or pass
--source-root) before running them.

The script must not depend on a universal runner, adapter registry, plugin
system, or editable checkout of an upstream repository. A helper is justified
only by a real scientific operation, input-format boundary, or invariant.

Public names use veloVI spelling. DeepVelo's C16 graph-convolutional
configuration and C18 reconstruction configuration remain distinct, as do
the C19 scVelo-input and C20 Dynamo-input GraphVelo configurations. Their
scientific settings are kept in configurations.tsv.

Mean-4 and CenterWeighted-4 are independent comparators beside the
20-method benchmark, not methods in this directory. Each is ranked against
the 20 methods separately and the two comparators are not compared with one
another. Their implementations are in ensemble/, which shares evaluation/cell/
and evaluation/cell/reference_neighbors/; formal fixed-neighbor all-14
ensemble results are available under ensemble/results.

DeepVelo (C16) requires an explicit paper analysis role and an empty output
directory.  Use `--evaluation-role cell` for the cell-level default run; this
sets `pred_unspliced=False`.  The gene and kinetics roles set it to `True` and
remain separate from the cell result family.  For example:

```bash
python methods/deepvelo.py \
  --evaluation-role cell \
  --output-root /path/to/results/deepvelo/results_default_cell
```

To reproduce the historical launcher’s RNG scope, add
`--reset-seed-per-dataset`; this resets seed 123 immediately before each
dataset and records the choice in `run_config.json`.
The canonical historical cell-level input used by the published cell metrics
is the existing central no-vu family `results/deepvelo/deepvelo_*.h5ad`;
newly generated role-specific outputs should be kept in a separate audit
directory until promoted.
