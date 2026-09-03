# Method reproduction

This directory contains one direct reproduction script per method or shared
official implementation. The 20 benchmark configurations remain explicit in
[`configurations.tsv`](configurations.tsv); configurations that share an
implementation are command-line modes of one script.

Each script accepts explicit input and output locations, keeps the upstream
tutorial sequence visible, validates the result keys consumed by evaluation,
and writes a method-native H5AD result. Upstream repositories, checkpoints,
Conda installations, logs, and generated H5AD files are not vendored here.

Methods that require raw counts own that input logic locally. In particular,
DeepKINET and cell2fate document the `capo_pre`, `mousebrain`, `neuron`, and
`rpe1_gen` `_counts.h5ad` selections, cell/gene alignment, integer-count
checks, and the rule preventing an accidental second preprocessing pass. No
count-specific branch is added to methods that do not consume counts.

The script must not depend on a universal runner, adapter registry, plugin
system, or editable checkout of an upstream repository. A helper is justified
only by a real scientific operation, input-format boundary, or invariant.
