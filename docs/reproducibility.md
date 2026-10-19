# Reproduction routes

Choose a route according to the result you need. Commands and exact numerical
settings live in the linked module guides.

| Route | Inputs | Guide |
| --- | --- | --- |
| Inspect scores or render CSV-only summaries | Committed compact tables and plotting dependencies | [Evaluation index](../evaluation/README.md) |
| Recalculate metrics | Frozen method H5ADs, references and evaluation environments | [Cell](../evaluation/cell/README.md), [Gene](../evaluation/gene/README.md), [kinetics](../evaluation/kinetics_recovery/README.md) |
| Retrain a method | Canonical inputs, registered method configuration and upstream software | [Methods](../methods/README.md) |
| Generate new simulations | Simulator and scenario configuration | [Simulation](../evaluation/gene/simulation/README.md) |
| Assemble a paper figure | Plot-data, source panels and any editable author source | [Figure map](figure_map.md) |

## Inputs and outputs

Start with the [data guide](../data/README.md). It distinguishes direct input
downloads, upstream-only references and packaged method outputs.
Map extracted archive members to the logical input paths expected by the
selected evaluator.

The [local-root table](../data/README.md#local-data-roots) defines those
layouts. Cell takes the parent of `velocity/`; Gene and ensemble take the
velocity directory; kinetics takes a root containing its reference and
method-result aliases. Python interpreter paths and environment names belong
to your installation; [environment roles](environments.md#evaluation)
identify the required software.

A method produces a native result; evaluation consumes its recorded layers,
cell/gene IDs and embeddings, then writes compact tables and panels.
Keep fields consumed by evaluation when retaining an H5AD. A final figure may
also require manual assembly.

For a new method output, point a separate evaluator input map/configuration
at its path and actual file identity. The committed input maps select the
frozen benchmark outputs. Plotting and summary commands state which tables
they read; selecting a new output directory alone changes where files are
written. Use the owning module's documented input options to select other
tables.

| Information | Source |
| --- | --- |
| Method identity, version, options and output keys | [methods/configurations.tsv](../methods/configurations.tsv) |
| Dataset identity and upstream origin | [data/datasets.tsv](../data/datasets.tsv) |
| Artifact identity, checksum and access | [artifacts.tsv](../artifacts.tsv) and [data guide](../data/README.md) |
| Files, layers and alignment used in a comparison | The evaluation module's input map and configuration |
| Numerical processing and aggregation | The producing module's README and code |

## Randomness and saved results

Method-specific randomness settings are listed in
[methods/configurations.tsv](../methods/configurations.tsv). To recalculate the
reported scores, use the saved method outputs and reference assets identified
by the module input maps and manifests. Simulation generation produces new
realizations; the registered historical H5ADs identify the benchmark inputs.

Kinetics uses the `central_results` inputs, the corrected organoid reference
direction and the labelled-data Dynamo control described in its module guide.
The [historical comparison](../evaluation/kinetics_recovery/plotting/README.md#historical-comparison)
documents cross-runtime numerical differences.

Use the [environment guide](environments.md) and fresh output directories.
Run relevant checks before promoting a new
result or figure; the [figure map](figure_map.md) distinguishes numerical
rendering from retained author exports.

## Historical evidence

The [fixed-neighbor reproduction record](fixed_neighbor_reproduction_20260912.md)
is a checksum-registered historical record.
The earlier simulation audit is retained at
[revision 8b4db85](https://github.com/Tomorin259/benchmark_repository/blob/8b4db85127b3333e6d60960d1e73359f5a5a81b6/docs/simulation_reproduction_audit_20260914.md).
Use the current module guides for runnable workflows.
