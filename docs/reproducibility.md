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
Materialize files at the logical paths expected by the selected evaluator;
archive grouping paths are not automatically valid evaluator roots.

A method produces a native result; evaluation consumes its recorded layers,
cell/gene IDs and embeddings, then writes compact tables and panels.
Keep fields consumed by evaluation when retaining an H5AD. A final figure may
also require manual assembly.

| Information | Source |
| --- | --- |
| Method identity, version, options and output keys | [methods/configurations.tsv](../methods/configurations.tsv) |
| Dataset identity and upstream origin | [data/datasets.tsv](../data/datasets.tsv) |
| Artifact identity, checksum and access | [artifacts.tsv](../artifacts.tsv) and [data guide](../data/README.md) |
| Files, layers and alignment used in a comparison | The evaluation module's input map and configuration |
| Numerical processing and aggregation | The producing module's README and code |

## Randomness and reproducibility limits

Use the [environment guide](environments.md) and record the configuration
used for a rerun. There is no common random seed across the benchmark:
DeepKINET uses its documented seed-123 deterministic configuration, while
new simulation generation intentionally leaves the stream unseeded.

Reusing a fixed canonical input and generating a new realization are different
operations. Structural agreement or matching rounded scores does not establish
byte identity. Kinetics uses the frozen `central_results` inputs, the corrected
organoid reference direction and the labelled-data Dynamo exceptions described
in its module guide. Small cross-runtime full-table differences in DeepKINET
have not been assigned a cause.

Use fresh output directories. Run relevant checks before promoting a new
result or figure; the [figure map](figure_map.md) distinguishes numerical
rendering from retained author exports.

## Historical evidence

The [fixed-neighbor reproduction record](fixed_neighbor_reproduction_20260912.md)
is a checksum-registered historical record.
The earlier simulation audit is retained at
[revision 8b4db85](https://github.com/Tomorin259/benchmark_repository/blob/8b4db85127b3333e6d60960d1e73359f5a5a81b6/docs/simulation_reproduction_audit_20260914.md).
Use the current module guides for runnable workflows.
