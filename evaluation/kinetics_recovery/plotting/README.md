# Kinetics plotting

These renderers consume frozen central-result H5ADs and their common evaluator
output. Scientific inputs, processing, commands, and scope are documented in
[the module README](../README.md).

| Entry | Input/output contract |
| --- | --- |
| `figure5.py` / `figure5_config.json` | Evaluator output root; renders B--D/F and exports exact plotted rows |
| `streams.py` / `stream_config.json` | Velocity root; validates and renders saved X/V for four RPE1 and three organoid methods |
| `compare_figure5.py` | Explicit archived/current canonical directories and a fresh comparison output directory |

B applies display-only z-scores. C uses selected-gene Spearman scores. D
classifies reference and predictions independently using cosine thresholds
-0.5/+0.5 and exports class counts and membership. F filters unified scores by
method, dataset, parameter, and branch, then displays selected-gene means.
DeepKINET supplies gamma only. TSvelo alpha uses two lineage files; its stream
uses the merged result. Dynamo A has no automatic renderer.

Streams validate cell IDs, finite arrays, dimensions, and any configured
reference embedding. Exact equality after casting reference coordinates to
the saved dtype is permitted for storage precision differences. Saved X/V
are explicitly supplied to scVelo; no neighbors, moments, graph, or projection
recomputation occurs. Missing independent provenance is disclosed in manifests.

RPE1 export uses a tight bounding box with 0.1-inch padding at 600 dpi, matching
the retained Dynamo canvas. Organoid export is unchanged. Render into fresh
directories; current assets under `../figures/fig5/` must not be overwritten
by an unreviewed reproduction run. Legacy-source validation additionally
requires explicit `--legacy-config`, `--legacy-artifact-root`, and
`--metrics-root`; it must not infer archived inputs from the current Git HEAD.
