# Kinetics plotting

These renderers consume frozen central-result H5ADs and their common evaluator
output. Scientific inputs, processing, commands, and scope are documented in
[the module README](../README.md).

| Entry | Input/output contract |
| --- | --- |
| `figure6.py` / `figure6_config.json` | Evaluator output root; renders B--D/F and exports exact plotted rows |
| `streams.py` / `stream_config.json` | Velocity root; validates and renders saved X/V for four RPE1 and three conventional organoid methods |
| `compare_figure6.py` | Explicit archived/current canonical directories and a fresh comparison output directory |

B applies display-only z-scores. C uses selected-gene Spearman scores. D
classifies reference and predictions independently using cosine thresholds
-0.5/+0.5 and exports class counts and membership. F filters unified scores by
method, dataset, parameter, and branch, then displays selected-gene means.
DeepKINET supplies gamma only. TSvelo alpha uses two lineage files; its stream
uses the merged result. Dynamo A and E are recorded labelled-data manual assets, with provenance stored beside the panels; neither has an automatic H5AD renderer in this module.

Figure 6 uses a fixed method palette: DeepVelo `#1f77b4`, cellDancer
`#ff7f0e`, Dynamo(label data) `#2ca02c`, DeepKINET `#B5BD61`, and TSvelo
`#e420ca`. The renderer applies this registry to both trajectory lines and
organoid boxplots; per-branch configuration colors must match the registry.

Streams validate cell IDs, finite arrays, dimensions, and any configured
reference embedding. Exact equality after casting reference coordinates to
the saved dtype is permitted for storage precision differences. Saved X/V
are explicitly supplied to scVelo; no neighbors, moments, graph, or projection
recomputation occurs. Missing independent provenance is disclosed in manifests.

RPE1 export uses a tight bounding box with 0.1-inch padding at 600 dpi, matching
the retained Dynamo canvas. Organoid export is unchanged. Render into fresh
directories; current assets under `../figures/fig6/` must not be overwritten
by an unreviewed reproduction run. Legacy-source validation additionally
requires explicit `--legacy-config`, `--legacy-artifact-root`, and
`--metrics-root`; it must not infer archived inputs from the current Git HEAD.
