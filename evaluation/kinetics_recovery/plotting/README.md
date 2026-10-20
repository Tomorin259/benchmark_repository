# Kinetics plotting

Use [the kinetics guide](../README.md) for inputs, processing and commands.

| Entry | Input and output |
| --- | --- |
| `figure6.py` / `figure6_config.json` | Evaluator output; renders B–D/F and exact plot-data |
| `streams.py` / `stream_config.json` | Saved velocity-root H5ADs; four RPE1 and four conventional organoid methods |
| `expanded_reference.py` | Expanded scores; three boxplots and a ten-page RPE1 trajectory PDF |

## Display settings

B independently z-scores processed reference/prediction curves.
C displays selected-gene Spearman scores and `Mean (n)`, sorted by unrounded
displayed-gene means. Its summary cells use the same heatmap scale and row
boundaries; color encodes mean, not count. Missing cells remain white
(`Not Calculated`).
D shows the three turnover-category Venn diagrams and Macro-F1.
F uses available selected-gene means by method, parameter and branch.
DeepKINET supplies gamma only; TSvelo rates use lineage files, while its
velocity field uses the merged file. The labelled Dynamo RPE1 panel has a
[Source Data renderer and inputs](../../../data/SOURCE_DATA.md); the organoid
panel is retained as an author-assembled image.

Individual/summary annotations and formal tables use four decimals, dense
heatmaps and C's `Mean (n)` use two, and D's Macro-F1 uses three.
CSV precision and sorting are unchanged. F inherits
`BENCHMARK_COLORS` from `evaluation/cell/plots.py`;
B retains a separate high-contrast line palette.
Stream colors encode clusters. RPE1 streams use a tight bounding box,
0.1-inch padding and 600 dpi.

Expanded boxplots retain the 18-by-8-inch F canvas without resizing/cropping,
finite-gene counts, smaller gray points and fixed method order.
Expanded trajectories reserve a separate legend slot at multiples of ten
and request a complete method legend. The selected-gene renderer retains
its original legend behavior.
[Expanded scope and outputs](../results/expanded_reference/README.md) remain
separate from selected Figure 6 scores.

Streams reuse saved X/V after alignment, dimension, finite-value and configured
reference checks. Casting the reference to the saved dtype may account for
storage precision; rendering does not recompute neighbors, moments or graphs.

## Historical comparison

<details>
<summary>Earlier kinetics sources and comparison commands</summary>

The retained `../archive/fig5/` contains unique historical sources.
These use the former kinetics Figure 5 numbering; current Figure 5 is
gene-level evaluation. Current Figure 7 hashes belong to its own composite
manifest, not the historical Figure 7 hash retained inside the kinetics
promotion record.

Materialize the full original bundle, then use `compare_figure6.py` with
explicit historical/current directories and a fresh comparison directory:

```bash
mkdir -p /path/to/historical-kinetics
git archive bed37b629dbc966109c8cb92dd22aca0810caf2f evaluation/kinetics_recovery/archive/fig5 \
  | tar -x -C /path/to/historical-kinetics
python evaluation/kinetics_recovery/plotting/compare_figure6.py \
  --old /path/to/historical-kinetics/evaluation/kinetics_recovery/archive/fig5/panels \
  --new evaluation/kinetics_recovery/figures/fig6 \
  --output /path/to/fresh-comparison
```

Archived-source rendering requires explicit `--legacy-config`,
`--legacy-artifact-root` and `--metrics-root`.
Use the source identities in the manifests.

</details>
