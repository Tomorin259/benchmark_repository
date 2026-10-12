# Kinetics plotting

These renderers consume frozen central-result H5ADs and their common evaluator
output. Scientific inputs, processing, commands, and scope are documented in
[the module README](../README.md).

| Entry | Input/output contract |
| --- | --- |
| `figure6.py` / `figure6_config.json` | Evaluator output root; renders B--D/F and exports exact plotted rows |
| `streams.py` / `stream_config.json` | Velocity root; validates and renders saved X/V for four RPE1 and four conventional organoid methods |
| `compare_figure6.py` | Explicit archived/current canonical directories and a fresh comparison output directory |
| `expanded_reference.py` | Isolated expanded scores; reuses Figure 6B/6F functions for three 600-dpi boxplot PNGs and one 10-page RPE1 PDF |

Expanded-reference materials retain the native 18-by-8-inch Figure 6F canvas,
without tight cropping or resizing. Labels give finite gene n; smaller gray
points reduce crowding. Method order is fixed, not an expanded-score ranking.
The shared boxplot function's default arguments still reproduce the selected
Figure 6F display. The trajectory function reserves a separate legend slot
at multiples of ten; expanded pages request a complete method legend rather
than inheriting coverage from their first gene. Existing 19-gene rendering
uses its original legend behavior. See [the expanded results contract](../results/expanded_reference/README.md).

B applies display-only z-scores. C uses selected-gene Spearman scores. D
classifies reference and predictions independently using cosine thresholds
-0.5/+0.5 and exports class counts and membership. F filters unified scores by
method, dataset, parameter, and branch, then displays selected-gene means.
DeepKINET supplies gamma only. TSvelo alpha uses two lineage files; its stream
uses the merged result. Dynamo A and E are recorded labelled-data manual assets, with provenance stored beside the panels; neither has an automatic H5AD renderer in this module.

The B/C/F renderers label the control `Dynamo (scEU)` without renaming
evaluator IDs or exported plot-data. The C heatmap shows
`Mean (n)`, where the mean uses available displayed-gene scores
and n counts valid genes. Rows remain sorted by these unrounded means, not by
full-reference scores or a common cross-method gene intersection. Summary cells
share the gene heatmap's row boundaries and color scale; their color encodes
the unrounded mean, not n. Their numeric labels use the same regular 10-point
font as the gene cells. These labels
do not rename evaluator IDs or alter scores. The shared colorbar is narrow;
missing cells retain the white `Not Calculated` encoding.

D uses aligned turnover-category headings and shows Macro-F1 below the three
Venn diagrams: 0.324 for cellDancer and 0.340 for DeepVelo. Cosine intervals
and the three-class Macro-F1 definition are documented in the module README;
balanced accuracy is not displayed. Venn counts and classification thresholds
are unchanged.

Individual scores and summary annotations (means, standard deviations, medians,
macro-averages, branch averages, and differences between summary values) use
four decimal places, except the C heatmap's `Mean (n)` column (two) and
D's Macro-F1 annotation (three).
Dense heatmap annotations also use two;
formal numerical tables use four.
Positive scores have no leading plus sign. Raw CSV files retain full
precision, and calculations and ranking use unrounded values.

Figure 6F organoid boxplots use the benchmark method colors in
`evaluation/cell/plots.py` (`BENCHMARK_COLORS`). The labelled Dynamo control
uses the Dynamo(RNA_seq) color. Figure 6B keeps a separate, fixed line palette
because several benchmark colors are difficult to distinguish when the
trajectories overlap. Its source values and method order are unchanged by this
display exception. Figure 6A/E streams color cells by cluster, not methods.

Streams validate cell IDs, finite arrays, dimensions, and any configured
reference embedding. Exact equality after casting reference coordinates to
the saved dtype is permitted for storage precision differences. Saved X/V
are explicitly supplied to scVelo; no neighbors, moments, graph, or projection
recomputation occurs. Missing independent provenance is disclosed in manifests.
The organoid TSvelo stream uses its registered merged result and the same
reference UMAP and cluster palette as the other E panels. This adds an E
velocity-field display only, not a TSvelo parameter score in F.

RPE1 export uses a tight bounding box with 0.1-inch padding at 600 dpi, matching
the retained Dynamo canvas. Organoid export is unchanged. Render into fresh
directories; current assets under `../figures/fig6/` must not be overwritten
by an unreviewed reproduction run. Legacy-source validation additionally
requires explicit `--legacy-config`, `--legacy-artifact-root`, and
`--metrics-root`; it must not infer archived inputs from the current Git HEAD.
