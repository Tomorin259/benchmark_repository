# Gene-level velocity evaluation

Nine method configurations are evaluated on seven real datasets and the three
fixed simulation scenarios. The canonical CBDir comparison uses the fixed
reference-30 graph. Native H5ADs are retained; numerical scoring reverses
veloVI's additional independent min-max axes for both state and velocity.

## Results and selection

| Source | Use |
| --- | --- |
| [results/gene_cbdir_reference30_all.csv](results/gene_cbdir_reference30_all.csv) | Canonical `new_cbdir` scores |
| [results/gene_cbdir_edge_details_reference30_all.csv](results/gene_cbdir_edge_details_reference30_all.csv) | Transition-level results |
| [results/simulation/gene_cbdir_cossim_reference.csv](results/simulation/gene_cbdir_cossim_reference.csv) | Simulation CosSim with restored veloVI units; its historical CBDir column is a different source |
| [results/unit_restoration/](results/unit_restoration/provenance.json) | Own-source scale factors, native/restored score pairs, input hashes and exact-forward checks |
| [results/common_genes/](results/common_genes/) | Finite nine-method common-gene comparison on seven real datasets |
| [selection/](selection/README.md) | Final genes, screening rules and candidate-table evidence |
| [config.py](config.py) | Input paths, layer keys, genes and transition edges |
| [figures/panels/assembly_manifest.tsv](figures/panels/assembly_manifest.tsv) | Paper panel paths, score sources and checksums |

Real datasets are `retina`, `DentateGyrus_lamanno`, `gastrulation`,
`capo_pre`, `bonemarrow`, `mousebrain` and `pancreas`;
simulation IDs are `branching`, `murk` and `all`.
[results/README.md](results/README.md) describes the dataset/scenario layout.

## Recalculate CBDir

Run from the repository root with the external artifacts arranged at their
recorded logical paths:

```bash
python evaluation/gene/run_reference30_gene_cbdir.py \
  --artifact-root /path/to/velocity-artifacts \
  --neighbor-mode reference \
  --datasets retina DentateGyrus_lamanno gastrulation capo_pre bonemarrow mousebrain pancreas branching murk all \
  --output-root outputs/gene-evaluation
```

The explicit list covers all ten datasets/scenarios; the script's default list
is a subset. `--methods veloVI` isolates the affected configuration; the default
still evaluates all nine. The output directory must be empty.
`VELOCITY_ARTIFACT_ROOT` can supply the artifact root.
The evaluator also reads the historical comparison-score files named in its
source/provenance; the saved H5ADs alone do not supply those comparisons.

`--neighbor-mode native` is a controlled comparison, with a separate output
root. Results include per-dataset scores, edge details, neighbor assets and
audit/provenance records. Every dataset reuses its committed reference-neighbor
NPZ; reevaluation does not rebuild neighbors.

For veloVI, the evaluator reconstructs its own pre-min-max moments from saved
abundances and the binary, self-connected, row-normalized graph using float32
moment arithmetic. `capo_pre` instead uses the actual preprocessed input
`data/capo_pre/capo_pre.h5ad`, because the saved final graph is not its moment
origin; identifiers and abundances must agree exactly. Forward min-max scaling
must reproduce saved `Ms`/`Mu` bit for bit. For each axis, `y = a*x + b` implies
`v_x = v_y/a`; the translation does not enter velocity. A failed exact-forward
check stops evaluation. Other methods retain their existing numerical units.

Simulation runs also write `simulation/velovi_cossim.csv`. This compares raw
restored veloVI velocities with unchanged `velocity_s_true`/`velocity_u_true`,
using the historical norm mask (`>1e-6` for prediction and truth) without
CBDir's velocity cleanup. Join this output by dataset/method/gene to replace
only the 70 veloVI CosSim fields in the compact 630-row table. Its other 560
CosSim values and all historical CBDir fields remain unchanged.

Check the common-gene tables without H5AD inputs:

```bash
python evaluation/gene/refresh_common_gene_summary.py
```

Add `--write --output-root outputs/common-genes` to regenerate them separately.
The finite intersections contain 6/6/6/7/2/4/4 genes in the real-dataset order
listed above; macro means give each dataset equal weight.

The benchmark uses the registered simulation H5ADs. Use the
[simulation generator](simulation/README.md) to create additional datasets.

## Render summaries and phase portraits

For CSV-only summaries, install the [plotting environment](../../docs/environments.md):

```bash
python evaluation/gene/plot_gene_assets.py --summaries-only \
  --output-root outputs/gene-summary
```

Phase/comparison panels additionally need the registered H5AD inputs and the
AnnData/Scanpy/scVelo plotting dependencies:

```bash
python evaluation/gene/plot_gene_assets.py \
  --artifact-root /path/to/velocity-artifacts \
  --output-root outputs/gene-panels
```

Use `--cases bonemarrow --genes LMNA` for focused phase panels.
Summaries still cover every registered gene in the selected dataset.
`--dpi` defaults to 600; output directories must be empty and cannot be the
committed panel directory. The renderer joins `new_cbdir` and CosSim by
dataset/method/gene, with registry method names and the Cell color palette.

### Display and source boundaries

[plotting/phase_portrait.py](plotting/phase_portrait.py) keeps linear
`Ms`/`Mu` coordinates and uses deterministic representative-cell sampling,
without neighborhood velocity smoothing:

| Setting | Value |
| --- | --- |
| Arrow length transform | `log1p` |
| Reference quantile / length | 0.90 / 0.05 |
| Maximum length | 0.10 |
| Grid / density percentile | 25 × 25 / 15 |

Real plotting copies replace NaN/Inf velocities by zero and clip to ±1e10;
simulation vectors keep their native values. Reference H5ADs supply cluster
labels/colors, including the retained equal-length label-assignment fallback.
Rendering uses copies of the source arrays. Summary jitter affects display only.

Exports include `scores_used.csv`, `provenance.json` and
`assembly_manifest.tsv`. They recreate scripted comparison panels and
summaries. The panel manifest records current hashes; panel provenance retains
the original generation identity. Final composites follow the
[assembly steps in the figure map](../../docs/figure_map.md).
Figure 4D's final paired Cell/Gene panel belongs to
[Simulation Cell](../cell/simulation/README.md).

Individual scores, summary annotations and formal tables use four decimals;
dense heatmaps two. CSVs and sorting retain full precision.
The sampling design is adapted from cellDancer 1.1.7; see
[THIRD_PARTY_NOTICES.md](plotting/THIRD_PARTY_NOTICES.md).

### Selected-edge and cross-dataset source panels

`plot_selected_edge_panels.py --artifact-root /path/to/velocity-artifacts
--output-root outputs/selected-edge-panels` redraws the three Figure 3 cases
from the native phase arrays, current transition-level Gene scores and
`results/selected_edge_cell_controls.csv`. The retained Cell controls refer
to the same annotated transitions. Each export includes plotted score rows.

`plot_real_comparisons.py --output-root outputs/real-gene-comparisons`
redraws the Figure 5 cross-dataset comparisons and the S9 common-gene panel.
The Figure 5 scatters weight dataset means equally; the adjacent branching
boxplot and heatmap retain their pooled-gene distribution. The S9 panel uses
the finite nine-method intersection in `results/common_genes/`.

Simulation phase colors retain the Figure 4 progression-bin gradient from
early to late. Phase arrays and velocities are plotting copies; this palette
and the canvas-bound checks do not change numerical evaluation.
