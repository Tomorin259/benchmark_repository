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
| [results/simulation/gene_cbdir_cossim_reference.csv](results/simulation/gene_cbdir_cossim_reference.csv) | Ground-truth cosine similarity on simulations |
| [results/unit_restoration/](results/unit_restoration/provenance.json) | veloVI coordinate scaling, source identities and paired scores |
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
recorded logical paths. Here `--artifact-root` is the velocity directory
itself: a configured `results/...h5ad` is read at
`<artifact-root>/results/...h5ad`. See the
[root layouts](../../data/README.md#local-data-roots) and
[evaluation dependencies](../../docs/environments.md#evaluation).
Real-data runs also need the `label_reference` file in `config.py`;
veloVI on `capo_pre` additionally reads `data/capo_pre/capo_pre.h5ad`.

```bash
python evaluation/gene/run_reference30_gene_cbdir.py \
  --artifact-root /path/to/velocity-artifacts \
  --neighbor-mode reference \
  --datasets retina DentateGyrus_lamanno gastrulation capo_pre bonemarrow mousebrain pancreas branching murk all \
  --output-root outputs/gene-evaluation
```

The explicit list covers all ten datasets/scenarios; the script's default list
is a subset. `--methods veloVI` evaluates veloVI alone; the default evaluates
all nine methods. The output directory must be empty.
`VELOCITY_ARTIFACT_ROOT` can supply the artifact root.
Historical comparison-score files named in the source/provenance supply
`old_cbdir`, old edge scores and their deltas. Missing historical files leave
those comparison fields undefined; current CBDir and veloVI simulation
CosSim are calculated from the H5ADs and fixed reference assets.

`--neighbor-mode native` is a controlled comparison, with a separate output
root. Results include per-dataset scores, edge details, neighbor assets and
provenance. Reference mode reuses the committed reference-neighbor NPZs.

The output root contains `gene_cbdir_reference30_all.csv`,
`gene_cbdir_edge_details_reference30_all.csv`, `neighbor_audit_all.csv` and
`provenance.json`, plus per-dataset tables. veloVI also writes
`unit_restoration/velovi_scale_factors.csv` and, for simulations,
`simulation/velovi_cossim.csv`. Compare current scores by dataset/method/gene,
including status and evaluable counts, against the corresponding formal rows.
The configured method paths and layer keys select the result family for a
new run; keep its output root separate from the committed results.

### veloVI coordinate units

[coordinate_units.py](coordinate_units.py) reconstructs veloVI's pre-min-max
moments from saved abundances and the binary, self-connected, row-normalized
graph in float32. `capo_pre` uses `data/capo_pre/capo_pre.h5ad`, its actual
preprocessed moment source, with matching identifiers and abundances. The
forward min-max transform must reproduce saved `Ms`/`Mu` exactly. Per axis,
`y = a*x + b` gives `v_x = v_y/a`; translation affects state only. Other
methods use their existing numerical units.

Simulation runs write `simulation/velovi_cossim.csv`, comparing raw restored
velocities with `velocity_s_true`/`velocity_u_true` using the historical norm
mask (`>1e-6` for both vectors), without CBDir's velocity cleanup. The compact
CosSim table combines these 70 veloVI values with the other eight methods by
dataset/method/gene; [the results guide](results/README.md) describes its fields
and source coverage.

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

These rendering commands read the committed CBDir and CosSim CSVs listed
above. `--output-root` selects the export destination. The resulting
`scores_used.csv` and `provenance.json` identify the actual score inputs.

### Common-gene comparison

Check the common-gene tables without H5AD inputs:

```bash
python evaluation/gene/refresh_common_gene_summary.py
```

Add `--write --output-root outputs/common-genes` to regenerate them separately.
The finite intersections contain 6/6/6/7/2/4/4 genes in the real-dataset order
listed above; macro means give each dataset equal weight.

### Display and source boundaries

[Source Data](../../data/SOURCE_DATA.md) includes native phase arrays and runnable
entries for the Figure 5A/E raw phase portraits and Figure 5D overviews.

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
Simulation colors follow the Figure 4 progression-bin gradient from early to late.

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

Run from the repository root:

```bash
python evaluation/gene/plot_selected_edge_panels.py \
  --artifact-root /path/to/velocity-artifacts \
  --output-root outputs/selected-edge-panels
python evaluation/gene/plot_real_comparisons.py \
  --output-root outputs/real-gene-comparisons
```

Figure 3 uses native phase arrays, current transition-level Gene scores and
`results/selected_edge_cell_controls.csv` for the same annotated transitions.
Each export includes plotted score rows.
The Figure 5 scatters weight dataset means equally; the adjacent branching
boxplot and heatmap retain their pooled-gene distribution. The S9 panel uses
the finite nine-method intersection in `results/common_genes/`.
