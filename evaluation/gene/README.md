# Gene-level velocity evaluation

Nine method configurations are evaluated on seven real datasets and the three
fixed simulation scenarios. Native `Ms`, `Mu` and spliced/unspliced velocities
are retained; the canonical CBDir comparison uses the fixed reference-30 graph.

## Results and selection

| Source | Use |
| --- | --- |
| [results/gene_cbdir_reference30_all.csv](results/gene_cbdir_reference30_all.csv) | Canonical `new_cbdir` scores |
| [results/gene_cbdir_edge_details_reference30_all.csv](results/gene_cbdir_edge_details_reference30_all.csv) | Transition-level results |
| [results/simulation/gene_cbdir_cossim_reference.csv](results/simulation/gene_cbdir_cossim_reference.csv) | Frozen simulation CosSim; its historical CBDir column is a different source |
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
is a subset. `VELOCITY_ARTIFACT_ROOT` can supply the artifact root.
The evaluator also reads the historical comparison-score files named in its
source/provenance; the saved H5ADs alone do not supply those comparisons.

`--neighbor-mode native` is a controlled comparison, with a separate output
root. Results include per-dataset scores, edge details, neighbor assets and
audit/provenance records. The committed reference-neighbor NPZs allow graph
inspection without the original H5ADs.

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
