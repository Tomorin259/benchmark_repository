# Bone marrow paired TIvelo root refits

## Scope

Paired root-sensitivity experiment using the archived TIvelo bone marrow output.
Both root conditions use matched inputs, optimization settings and seeds.

The intervention starts **after main-path orientation**: compare the archived
Mono_1 DTI root with a biologically informed HSC_1 DTI root. It changes the
root-dependent child dictionary, DPT and directed neighbors together. The
upstream main-path orientation and root-selection outputs are retained.

## Controlled design

- 5,780 cells and 473 archived selected genes; identical cell/gene ordering.
- Frozen Mu/Ms, undirected KNN/connectivities, archived PAGA and diffusion map.
- Original-root replay exactly matches the archived child dictionary and d_nn.
- DTI t1=0.1, t2=2.0; t2=2 reproduces the archived graph.
- Upstream root_select='connectivities' selects the **first cell** in the root
  cluster; the same selection rule is used in both conditions.
- Paired seeds 20261001, 20261002, 20261003; 100 epochs each; MSE, both velocity
  components, constrain=True, alpha_1=1, alpha_2=0.1, batch_size=1024,
  filter_genes=False.
- Native velocity graphs/projections regenerated independently after refitting.
- Canonical reference-30 cell evaluator, frozen evaluation moments/neighbors,
  dataset embedding and configured edges; official TIvelo formula parity checks.
- The three paired fits measure variability across optimization seeds.

## Results

Means over the three paired seeds:

| Metric | Mono_1 root refit | HSC_1 root refit | Paired mean delta |
| --- | ---: | ---: | ---: |
| CBDir | 0.109349 | 0.662020 | +0.552671 |
| BATC | 0.283570 | 0.745330 | +0.461759 |
| High-dimensional CBDir | 0.298701 | 0.485683 | +0.186982 |
| TransCosine | 0.410822 | 0.539704 | +0.128882 |
| Consistency | 0.563642 | 0.581049 | +0.017406 |
| ICCoh | 0.613680 | 0.621251 | +0.007571 |

Every paired delta is positive for all six metrics and all three seeds.
CBDir paired deltas range from +0.551091 to +0.554339; BATC paired deltas range
from +0.460534 to +0.462999.

| Evaluated edge | Boundary source cells | Mono_1 refit CBDir | HSC_1 refit CBDir |
| --- | ---: | ---: | ---: |
| HSC_1 -> HSC_2 | 1007 | -0.483490 | +0.541545 |
| HSC_1 -> Ery_1 | 157 | -0.008433 | +0.597987 |
| Ery_1 -> Ery_2 | 291 | +0.819970 | +0.846528 |

Boundary-cell counts are unchanged across all seven scored cases. Most of the
overall CBDir improvement comes from the two HSC_1-origin edges; erythroid
Ery_1 -> Ery_2 directionality remains positive. High-dimensional CBDir also
improves, so the rescue is not confined to the stored 2D projection.

The archived fit's CBDir=0.192291 and BATC=0.296348 differ from the original-root
refits. Intervention effects are measured by the within-seed difference between
the two root refits. Scores recomputed on the unchanged
archive match the canonical tables (maximum absolute difference <3e-17).

DPT versus the independent inferred Palantir pseudotime:

| Cell subset | Cells | Mono_1 root Spearman | HSC_1 root Spearman |
| --- | ---: | ---: | ---: |
| Myeloid-associated subset | 4224 | -0.694711 | +0.891176 |
| Erythroid-associated subset | 2125 | +0.796731 | +0.920618 |

Palantir supplies an independently inferred pseudotime ordering. Both lineage
subsets include HSC_1 cells.

## Interpretation

This controlled sensitivity experiment supports a root-dependent DTI mechanism
contributing to poor directionality in the real bone marrow dataset: changing
the DTI root alters ordering/directed neighbors and consistently improves the
configured cell-level directions after otherwise matched refitting.

The HSC_1 root uses biological prior information. Residual cycles and
disconnections in the root-dependent child graph are retained in both
conditions.

## Files

- `cell_scores.csv`: archived output plus six refits, all six cell metrics.
- `paired_metric_deltas.csv`: 18 within-seed, within-metric comparisons.
- `paired_summary.csv`: paired means, seed SD and paired-delta ranges.
- `edge_scores.csv`, `edge_summary.csv`: CBDir, high-dimensional CBDir,
  TransCosine and BATC candidate edge details. Candidate edge BATC is not the
  overall BATC aggregation and should not be presented as interchangeable.
- `directed_neighbor_pairs.csv`: cross-cluster directed-neighbor counts.
- `dpt_reference_correlations.csv`: the two prespecified ordering checks.
- `provenance.json`: portable source and evaluator hashes, settings, package versions,
  graph replay, output hashes and read-only input verification.
- `manifest.json`: generated in each run's output directory, recording inputs,
  fits and evaluation progress for `--resume`.
- `native_outputs/*.h5ad`: six refitted outputs retained remotely for plotting.
- `evaluations/` and `logs/`: per-case evaluation inputs, edge evidence,
  parity checks and runtime logs generated alongside the refits.

Native X_umap equals the archived X_tsne in this dataset: future visualizations
should identify these coordinates as t-SNE, not newly computed UMAP.

## Reproduction

Use your own TIvelo 0.1.4 installation with the package versions in
`provenance.json`; [Cell dependencies](../../../docs/environments.md#evaluation)
describe the separate preparation worker. The fitting interpreter also runs
the official TIvelo scorer. Choose local paths, for example:

```bash
TIVELO_PYTHON=/path/to/tivelo-analysis/bin/python
CELL_METRICS_PYTHON=/path/to/cell-evaluation/bin/python
ARTIFACT_ROOT=/path/to/artifacts
VELOCITY_ROOT="$ARTIFACT_ROOT/velocity"
REPOSITORY_ROOT=/path/to/benchmark_repository
REFERENCE_MOMENTS_ROOT=/path/to/reference-moments
NEW_EXPERIMENT_DIR="$VELOCITY_ROOT/outputs/tivelo-bonemarrow-control"
```

The evaluation stage converts native H5AD paths to Cell `velocity:` locators.
Both `--source` and the new output directory must therefore be below
`<artifact-root>/velocity/`. Refits are written to
`<output-dir>/native_outputs/`, with evaluations and logs in sibling
subdirectories. Use a **new output directory** for an independent rerun;
`--resume` resumes only hash-verified outputs.

```bash
"$TIVELO_PYTHON" -u supplement/tivelo_gene/bonemarrow_root_control.py \
  --artifact-root "$ARTIFACT_ROOT" \
  --repo-root "$REPOSITORY_ROOT" \
  --source "$VELOCITY_ROOT/results/TIvelo/tivelo_bonemarrow.h5ad" \
  --reference-moments-root "$REFERENCE_MOMENTS_ROOT" \
  --indicators-python "$CELL_METRICS_PYTHON" \
  --output-dir "$NEW_EXPERIMENT_DIR" \
  --seeds 20261001 20261002 20261003 --threads 8
```

Canonical input result files and benchmark tables were not modified. Native
outputs and their SHA256 hashes were individually checked after completion.
