# Exploratory gene–cell directionality cases

Early exploratory analysis of 15 gene cases and 28 displayed gene–transition
pairs. The cases compare individual-gene velocity directions with cell-level
directions on annotated transitions. Full benchmark results are in the
[Cell](../../evaluation/cell/README.md) and [Gene](../../evaluation/gene/README.md)
modules; the formal Figure 3 panels have their own
[selected-edge rendering entry](../../evaluation/gene/README.md#selected-edge-and-cross-dataset-source-panels).

## Materials and reproduction

| Material | Contents |
| --- | --- |
| [figures/](figures/) | Masked spliced–unspliced phase portraits and case manifests |
| [analysis/](analysis/) | Paired observations, aggregations, correlations, ranks and sign comparisons |
| [Source mapping](analysis/source_mapping_audit.csv) | Method H5ADs, velocity-layer aliases and case-source pairing |
| [results/](results/) | Matched Cell results and the additional Diaph3 transition |
| [Analysis report](analysis/gene_cell_edge_analysis.md) | Descriptive results and aggregation units |
| [Analysis provenance](analysis/provenance.json) | Input, generator and output checksums |

Run from the repository root with a fresh output directory. Summary analysis
uses the Python standard library and committed tables:

```bash
python supplement/gene-kinetic/analysis/gene_kinetic_analysis.py \
  --repo-root . --output-dir outputs/gene-edge-analysis
```

Phase rendering additionally needs the registered method H5ADs, reference
labels and the [Gene plotting environment](../../docs/environments.md#evaluation).
The artifact root is the velocity directory containing `results/`:

```bash
python supplement/gene-kinetic/gene_kinetic_pilot.py \
  --artifact-root /path/to/velocity-artifacts \
  --case dentate_ncald --output-dir outputs/gene-cases
```

Use `--case all` to render all 15 cases. Each case exports a PNG, PDF and manifest
with its source locators, cluster mask, annotated scores and display diagnostics.

Check the saved analysis, source checksums and all 252 score annotations without
H5AD inputs:

```bash
python -m unittest discover -s tests -p test_gene_case_results.py
```

## Score pairing and aggregation

Gene scores use the current [reference-30 edge table](../../evaluation/gene/results/gene_cbdir_edge_details_reference30_all.csv).
The Cell controls use the recorded matched outputs for DeepVelo, UniTVelo and
veloVI, and the canonical Cell outputs for the other six methods. MultiVelo
uses its official normalized native spliced velocity and matching projection.
The [current edge and source tables](results/current/) supply these controls.
Source H5ADs, velocity layers, method configurations and source-cell sets are recorded
in the mapping table. Valid pairs require finite scores; missing values remain
missing.

Genes sharing a dataset, source, target and method are averaged before
edge-level association analysis. This avoids counting shared cell directions
as independent gene replicates. Method-specific and within-edge-centered
correlations remain descriptive analyses of this small case set.
`batc_candidate_source_mean` is the mean candidate cosine before branch
maximization; `branch_winner_share` records how often that edge wins. These
candidate diagnostics have a different aggregation from the final BATC score.

veloVI numerical annotations use restored pre-min–max coordinates and velocity
units, following the [Gene scoring contract](../../evaluation/gene/README.md#velovi-coordinate-units).
Phase portraits retain method-native coordinates and velocities.

## Case definitions and biological context

The primary-case subset and additional cases below retain their original
biological questions and control transitions. Coverage is the number of
methods with an available gene score out of nine.

| Group | Dataset and genes | Transitions and controls | Biological context |
| --- | --- | --- | --- |
| Primary: late kinetic switching | Gastrulation: Hba-x (7/9), Gypa (8/9), Alas2 (7/9), Fam210b (8/9); Klf3, Hemgn and Fech (9/9 each) | Erythroid2 → Erythroid3; Erythroid1 → Erythroid2 is the early control | The original mouse MURK analysis uses an increased u–s slope in Erythroid3 relative to Erythroid2. These genes occur in its MURK list; Hba-x and Fam210b have illustrated phase/experimental examples. [Barile et al.](https://doi.org/10.1186/s13059-021-02414-y) |
| Primary: granule-cell maturation | DentateGyrus_lamanno: Ncald (9/9) | ImmGranule1 → ImmGranule2; ImmGranule2 → Granule; CA branches provide additional comparisons | Ncald is a cellDancer branching gene; functional work links it to dentate-gyrus neurogenesis. [cellDancer](https://doi.org/10.1038/s41587-023-01728-5), [functional study](https://doi.org/10.3389/fnmol.2019.00019) |
| Primary: progenitor exit | DentateGyrus_lamanno: Diaph3 (8/9) | nIPC → Nbl1; Nbl1 → Nbl2 is the secondary check | Diaph3 is a cellDancer branching gene. Knockout evidence concerns progenitor division/neurogenesis in cortex. The added nIPC → Nbl1 comparison has 304 reference-neighbor source cells. [cellDancer](https://doi.org/10.1038/s41587-023-01728-5), [knockout study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8102060/), [transition annotation](https://tivelo.readthedocs.io/en/latest/notebooks/notebooks/DentateGyrus2.html) |
| Primary: cell-cycle context | capo_pre: E2F7 and FAM111A (9/9 each) | G1 → S and S → G2M; G2M → G1 supplies cycle context | E2F7 regulates cell-cycle genes; FAM111A promotes replication initiation. [E2F7](https://pmc.ncbi.nlm.nih.gov/articles/PMC3333892/), [FAM111A](https://pmc.ncbi.nlm.nih.gov/articles/PMC10551639/) |
| Additional: migration/localization | DentateGyrus_lamanno: Dab1 (8/9) | CA → CA1-Sub | CA1 migration phenotypes provide functional context rather than direct evidence of unusual Dab1 RNA kinetics. [Dab1 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10089055/) |
| Additional: retinal expression | retina: Celf4 (9/9) | Neuroblast → RGC and Neuroblast → AC/HC | In situ data support stage-dependent expression; the transitions have 52 and 41 evaluable source cells. [Celf4 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3825781/) |
| Additional: cortical migration | mousebrain: Flrt2 (9/9) | V-SVZ → Upper Layer | FLRT2 supplies a cortical migration signal. Its own RNA velocity is a separate measurement; Upper Layer → Deeper Layer is outside this case's directional interpretation. [FLRT2/Unc5D study](https://doi.org/10.1038/emboj.2011.189) |
| Additional: terminal function | pancreas: G6pc2 (9/9) | Pre-endocrine → Beta; Pre-endocrine → Alpha is the control | Knockout evidence provides β-cell functional context. [G6pc2 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7331801/) |

## Phase masks and display

| Case group | Cluster mask | Annotated transitions |
| --- | --- | --- |
| Late kinetic switching | Erythroid1, Erythroid2, Erythroid3 | Erythroid1 → Erythroid2; Erythroid2 → Erythroid3 |
| Granule-cell maturation | Nbl2, ImmGranule1, ImmGranule2, Granule | ImmGranule1 → ImmGranule2; ImmGranule2 → Granule |
| Progenitor exit | nIPC, Nbl1, Nbl2 | nIPC → Nbl1; Nbl1 → Nbl2 |
| Cell-cycle context | G1, S, G2M | G1 → S; S → G2M |
| Migration/localization | CA, CA1-Sub, CA2-3-4 | CA → CA1-Sub |
| Retinal expression | Neuroblast, RGC, AC/HC | Neuroblast → RGC; Neuroblast → AC/HC |
| Cortical migration | V-SVZ, Upper Layer, Deeper Layer | V-SVZ → Upper Layer |
| Terminal function | Pre-endocrine, Beta, Alpha | Pre-endocrine → Beta; Pre-endocrine → Alpha |

Each case has nine method panels and a shared cluster legend in a 2 × 5 layout.
Method titles show transition-level Gene CBDir, with `N/A` for unavailable
scores. Diaph3 nIPC → Nbl1 uses the independent
[supplementary calculation](results/diaph3_nipc_to_nbl1.csv); its missing
DeepVelo_scv gene remains `N/A`.

The [phase API](../../evaluation/gene/plotting/phase_portrait.py) uses native
`Ms`/`Mu` and maps the method's two velocity components to
`velocity`/`velocity_u`. Display settings remain a 25 × 25 grid, density
percentile 15, log1p arrow lengths, reference quantile 0.90, reference length
0.05, maximum length 0.10, point size 12 and alpha 0.6, with hidden axes/ticks.
These settings affect display; numerical annotations come from the recorded
score tables.
