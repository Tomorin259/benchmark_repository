# RNA-only MultiVelo phase-assignment diagnostics

Post-hoc controls of the simulation benchmark. The diagnostic tests
whether a gene-level induction-only prefit constrains downstream fitting in a
way that favors directional recovery in the evaluated branching subsets.

## Scientific scope

MultiVelo 0.1.0 in RNA-only mode uses the RNA transcription--splicing--degradation
equations with chromatin activity fixed to one. A gene-level partial-trajectory
decision precedes continuous parameter fitting and cell-by-gene time assignment.
An `on` decision excludes repression candidates in `calculate_dist_and_time`
before continuous parameter fitting and cell-specific time assignment.

The principal intervention keeps the automatic `on` prefit and changes only
the conditional excluding repression candidates. It retains initial parameters,
fitting inputs, neighbor graphs, optimizer settings and output processing.
Final rates and times may change as downstream results of refitting. The
installed package is untouched; the driver verifies that the in-process function
differs by one AST condition. Native `off` and `complete` paths are checked for
numerical equality before fitting.

Default refits reproduce archived velocities and scores exactly for the 40
branching dataset--gene pairs. Opening candidates decreases CBDir-gene and
CosSim in all 40 pairs, while the fitting objective decreases in all 40. Each
condition retains 2,000 valid velocity comparisons per gene. This isolates the
candidate constraint as a contributor in these frozen realizations, rather
than attributing the difference to richer lineage-specific rates.

Other retained controls answer distinct questions:

| Control | Fixed | Changed |
| --- | --- | --- |
| scVelo induction-only | fitted rates and scaling | selected phase and corresponding smoothed induction time |
| MultiVelo gate open | `on` prefit, initialization, fitting and output settings | availability of repression candidates during refitting |
| MultiVelo no output smoothing | fitted parameters, times and fitting loss | final velocity smoothing/rescaling |
| MultiVelo complete trajectory | data and ordinary fitting settings | trajectory scope, gap handling and refitted parameters |

`complete` varies trajectory scope, gap handling and refitted parameters
together. The separate constant-positive-direction comparison characterizes
the induction-heavy ground-truth direction distribution.

## Retained materials

- [Native candidate-gate outputs](https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark-results/resolve/v1.0.0/v1.0.0/diagnostics/multivelo/candidate_gate_outputs.tar.gz?download=true):
  four H5ADs for the default and gate-open conditions in normal + branching
  and hybrid. The archive contains `multivelo/{branching,all}_{refit_default,on_allow_repression}.h5ad`.
  [The member index](../../data/result_member_manifest.tsv) maps them to their
  original locators and scientific hashes in `provenance.json`.
- `results/scvelo_controls.csv`: 490 full-precision rows for seven valid controls
  across the 70 evaluated dataset--gene pairs.
- `results/multivelo_refits.csv`: 250 non-duplicated rows: 70 default, 70 complete,
  70 no-output-smoothing, and 40 gate-open results.
- `results/condition_summary.csv`: per-subset gene means and sample SDs. The SD
  describes variation across the evaluated genes.
- `results/candidate_gate_effects.csv`: paired counts for the two branching subsets.
- Branch-level, pilot cell-state, native scVelo state, source-input identity,
  default reproduction and initialization checks remain under `results/`.
- `inputs/phase_case.npz`: exact coordinates and velocity arrays for the selected
  BoostALL_5 case. The archive has no pickled objects and can be read safely with
  `allow_pickle=False`; a clean checkout can re-render without large H5ADs.
- `provenance.json`: frozen input identities, historical execution records,
  versions, code and result hashes, and external native output locators.
- `figures/`: three phase portraits, one paired-score row and an early--late key.
- `tables/tableS9.tex`: the booktabs source for Supplementary Table S9, with
  four conditions and four-decimal scores for the two branching subsets.
  MURK, smoothing, complete-trajectory and paired-count results remain in
  the full CSVs. The image files in `tables/` retain that full diagnostic
  preview; they are not the compact manuscript Table S9.
- `materials_manifest.json`: rendering inputs, sampled arrows, output
  checksums, font and resolution for Supplementary Figure S15 A/B and Table S9.
- `../../figures/assembly/suppfig15/`: the retained editable author assembly,
  composed A/B materials and their assembly entry point. The final author
  export and its source map are registered in the paper-level figure manifest.

CBDir-gene uses the fixed reference-30 assets owned by `evaluation/gene/`.
CosSim excludes either predicted or true vectors with norm at or below 1e-6;
the scVelo archive also retains the fixed truth-eligible-cell control. Means and
SDs are computed from unrounded gene scores. Figure/table display uses four
decimal places, positive values have no leading plus sign, and source CSVs
retain full precision. All controls use the registered simulation inputs.

S15A uses BoostALL_5 from normal + branching. Its annotated CosSim averages
all 2,000 valid cells (default 0.7534, Gate open 0.2074). The supplementary
text's branch-specific values (0.8037 and -0.2920) use the 1,000 BranchB cells;
they are a different statistic, not alternative whole-gene scores. In Gate
open, 722 of those 1,000 cells receive repression assignments. The gene-level
`on` classification is retained; no cells are removed by the intervention.

## Reproduce

Run from repository root. Use Python 3.9.24, MultiVelo 0.1.0, scVelo 0.3.3 and
NumPy 2.0.2 for the verified fitting controls. The expected MultiVelo source
hash is in `provenance.json`. Rerendering needs NumPy, pandas, AnnData,
Matplotlib, SciPy and scikit-learn; it performs no fitting.
The vector/raster A/B assembly additionally requires pypdf and Pillow.
Select those packages in a local environment of your choice. For refits,
`--artifact-root` is the velocity directory containing the relative H5AD
paths in `evaluation/gene/config.py`; `--repository-root` is this checkout.
The rendering and table commands use the committed compact inputs by default.

```sh
python supplement/multivelo/plot_materials.py --output-dir outputs/multivelo-materials
python supplement/multivelo/plot_materials.py --paired-only --output-dir outputs/multivelo-paired-materials
python supplement/multivelo/summarize_controls.py --output-root outputs/multivelo-summary
python supplement/multivelo/summarize_controls.py --table-only --output-root outputs/multivelo-table
python figures/assembly/suppfig15/assemble_materials.py \
  --output-dir outputs/multivelo-assembly
python supplement/multivelo/run_scvelo_controls.py \
  --repository-root . --artifact-root /path/to/velocity-artifacts \
  --output outputs/multivelo-scvelo-controls
python supplement/multivelo/run_multivelo_controls.py \
  --repository-root . --artifact-root /path/to/velocity-artifacts \
  --candidate-gate-control --branching-only --output outputs/multivelo-gate-control
python supplement/multivelo/run_multivelo_controls.py \
  --repository-root . --artifact-root /path/to/velocity-artifacts \
  --all-genes --output outputs/multivelo-complete-control
python supplement/multivelo/run_multivelo_controls.py \
  --repository-root . --artifact-root /path/to/velocity-artifacts \
  --all-genes --no-output-smoothing --output outputs/multivelo-smoothing-control
python supplement/multivelo/inspect_phase_rule.py \
  . /path/to/velocity-artifacts outputs/phase_rule_inputs.csv
python3 -m unittest discover -s tests -p test_multivelo_diagnostics.py -v
```

Use a fresh output directory for each fit. Model inputs resolve through
`evaluation/gene/config.py`; reference neighborhoods and annotation edges are
not recomputed. `--case-gene BoostALL_5` can replace `--branching-only` for the
single-case gate pilot. Each diagnostic reads the registered simulation inputs.

`--paired-only` exports one `paired_controls` row containing both branching
subsets, following the accepted paired-dot layout in
`supplement/veloae/plot_pca_summary.py`. Filled points show default outputs;
open points show modified conditions (scVelo induction-only and MultiVelo
gate open), keyed as Default/Modified at the lower right. Method names are neutral axis labels; method colors are reserved
for data marks. `scores_used.csv` retains all 160 plotted gene-condition rows.

Figure rendering and summary verification use the compact files in Git.
Refitting additionally requires the frozen method-native H5ADs identified in
[provenance.json](provenance.json). Its `repository:`, `artifact-root:` and
`environment:` locators resolve against user-supplied roots; they are not
public download URLs. Large native outputs remain external to Git.
The provenance record also documents excluded controls and reconstruction attempts.
