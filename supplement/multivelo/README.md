# RNA-only MultiVelo phase-assignment diagnostics

Post-hoc controls of the frozen simulation benchmark. These analyses do not
replace canonical velocities, scores or method rankings. The diagnostic tests
whether a gene-level induction-only prefit constrains downstream fitting in a
way that favors directional recovery in the evaluated branching subsets.

## Scientific scope

MultiVelo 0.1.0 in RNA-only mode uses the RNA transcription--splicing--degradation
equations with chromatin activity fixed to one. A gene-level partial-trajectory
decision precedes continuous parameter fitting and cell-by-gene time assignment.
An `on` decision excludes repression candidates in `calculate_dist_and_time`;
it is not a replacement for parameter estimation or a lineage-specific model.

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

`complete` is a whole-fitting-scheme sensitivity control, not a pure candidate
ablation. The constant-positive-direction comparison is retained separately as
an exploratory characterization of the induction-heavy true direction
distribution; it is not a biological method or a primary ranking entry.

## Retained materials

- `results/scvelo_controls.csv`: 490 full-precision rows for seven valid controls
  across the 70 evaluated dataset--gene pairs.
- `results/multivelo_refits.csv`: 250 non-duplicated rows: 70 default, 70 complete,
  70 no-output-smoothing, and 40 gate-open results.
- `results/condition_summary.csv`: per-subset gene means and sample SDs. The SD
  describes genes, not repeated simulations or optimizer seeds.
- `results/candidate_gate_effects.csv`: paired counts for the two branching subsets.
- Branch-level, pilot cell-state, native scVelo state, source-input identity,
  default reproduction and initialization checks remain under `results/`.
- `inputs/phase_case.npz`: exact coordinates and velocity arrays for the selected
  BoostALL_5 case. The archive has no pickled objects and can be read safely with
  `allow_pickle=False`; a clean checkout can re-render without large H5ADs.
- `provenance.json`: frozen input identities, historical execution records,
  versions, code and result hashes, and external native output locators.
- `figures/`: three phase portraits, one paired-score row and a short early--late key. No
  theoretical flowchart, branch-state fraction plot or fitted curve is included.
- `tables/tableS9.tex`: the booktabs source for Supplementary Table S9, with
  four-decimal means in A/B and paired-decrease counts in C. The image files
  in `tables/` are numerical previews, not the typeset manuscript table.
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
retain full precision. The frozen independent realizations are not regenerated.

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

```sh
python supplement/multivelo/plot_materials.py --output-dir outputs/multivelo-materials
python supplement/multivelo/plot_materials.py --paired-only --output-dir outputs/multivelo-paired-materials
python supplement/multivelo/summarize_controls.py --output-root outputs/multivelo-summary
python supplement/multivelo/archive_results.py --materials-only \
  --materials-dir outputs/multivelo-materials
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
single-case gate pilot. The input generation policy remains intentionally
unseeded; no new simulation is generated by any diagnostic command.

`--paired-only` exports one `paired_controls` row containing both branching
subsets, following the accepted paired-dot layout in
`supplement/veloae/plot_pca_summary.py`. Filled points show default outputs;
open points show modified conditions (scVelo induction-only and MultiVelo
gate open), keyed as Default/Modified at the lower right. Method names are neutral axis labels; method colors are reserved
for data marks. `scores_used.csv` retains all 160 plotted gene-condition rows.

Large native outputs remain on the execution host and are ignored by Git.
`repository:`, `artifact-root:` and `environment:` locators are relative to
user-supplied roots, not public download URIs. `archive_results.py` is the
explicit compact-result import boundary for the original completed runs;
it validates and de-duplicates defaults, rather than selecting favorable scores.
Runtime manifests may contain execution-local paths; their public counterpart
is the portable `provenance.json`.

To archive identified completed runs on the execution host, use
`archive_results.py --diagnostics-root /path/to/completed-diagnostics
--artifact-root /path/to/velocity-artifacts`. After checking a fresh rendering,
add `--materials-dir /path/to/accepted-rendering` to publish its figure materials,
shared key and Table S9 source under this module. Input and output hashes must
match before the copy is accepted.

For rendering-only updates, use `--materials-only` as above. It validates
the current plotting script, source inputs and accepted export hashes without
re-importing or reformatting any frozen scientific CSV. Clean-checkout figure
rendering and summary verification use the compact files in Git. Repeating the
model fits additionally requires the frozen method-native H5ADs identified in
`provenance.json`; unresolved external locators are not public download URLs.

Hand-reconstructed MultiVelo velocities failed output reproduction and are
excluded from attribution. Two attempted scVelo controls disabling steady-state
candidates raise an upstream IndexError and have no reported scores. The
excluded attempts remain documented in provenance, not silently treated as
negative scientific results. General superiority, true rate recovery and
behavior on real datasets are not established by these controls.
