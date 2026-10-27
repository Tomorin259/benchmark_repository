# TIvelo root and ordering controls

Supplementary Figure 13 examines root-dependent direction in frozen branching
and hybrid simulations and a paired bone marrow intervention.
The paired refits measure sensitivity to root selection and trajectory ordering.

## Final materials

| Panel | Inputs / source | Output material |
| --- | --- | --- |
| A/B, lineage and ordering | [stage_audit/](stage_audit/), [branch_dpt_export_20260930/](branch_dpt_export_20260930/) | [figures/](figures/), compact `ab_right_legend_material` |
| C, simulation root pairs | [root_stage_intervention/](root_stage_intervention/), [independent seed control](root_stage_intervention_seed20260927/) | `tivelo_root_control_gene_metrics` and plotted CSVs |
| D, bone marrow pairs | [bonemarrow_root_control_20261001/](bonemarrow_root_control_20261001/README.md) | `bonemarrow_root_control_row` and source assets |

[Editable assembly and final PDF](../../figures/assembly/suppfig13/README.md)
carry panel letters and whole-figure layout.
Input/result identities and current material hashes remain in the experiment
and assembly manifests.

## Re-render

Run from the repository root with AnnData, Matplotlib, NumPy, pandas and
SciPy for the saved-data plots. PDF assembly additionally needs pypdf, reportlab
and the original Arial TTF. Set `VELOCITY_ROOT`, `EXPORT_DIR` (fresh) and
`ARIAL_TTF` to your local paths.

```bash
MODULE=supplement/tivelo_gene
python "$MODULE/figures/scripts/plot_tivelo_full_lineage_shared_dpt.py" \
  --artifact-root "$VELOCITY_ROOT" --audit-root "$MODULE" \
  --branch-data "$MODULE/branch_dpt_export_20260930" \
  --output-dir "$EXPORT_DIR/ab_assets"
python "$MODULE/figures/scripts/assemble_ab_material.py" \
  --materials-root "$EXPORT_DIR" --output-dir "$EXPORT_DIR" \
  --font-file "$ARIAL_TTF"
python "$MODULE/figures/scripts/plot_tivelo_compact_control.py" \
  --data-root "$MODULE" --output-dir "$EXPORT_DIR"
python "$MODULE/figures/scripts/plot_tivelo_bonemarrow_dpt_assets.py" \
  --data-root "$MODULE/bonemarrow_root_control_20261001" \
  --output-dir "$EXPORT_DIR" --seed 20261001
```

Native H5ADs are external. Locators prefixed `repository:`, `velocity:`,
`artifact-root:` or `environment:` use caller-supplied roots.
Bone marrow `provenance.json` records the published experiment. Independent
reruns use a fresh output root; each run saves its own manifest and logs for
subsequent `--resume`.

The lineage plot's `--artifact-root` is the velocity directory containing
the configured native H5AD paths. Bone marrow refitting instead uses the
Cell parent-root convention and separate fitting/metric Python roles;
its [local setup](bonemarrow_root_control_20261001/README.md#reproduction)
defines those paths. Refits require TIvelo 0.1.4 and the recorded scientific
dependencies, separately from saved-data rendering.

## Controls and interpretation

The simulation pairs compare corrected-path roots with known early roots:
branching CBDir-gene/CosSim changes from -0.282741/-0.476315 to
0.513026/0.535458; hybrid changes from -0.400341/-0.463156 to
0.686810/0.745974. The paired bone marrow comparison changes roots Mono_1/HSC_1
across three optimizer seeds. Mean Cell CBDir changes from 0.109349 to 0.662020
and BATC from 0.283570 to 0.745330.
The root intervention jointly recomputes the child graph, DPT and directed
neighbors before velocity fitting.

Simulation gene-wise CosSim uses the shared noise-free ODE reference in
`evaluation/gene/simulation/truth/`, matched by cell and gene names.
The full-expression upstream and truth-graph diagnostics require the
all-gene reference generated with `make_truth_reference.py --all-genes`;
pass its directory with `--truth-root`. Their high-dimensional direction
comparisons and root/graph interventions are distinct from per-gene CosSim.

For a paired replay of retained root interventions, pass the corresponding
`root_stage_gene_scores.csv` with `root_stage_intervention.py --root-table`
and its recorded seed. `stage_audit.py --root-table` likewise fixes the root
used for velocity scoring while retaining the inferred upstream diagnostics.
Without this option, the scripts run the upstream root-selection procedure;
its selected root can vary even when the velocity optimizer seed is fixed.

Supporting root-decision, upstream-scVelo, ordering-trigger, flux and
orientation interventions remain linked by their manifests and top-level
scripts. Each records its intervention stage and controlled inputs. In the recorded runs, branching and hybrid fail at different ordering
stages; the fixed-expression hybrid direction score changes from -19.594465
to 77.043071. Use the relevant script's `--help` and a fresh output directory
for a new diagnostic.
