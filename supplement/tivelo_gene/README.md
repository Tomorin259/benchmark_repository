# TIvelo root and ordering diagnostics

Mechanism sensitivity analyses for branching and hybrid simulations, with a
paired bone marrow root intervention. These biologically informed refits explain
failures; they must not replace original benchmark scores or method rankings.

## Evidence

| Evidence | Data | Entry script |
| --- | --- | --- |
| Archived gene scores | `results/` | `analyze_tivelo_gene.py` |
| Main-path selection and orientation | `stage_audit/` | `stage_audit.py` |
| Root decision and upstream diagnostics | `root_decision_audit/`, `upstream_scvelo_audit/` | `root_decision_audit.py`, `upstream_scvelo_audit.py` |
| Fixed-expression ordering intervention | `ordering_trigger_audit_20260929_threads1/` | `tivelo_ordering_trigger_audit.py` |
| Secondary-path DPT | `branch_dpt_export_20260930/` | `export_tivelo_branch_dpt.py` |
| Simulation root pairs used in panel C | `root_stage_intervention/` | `root_stage_intervention.py` |
| Independent optimizer-seed control | `root_stage_intervention_seed20260927/` | `root_stage_intervention.py` |
| Bone marrow root pairs used in panel D | `bonemarrow_root_control_20261001/` | `bonemarrow_root_control.py` |

Other retained audits (`truth_graph_root_audit/`, `transition_flux_audit/`,
`orientation_order_audit/`, `orientation_path_intervention/`,
`orientation_stat_intervention/`, `root_mode_intervention/`, and the earlier
`results/root_intervention*`) test different interventions, not duplicates of
the final root-stage pairs. Their top-level scripts remain available; use
`--help` for explicit inputs and outputs.

Branching initially selects BranchA_9, then BranchB_9 after orientation
correction, rather than an early root. Hybrid initially selects BranchA_0;
path-wise DPT subsequently changes the root to BranchB_9. Its fixed-expression
ordering audit changes the u/s direction score from -19.594465 to 77.043071.
These failures occur at different stages in the recorded runs.

| Final paired simulation | Metric | Corrected-path root | Known early root |
| --- | --- | ---: | ---: |
| Branching (20 genes) | CBDir-gene | -0.282741 | 0.513026 |
| Branching (20 genes) | CosSim | -0.432665 | 0.480011 |
| Hybrid (40 genes) | CBDir-gene | -0.400341 | 0.686810 |
| Hybrid (40 genes) | CosSim | -0.403071 | 0.672195 |

Bone marrow compares Mono_1 and HSC_1 DTI roots across three paired optimizer
seeds. Mean cell-level CBDir changes from 0.109349 to 0.662020, and BATC from
0.283570 to 0.745330. The intervention changes the child graph, DPT and directed
neighbors together; it does not isolate DPT alone or establish real-data
ground-truth gene velocities. See the experiment README for controls.

## Final materials and reproduction

`figures/` is the sole active material bundle: compact A/B
`ab_right_legend_material.{pdf,png}`, simulation C
`tivelo_root_control_gene_metrics.{pdf,svg,png}`, bone marrow D
`bonemarrow_root_control_row.{pdf,svg,png}`, source panels, plotted CSVs and
historical source manifests. Paper composites, panel letters and editable PPT
belong in [`figures/assembly/suppfig13/`](../../figures/assembly/suppfig13/).
The assembly manifest records current hashes; historical generator hashes are
not relabeled as hashes of the consolidated scripts.

Run from repository root in the TIvelo environment (anndata, scanpy,
matplotlib, numpy, pandas). PDF assembly additionally needs pypdf, reportlab
and the original Arial TTF. Use a fresh export directory.

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

Large H5AD outputs are external artifacts retained on the execution host, not
committed binaries. The small manifests retain their checksums. Public path
prefixes `repository:`, `velocity:`, `artifact-root:` and `environment:` refer
to user-supplied roots. No fitting is needed to re-render when the native
artifacts have been materialized.

The raw bone marrow `manifest.json`, evaluation workspaces and logs remain on
the execution host and are ignored by Git. `provenance.json` is its portable
published counterpart; the D renderer accepts either. Historical `--resume`
requires the original runtime manifest and matching inputs. Independent reruns
should use a new output directory.

Cleanup removed superseded preview layouts, duplicate renderers, empty failed
runs and caches from the active module. Result CSVs, native outputs and approved
panel binaries are unchanged.
