# Kinetic-reference gene usage tables

These tables record the current expanded reference-gene evaluation at repository commit 2402e9a7ff300162ae9d573348a535d502c77b28.

- rpe1_reference_gene_usage.csv: 250 retained genes from the 528 genes in Battich et al. Table S1, in the frozen benchmark candidate order. Figure 6 B/C display flags identify the 19 selected examples.
- organoid_reference_gene_usage.csv: 270 retained differentiation-associated genes from the 295 selected in Table S3 and profiled in Table S4, in the frozen candidate order. The same reference set is used in both branches. Figure 6 branch display flags identify the 6 and 8 selected examples. The 25 unavailable selected genes and 6 housekeeping controls are outside this retained table.

## Column meanings

- gene_symbol: original reference symbol, retaining human/mouse capitalization.
- ensembl_gene_id: Table S1 gene_id for RPE1; original canonical organoid input var.ID for organoid, also checked against all 37 registered TSvelo gene-ID mappings.
- reference_order: original reference-table order (Table S1 or S4), not a score rank.
- alpha_reference and gamma_reference: TRUE indicates a retained reference gene for that parameter. Benchmark alpha corresponds to kappa (transcription) in the original scEU tables; gamma is degradation.
- branch1_in_reference_scope / branch2_in_reference_scope: TRUE indicates membership in that branch's reference evaluation universe. This is independent of each method's output availability.
- figure6_*_display: TRUE identifies a selected gene displayed in the corresponding Figure 6 parameter panels. These flags do not describe streamline recomputation or the entire reference evaluation.
- *_status: per-method, per-parameter and, for organoid, per-branch Spearman score status:
  - valid: an evaluation row with a finite Spearman score; included in the available-gene boxplot and numerical summary.
  - undefined: an evaluation row exists, but its Spearman score is non-finite; excluded from numerical summaries, never replaced by zero.
  - unavailable: no evaluation row for this method/parameter/gene in the current reference universe.

The tables describe current parameter evaluation, not every gene used during model training. cellDancer and DeepVelo have alpha and gamma columns; DeepKINET has gamma only; TSvelo and Dynamo(label data) have alpha only in this comparison. Omitted method/parameter combinations were not evaluated. Dynamo_label denotes the registered Dynamo(label data) comparator. Organoid TSvelo undefined scores arise from constant predicted profiles (17 per branch).

## Sources and checks

Original references: Battich et al., Science (2020), https://doi.org/10.1126/science.aax3072, Tables S1/S3/S4. Rate trajectories remain in the original reference CSVs; these gene tables do not duplicate them.

Repository sources: evaluation/kinetics_recovery/results/expanded_reference/gene_scope.csv; results/central_results/gene_metrics.csv; results/expanded_reference/tsvelo_organoid_gene_metrics.csv; results/expanded_reference/tsvelo_organoid_gene_mapping.csv; and figures/fig6/plot_data/ selected-gene records, all within evaluation/kinetics_recovery/.

All 21 method/parameter/branch coverage rows match coverage_summary.csv exactly. All 18 finite-gene comparison memberships match common_gene_membership.csv exactly. Gene IDs and symbols are unique within each table. File checksums and source checksums are recorded in kinetics_reference_gene_usage_manifest.json.
