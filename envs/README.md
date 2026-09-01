# Environment evidence

Supplementary Table S1 reports the package and selected dependency versions used in the benchmark. `manifests/environments.tsv` connects those statements to Conda environments still present on the audited server and records two separate audit dimensions: whether the installed entry point imports, and whether a clean public specification has been prepared.

Most assignments are inferred from matching environment names and package tuples. They are not lock files. The Phase 2 audit found an importable candidate for each of the 16 configuration environments; these checks were package or entry-point import smoke tests and did not rerun model training. DeepVelo GCN imports with a NumPy/SciPy compatibility warning. The SDEvelo environment used for the original check is not import-compatible, so the audited C17 candidate is the working `sdevelo_py310` environment; its release specification still requires curation. DeepKINET, veloAE, and DeepVelo_scVelo also require specification cleanup because of version drift, editable or historical metadata, or missing revision information.

For that reason, Phase 2 does not export the current Conda environments. A raw `conda env export` would capture unrelated build strings, editable installs, local paths, and historical contamination without proving which packages produced the paper results. A passing import smoke test must not be read as a complete reproduction claim.

In `environments.tsv`, `evidence_status` describes historical assignment evidence, `runtime_status` describes the observed import smoke result, and `release_spec_status` describes public-release readiness. These statuses are intentionally independent.

## Requirements for a reviewed specification

Each future environment file must:

- identify exactly one or an explicitly shared set of configuration IDs;
- pin the benchmarked method version or upstream commit;
- reproduce the manuscript-reported Python and key dependency versions, or document a validated deviation;
- avoid local filesystem paths and unreviewed editable installs;
- record GPU/CUDA constraints where required;
- install successfully in a clean environment;
- pass the corresponding runner smoke test;
- preserve third-party license notices without vendoring unapproved assets.

The audited Conda paths are provenance evidence only and will not appear in the final public execution interface.
