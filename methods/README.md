# Method configuration registry

This directory will contain only benchmark-owned thin runners and adapters after they are reviewed. Upstream method packages, cloned repositories, checkpoints, and generated results are not vendored here.

The manuscript evaluates 20 **configurations**, not 20 distinct software packages. The three scVelo modes, two UniTVelo modes, and two GraphVelo input variants are separate benchmark rows. The nine-method gene-level subset, the metabolic-labeling Dynamo positive control, and the optional geometric-median ensemble do not add rows to this registry.

The machine-readable sources are:

- `manifests/method_configurations.tsv`: manuscript identity, aliases, configuration, output keys, and release strategy;
- `manifests/runner_inventory.tsv`: audited candidate runners and their input/output contracts;
- `manifests/environments.tsv`: manuscript-reported versions and evidence for the historical environments.

## Evidence boundary

`candidate_runner` means that a script structurally matches the method and at least part of the observed output contract. It does not mean the script has been proven to regenerate all 14 final H5AD files. Several scripts have partial active dataset lists, dataset-specific wrappers, stale source aliases, or an unrecorded copy step between method-local output and the curated `results/` tree.

The canonical cell-level output contract is defined by the 280-row render manifest: each of 20 configurations contributes one result for each of 14 datasets. The method-specific velocity and embedding keys recorded in `method_configurations.tsv` describe those final inputs. Forebrain and a few method-specific projections have explicit exceptions and must not be forced into a false universal schema.

## Promotion criteria

A candidate becomes a public runner only after it has:

1. a manifest-driven 14-dataset input mapping;
2. no commented or manually edited active-dataset list;
3. explicit output paths and layer keys matching the publication contract;
4. a reviewed environment specification tied to an upstream version or commit;
5. at least one smoke test and one numerical regression check;
6. no bundled upstream source or asset whose redistribution terms are unresolved.

Until those conditions are met, this Phase 2 registry is an audit record rather than an executable method suite.
