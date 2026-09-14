# Repository design

This document defines the stable design rules of the benchmark repository. It is normative for repository structure, ownership, interfaces, and scientific-script style. Current completion status, available artifacts, and known gaps belong in the top-level `README.md`, not here.

## Purpose

The repository supports audit and reproduction of the paper's benchmark results. It should make the following chain explicit:

```text
registered method and dataset
    -> method-native result
    -> benchmark evaluation
    -> compact result table or panel
    -> paper figure or supplementary table
```

The repository is deliberately a collection of inspectable scientific workflows, not a general workflow framework. A contributor should be able to identify the scientific inputs, transformations, and outputs without first learning a large internal abstraction layer.

## Development reading rule

Before changing a module:

1. Read this document and the module's `README.md` in full.
2. Inspect the actual entry script, registry rows, and tests that define the relevant contract.
3. Check `git status` and preserve unrelated worktree changes.
4. Keep module-specific instructions in the module, and update cross-module documentation only when a shared contract changes.

If prose conflicts with executable code or a machine-readable registry, first determine which reflects the paper's intended configuration. Fix the drift rather than adding another compatibility description.

## Ownership model

Each fact should have one authoritative home.

| Location | Responsibility |
| --- | --- |
| `README.md` | Public entry point, current repository capabilities, limitations, and reproduction routes |
| `methods/` | Method registry, paper-facing names, method launchers, and method-native outputs |
| `envs/` | Release environment definitions and installation checks when those environments are published |
| `data/` | Dataset identity, source description, canonical-input contract, and simulation generation |
| `artifacts.tsv` | Logical IDs, checksums, stable locations, and materialization metadata for external artifacts |
| `evaluation/` | Metric definitions, evaluation workflows, compact results, and module-level panels |
| `ensemble/` | Ensemble definitions, inputs, outputs, evaluation, and ensemble-specific figures |
| `figures/` | Paper-level composite figures and deterministic assembly metadata or scripts |
| `docs/` | Cross-module contracts and paper-to-repository maps |
| `provenance/` | Durable historical evidence that cannot be expressed in code, registries, or artifact metadata |
| `tests/` | Fast contract checks for names, schemas, mappings, and lightweight deterministic behavior |

The owning module keeps its own compact tables and generated panels. The root `figures/` directory is reserved for paper-level composites and their assembly records; it is not a second home for module results.

## Scientific-script coding standard

### Start direct

A scientific script should begin as a direct, linear account of the computation:

1. declare scientific configuration;
2. load and validate inputs;
3. perform the scientific transformation;
4. save the result and essential metadata.

Do not introduce a class, registry object, generic runner, plug-in system, or deep helper stack before there are multiple real callers that need it.

### A function must earn the split

Create a helper when it provides at least one of these benefits:

- it is reused by multiple real call sites;
- it names a coherent scientific operation;
- it isolates nontrivial validation or file-format handling;
- it has a stable contract that can be tested independently;
- it substantially reduces repeated code without hiding the scientific sequence.

Do not create helpers merely to move a few lines away from the main flow. Thin wrappers that only rename a library call, forward arguments, or mutate hidden state make audit harder.

### Keep configuration scientific

Parameters that define the paper's computation should remain visible near the entry point or in a small, explicit configuration file. Examples include dataset IDs, method modes, neighbor definitions, randomization policy, metric choices, and output names.

Avoid large generic configuration systems when a short table or command-line interface communicates the experiment more clearly. Operational details may be configurable, but they must not obscure the formal scientific setting.

### Preserve visible scientific flow

The main path should expose the order in which data are transformed. A reviewer should not have to jump through several files to discover which matrix was evaluated, how cells were aligned, or which neighborhood was used.

When a library tutorial is the authoritative method workflow, keep the repository script structurally close to that tutorial. Add only the validation, path handling, metadata capture, and output materialization needed for benchmark integration.

### Keep I/O and effects obvious

- Inputs and outputs must be explicit.
- File writes should occur at visible boundaries.
- Randomness must be described: either fixed and recorded, or intentionally unfixed as part of the experiment definition.
- Server-specific absolute paths must not appear in public code or documentation.
- Failed validation should stop with an actionable error instead of silently falling back.
- Generated files must not overwrite an unrelated scenario or method output.

### Test behavior, not helper choreography

Tests should protect scientific and repository contracts: registry completeness, canonical names, table schemas, method-to-script mappings, artifact hashes, figure inputs, and small deterministic calculations. They should not freeze an incidental helper decomposition.

Full method training, GPU workloads, large downloads, and complete figure regeneration do not belong in routine CI. Provide smoke tests or validation modes for those workflows where practical.

### Refactoring acceptance rule

A refactor is acceptable only when it preserves the scientific result contract, keeps provenance visible, and leaves the final implementation no larger or more indirect without a demonstrated need. Temporary migration code must have a clear removal condition.

## Cross-cutting contracts

### Registries and names

`methods/configurations.tsv` is the authoritative registry of benchmark method configurations and their paper-facing names. `data/datasets.tsv` is the authoritative registry of dataset IDs and display names.

Operational IDs are stable machine identifiers used in paths and joins. Display names are presentation-layer labels used in figures, tables, and prose. Code must join through the registries rather than maintaining private rename dictionaries. Prose documents should explain the rule and link to the registry, not copy the full list and allow it to drift.

### Method reproduction

Each registered configuration must resolve to a runnable method entry point or an explicitly documented source-backed launcher. The module documentation must state:

- required input and expected schema;
- scientific parameters or modes used in the paper;
- invocation pattern;
- output location and format;
- validation or smoke-test route;
- upstream release or commit when external source is required.

Multiple configurations may share a script when the mode selection is explicit and auditable. The repository does not require one universal runner.

### Data and artifacts

Small registries, checksums, configuration files, and compact results belong in Git. Large inputs and model outputs may live in an external artifact store.

Every external artifact needed for reproduction should have a logical ID, SHA256 checksum, media or object type where useful, and a stable URI once published. Code should resolve logical IDs or user-supplied roots; it must not depend on a private server path.

The dataset registry defines identity. Dataset-specific documentation defines the canonical input contract, including required layers, annotations, and transformations. If the project distributes a canonical real-data input directly, it need not also provide a redundant raw-data construction workflow, but its source and checksum must remain auditable.

Simulation code must record the simulator source/version, scenario parameters, output mapping, and randomness policy. Intentional stochasticity is valid when it is stated explicitly; it must not be mislabeled as byte-exact reproduction.

### Result boundary

Method execution produces method-native result files. Evaluation code consumes those files through explicit mappings and writes compact metric tables, edge-level audit data where required, or figure-ready summaries.

Method runners must not silently redefine benchmark metrics. Evaluation modules must not become hidden method launchers. Shared validation and file-format utilities are acceptable when their contracts are small and stable.

Native H5AD outputs should retain the layers, embeddings, graph objects, velocity fields, and metadata required by downstream evaluation or audit. Optional exploratory state may be omitted, but omission must not prevent verification of a reported result.

### Results and figures

An evaluation module owns the compact numerical results and panel renderers for its scientific question. Paper-level figure assembly belongs in `figures/`.

The assembly manifest is the authoritative mapping from a paper panel to its source artifact. It should record source paths or artifact IDs, output names, and hashes where appropriate. A final composite may be stored directly or built deterministically; either route must make panel ordering and provenance unambiguous.

Do not maintain hand-written inventories of panel counts or completion percentages in multiple documents. Generate inventories from manifests when needed, and keep the current public status only in the top-level `README.md`.

### Documentation and provenance

Documentation is layered by responsibility:

- the top-level `README.md` tells a new reader what is available now and where to start;
- module `README.md` files describe local inputs, commands, outputs, validation, and limitations;
- `docs/` defines shared evaluation, result, naming, reproduction, and figure-mapping contracts;
- `provenance/` records durable evidence for exceptional historical or manual steps.

Provenance is not a project diary. Do not commit agent notes, temporary audit reports, private absolute paths, duplicated status tables, or narratives that can be regenerated from Git history and machine-readable metadata.

### Tests and CI

Routine CI should be fast and environment-light. It should verify repository contracts without downloading large artifacts or importing every method environment. Typical checks include:

- registry uniqueness and completeness;
- referenced script and result paths;
- table schemas and canonical labels;
- artifact metadata and checksums for tracked files;
- figure manifest consistency;
- small deterministic calculations and import smoke tests where dependencies permit.

Expensive end-to-end runs should be documented and may be executed as release validation outside routine CI.

## Change policy

This design is stable but not frozen. Add or reorganize structure only when required by a scientific distinction, an upstream method constraint, a durable file-format boundary, or a testable cross-module contract.

Avoid framework-style infrastructure, duplicate registries, speculative compatibility layers, and documents whose only purpose is to repeat current status. Prefer deleting obsolete descriptions over marking them historical inside active documentation.

Every public change should leave the path from registered input to reported result at least as easy to inspect as before.
