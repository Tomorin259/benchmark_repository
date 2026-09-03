# Frozen publication-repository design

**Status:** frozen baseline, 2026-09-03

**Scope:** final public repository for the RNA velocity benchmark

**Internal baseline:** private audit repository at commit `5284bdd`

This document fixes the repository architecture for the remainder of the
publication-code curation. The top-level design is no longer an open question.
Concrete filenames may be adjusted only when required by an upstream method,
an unavoidable scientific distinction, or a deletion/merge discovered during
implementation.

## Purpose

The public repository records the shortest path needed to understand and
reproduce the experiments reported in the paper. It is not a copy of the
server project, a history of the release audit, or a general-purpose benchmark
framework.

The governing contract is:

```text
methods/                         produce method-result H5AD files
evaluation/ and ensemble/       consume those result files
results/                         contains compact canonical paper results
figures/                         contains generated publication artifacts only
large data and H5AD files       live outside GitHub
```

Complexity is allowed when it comes from the scientific experiment. Complexity
introduced only by repository curation is not part of the public design.

## Four reader-facing layers

```text
1. Method reproduction
   methods/
   envs/
   preprocessing/

2. Benchmark evaluation
   evaluation/cell/
   evaluation/gene/
   evaluation/kinetics/
   simulation/

3. Added utility
   ensemble/

4. Publication artifacts
   results/
   figures/
   data/
   docs/
```

A reader should be able to answer four questions immediately:

- How were method results produced? `methods/`
- How were those results evaluated? `evaluation/`
- What is the optional post-hoc ensemble? `ensemble/`
- Where are the reported numerical results? `results/`

## Final public layout

```text
rna-velocity-benchmark/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
│
├── data/
│   ├── README.md
│   └── datasets.tsv
│
├── preprocessing/
│   ├── README.md
│   ├── prepare_counts_inputs.py
│   └── <only preprocessing scripts actually used by the paper>
│
├── methods/
│   ├── README.md
│   ├── configurations.tsv
│   ├── scvelo.py
│   ├── multivelo.py
│   ├── unitvelo.py
│   ├── dynamo.py
│   ├── tivelo.py
│   ├── cell2fate.py
│   ├── tfvelo.py
│   ├── veloae.py
│   ├── velovi.py
│   ├── deepkinet.py
│   ├── tsvelo.py
│   ├── celldancer.py
│   ├── deepvelo.py
│   ├── sdevelo.py
│   ├── deepvelo_scvelo.py
│   └── graphvelo.py
│
├── envs/
│   ├── README.md
│   ├── evaluation.yml
│   └── <one reviewed environment per dependency-compatible method group>
│
├── evaluation/
│   ├── cell/
│   │   ├── README.md
│   │   ├── config.json
│   │   ├── metrics.py
│   │   └── plots.py
│   ├── gene/
│   │   ├── README.md
│   │   ├── config.json
│   │   ├── metrics.py
│   │   └── plots.py
│   └── kinetics/
│       ├── README.md
│       ├── config.json
│       ├── metrics.py
│       └── plots.py
│
├── simulation/
│   ├── README.md
│   ├── config.json
│   └── generate.py
│
├── ensemble/
│   ├── README.md
│   ├── config.json
│   └── ensemble.py
│
├── results/
│   ├── README.md
│   ├── cell/
│   ├── gene/
│   ├── kinetics/
│   └── ensemble/
│
├── figures/
│   ├── README.md
│   ├── main/
│   └── supplementary/
│
├── docs/
│   ├── benchmark_design.md
│   ├── result_contract.md
│   ├── reproducibility.md
│   └── figure_map.md
│
├── tests/
│   ├── test_metadata.py
│   └── test_kinetics_evaluation.py
│
└── .github/
    └── workflows/
        └── ci.yml
```

Directories are added only when they contain release material. Empty
`.gitkeep` scaffolding is not part of the final repository.

## Method-reproduction contract

The sixteen method scripts represent the twenty benchmark configurations.
Multiple configurations of the same software remain command-line modes of one
script when they share the same official implementation, for example the
three scVelo modes, two UniTVelo modes, and two GraphVelo input variants.

Each method script must:

1. expose the necessary dataset/key mapping;
2. keep the upstream tutorial sequence visible in the main execution path;
3. accept explicit input and output locations;
4. validate the result keys needed by downstream evaluation;
5. write a method-native H5AD result.

The public repository does not contain `BaseRunner`, plugin, registry,
provider, workflow, or universal-runner layers. It also does not vendor full
upstream repositories. A helper or additional method file is allowed only
when it owns a real upstream or scientific distinction.

## Input and preprocessing contract

`data/README.md` explains access. `data/datasets.tsv` records the dataset ID,
paper label, role, accession or stable source, canonical processed input, and
checksum when redistribution is possible.

Original data that already have a stable GEO, Zenodo, Figshare, or official
package source are linked rather than redundantly archived. Our own external
archive is limited to:

1. canonical processed inputs actually used by the paper;
2. canonical method-result H5AD files;
3. reference or intermediate data without a stable substitute.

It is not intended to become a complete data lake. Raw server trees, Conda
installations, logs, exploratory outputs, obsolete result generations, and
checkpoints that are not required for a publication claim remain remote-only.

For count-requiring methods, `prepare_counts_inputs.py` owns the shared input
contract. In particular, `capo_pre`, `mousebrain`, `scnt_neuron`, and
`rpe1_gen` retain the benchmark cell/gene/metadata contract while their
`spliced` and `unspliced` layers are aligned from native integer counts. These
fixed objects are not passed through a second
`scv.pp.filter_and_normalize()` call. Inputs beginning from full raw matrices
receive the declared HVG preprocessing at most once.

Before an integer dtype conversion, the script reports and checks finite,
non-negative, integer-valued counts. The conversion must preserve every
value. Count inputs use the `_counts.h5ad` suffix.

## Result-H5AD boundary

Method H5AD files retain their native layer and embedding names. The benchmark
does not rewrite every result into a synthetic universal schema. Each
evaluation module translates the heterogeneous keys at its own input boundary
through a small explicit mapping in `config.json`.

At minimum, a canonical result preserves cell and gene identifiers and stores
the method-specific quantities consumed by the corresponding evaluation.
Large canonical result files are distributed externally; compact numerical
tables used in the paper are committed under `results/`.

## Kinetic-parameter recovery module

The public repository does not use the internal name `Case 3`. Its final name
and layout are fixed as:

```text
evaluation/kinetics/
├── README.md
├── config.json
├── metrics.py
└── plots.py
```

`metrics.py` retains the scientific data flow:

```text
method H5AD and reference tables
  -> parameter extraction and alignment
  -> RPE1 or organoid processing
  -> per-gene Spearman and cosine scores
  -> component, system, and method summaries
  -> compact CSV and trajectory outputs
```

The following method-specific extraction is required and remains explicit:

- cellDancer: `celldancer_alpha` and `celldancer_gamma`;
- DeepVelo: `cell_specific_alpha` and `cell_specific_gamma` from the
  unspliced-result family;
- DeepKINET: `degradation_rate`;
- Dynamo: sparse `cell_wise_alpha` from the labelled RPE1 H5AD;
- TSvelo: `alpha` from two RPE1 lineage H5AD files, with duplicate cells
  averaged.

The established RPE1 and organoid selection, binning, smoothing,
normalization, branch direction, negative clipping, and gene lists also remain
part of the scientific code.

The public configuration removes `source_generation`,
`historical_regression`, server generation names, unused modes, and absolute
paths. `validate_historical_extraction.py` stays in the private audit history
and is not a public runtime entry point. Its enduring claims are covered by a
small synthetic regression test.

`plots.py` consumes metric outputs and writes plot-ready CSV files and
descriptively named panels such as `rpe1_trajectories.png` or
`organoid_branch1_boxplot.png`. It does not write another manifest and does not
claim to assemble an exact manuscript Figure PDF unless it actually does so.

## Results and figures

`results/` contains only compact, canonical, machine-readable values reported
or directly plotted in the paper. Full method H5ADs and large intermediate
matrices are external artifacts.

`figures/` is an output directory, not a code layer. Plotting logic stays with
the experiment that owns it in `evaluation/*/plots.py` or `ensemble.py`.
Final main and supplementary figures may be committed when they are genuine
outputs of the documented workflow. Manual assembly, where unavoidable, is
stated in `docs/figure_map.md`.

## Documentation boundary

Module-specific instructions belong in that module's `README.md`. The `docs/`
directory is reserved for cross-module facts only:

- overall experimental design;
- the common H5AD/result contract;
- cross-method reproducibility and randomness policy;
- mapping from manuscript panels to code, inputs, and outputs.

No public `phase`, `audit`, `provenance`, `historical-path`, or agent-oriented
documentation hierarchy is recreated.

## Tests and CI

Tests protect cheap, high-value publication contracts; they are not a software
product coverage programme. The baseline suite checks metadata consistency,
absence of private absolute paths, and the heterogeneous kinetic extraction
and numerical processing on synthetic data.

CI runs only these lightweight checks. GPU training, all-method environment
creation, and complete benchmark execution remain documented HPC workflows
and are not required in GitHub Actions.

## Internal/public boundary

The current private repository remains the audit record. It may keep manifests,
source hashes, path aliases, historical comparisons, and decision notes used
to establish the public result. Those materials are evidence for curation, not
the public execution interface.

The final public repository is created from a clean history. It must not
contain `/mnt/...`, `/Users/...`, `velocity:`, `project:`, rerun directory
names, local environment prefixes, or historical result locations.

## Allowed future changes

The top-level architecture is frozen. Later changes are allowed only when:

1. an upstream method genuinely requires an additional method-specific file;
2. an experiment has an unavoidable scientific special case;
3. implementation shows that two files should be merged or an unused
   directory/file should be deleted.

Changes justified only by making the repository look more framework-like are
out of scope. In particular, new `runner/`, `adapter/`, `core/`, `workflow/`,
`registry/`, `provenance/`, `phase/`, or `audit/` layers are not part of the
frozen public design.

## Implementation order

1. Freeze and retain the current private audit baseline.
2. Curate `evaluation/kinetics/` as the first complete public module.
3. Use that module as the direct template for cell- and gene-level evaluation.
4. Curate method scripts and reviewed environments one method at a time.
5. Freeze canonical processed inputs and result H5AD files externally.
6. Commit compact results and publication figures.
7. Create the clean-history public repository and complete citation/license
   metadata.

No later implementation phase reopens the repository's top-level design.
