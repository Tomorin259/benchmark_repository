# Final publication-repository design

**Status:** final and frozen, 2026-09-03

**Scope:** final public repository for the RNA velocity benchmark

**Internal baseline:** private audit repository at commit `8cff281`

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
methods/                         reproduce method-native result H5AD files
evaluation/ and ensemble/       own each experiment's analysis, compact results, and panels
artifacts.tsv                    locates and verifies large external inputs and result H5AD files
data/ and docs/                  contain cross-cutting dataset metadata and documentation
```

Complexity is allowed when it comes from the scientific experiment. Complexity
introduced only by repository curation is not part of the public design.

## Development reading rule

Before modifying a method, data-preparation, evaluation, simulation, or ensemble
script, the developer must read this document and the README for the module being
changed. The existing entry point must then be read from input to output before
new helpers, files, configuration options, or abstractions are proposed.

The private curation repository keeps this requirement in its root `AGENTS.md` so
that automated development sessions encounter it before editing. The durable rule
is this document; `AGENTS.md` is only its short development entry point.

## Scientific-script coding standard

The repository contains research reproduction and data-analysis scripts. They are
not production services, security-critical model-serving code, or a reusable
benchmark framework. Code should therefore resemble a carefully reviewed analysis:
the data source, scientific transformations, metric calculation, and written
outputs should be visible in execution order.

### Start direct

Begin with the shortest sequential implementation that exposes the scientific
workflow. A dataset-level or figure-level function may own a complete analysis when
that makes the experiment easier to read. Do not begin by designing classes,
interfaces, runners, adapters, registries, factories, dependency injection, or a
universal execution framework.

Plain Python functions, dictionaries, pandas objects, NumPy arrays, and AnnData
objects are the default vocabulary. A class is justified only when the analysis has
real persistent state or a resource lifecycle that is clearer as an object; grouping
stateless helper functions is not sufficient justification.

### A function must earn the split

Extract a function when at least one of these is true:

1. it names a scientific operation that should be understood or tested on its own;
2. it is used in more than one meaningful place;
3. it translates a genuinely different external file or method-result format;
4. it enforces an invariant that would otherwise be duplicated or easy to miss;
5. keeping it inline would make the owning experimental step materially harder to
   read.

A function used once is not automatically wrong, but a single-use function that
only renames one expression, forwards arguments, reads one JSON/CSV file, converts
one array, or constructs one path should normally be inlined. Do not split one
linear analysis into a chain of mechanical helpers that must all be opened to learn
what happens.

Method- or dataset-specific duplication is acceptable when the two workflows have
different scientific meanings. Do not merge RPE1 and organoid processing, or two
method-specific storage formats, merely because a few statements look similar.

### Keep configuration scientific

`config.json` contains paths, dataset keys, method layer mappings, gene sets, and
parameters that define or select the experiment. It must not become an inventory of
private implementation details.

Do not add a configurable strategy, mode, aggregation rule, reader type, or boolean
for hypothetical flexibility. If the code supports exactly one behavior and that
behavior is an established method-specific rule, implement it directly and explain
it in the module README. Keep a config field only when a current dataset or method
actually varies, or when changing the value is a legitimate experiment choice.

### Preserve visible scientific flow

For method scripts, the upstream tutorial sequence should remain recognizable in
the main execution path. For evaluation scripts, a reader should be able to follow:

```text
read inputs -> extract quantities -> dataset-specific processing
-> align references -> calculate metrics -> summarize -> write outputs
```

Names should describe scientific meaning, such as `process_rpe1`,
`select_organoid_branch`, or `calculate_metrics`. Avoid generic names such as
`Manager`, `Handler`, `Processor`, `Service`, or `run_step` when they conceal which
analysis is being performed.

Input validation should protect reproducibility claims: required keys, shapes,
identifiers, finite values, count properties, and assumptions used by the metric.
Do not add enterprise-style validation layers, custom exception hierarchies, or
defensive branches for states that the script cannot receive through its documented
interface.

### Keep I/O and effects obvious

File reads and writes should remain close to the experimental step that owns them,
or behind one small helper when a repeated naming convention or path invariant is
being enforced. Avoid multi-stage writer/helper chains and output registries.

Outputs have deterministic names documented in the module README. Public scripts
do not create timestamps, manifests, provenance objects, run classes, or logging
frameworks unless a publication contract explicitly requires them. A few clear
`print` statements are sufficient for long-running scientific scripts.

### Test behavior, not helper choreography

Tests should protect publication-relevant behavior through the script entry point
or through a small number of genuine scientific transformations. Do not retain or
export a helper solely because a test imports it. Tests must not require the
implementation to stay fragmented into the same call sequence.

Full GPU workflows remain outside CI. Small synthetic arrays or AnnData files are
preferred for extraction, alignment, normalization, and metric regression checks.

### Refactoring acceptance rule

Before simplifying a working script, record the observable contract: output files,
matrix shapes and labels, selected cells and genes, transformation parameters, and
reported metrics. After the change, compare those outputs rather than private
function calls.

The default refactor must leave the final module the same size or smaller in total
files, functions, and lines. An increase is accepted only for a present scientific
capability, a real external-format boundary, or a regression check that protects a
paper claim. “Cleaner architecture”, future flexibility, and consistency with a
software framework are not sufficient reasons.

During review, ask:

- Can the scientific sequence be understood without jumping through mechanical
  helpers?
- Which single-use wrappers can be inlined?
- Does every config option have a current second value or a real user choice?
- Does every abstraction remove more reader burden than it adds?
- Can obsolete compatibility, audit, or output-writing scaffolding be deleted?
- Did the change preserve the numerical and file-output contract?

## Paper-aligned reader path

```text
1. Method reproduction
   methods/
   envs/

2. Benchmark evaluation
   evaluation/cell/
   evaluation/gene/
   evaluation/gene/simulation/
   evaluation/parameter_recovery/
   evaluation/summary/

3. Added utility
   ensemble/

4. Cross-cutting records
   data/
   docs/
   artifacts.tsv
```

A reader should be able to answer four questions immediately:

- How were method results produced? `methods/`
- How were those results evaluated? `evaluation/`
- What is the optional post-hoc ensemble? `ensemble/`
- Where are the reported values and panels? Beside the experiment that owns
  them, under its `results/` and `figures/` subdirectories.
- Where are large inputs and method-result H5AD files? In external storage,
  indexed by `artifacts.tsv`.

## Final public layout

```text
rna-velocity-benchmark/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── artifacts.tsv
│
├── data/
│   ├── README.md
│   └── datasets.tsv
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
│   ├── README.md
│   ├── cell/
│   │   ├── README.md
│   │   ├── config.json
│   │   ├── metrics.py
│   │   ├── plots.py
│   │   ├── results/
│   │   └── figures/
│   ├── gene/
│   │   ├── README.md
│   │   ├── config.json
│   │   ├── metrics.py
│   │   ├── plots.py
│   │   ├── simulation/
│   │   │   ├── README.md
│   │   │   ├── config.json
│   │   │   └── generate.py
│   │   ├── results/
│   │   │   ├── simulation/
│   │   │   └── real/
│   │   └── figures/
│   ├── parameter_recovery/
│   │   ├── README.md
│   │   ├── config.json
│   │   ├── metrics.py
│   │   ├── plots.py
│   │   ├── results/
│   │   └── figures/
│   └── summary/
│       ├── README.md
│       ├── build.py
│       ├── results/
│       └── figures/
│
├── ensemble/
│   ├── README.md
│   ├── config.json
│   ├── geometric_median.py
│   ├── ensemble.py
│   ├── benchmark.py
│   ├── plots.py
│   ├── results/
│   └── figures/
│
├── docs/
│   ├── benchmark_design.md
│   ├── result_contract.md
│   ├── reproducibility.md
│   └── figure_map.md
│
├── tests/
│   ├── test_metadata.py
│   ├── test_cell_evaluation.py
│   ├── test_gene_evaluation.py
│   └── test_parameter_recovery.py
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

A method that requires raw counts must additionally make that requirement
visible in its own script and README. For DeepKINET and cell2fate this includes
selecting the count-bearing input, documenting any dataset-specific mapping,
checking count invariants, and preventing accidental repeated preprocessing.
Methods that do not require raw counts do not inherit this machinery.

The public repository does not contain `BaseRunner`, plugin, registry,
provider, workflow, or universal-runner layers. It also does not vendor full
upstream repositories. A helper or additional method file is allowed only
when it owns a real upstream or scientific distinction.

## Data and method-input contract

`data/README.md` explains access. `data/datasets.tsv` records the dataset ID,
paper label, role, accession or stable source, canonical processed input, and
checksum when redistribution is possible.

Original data that already have a stable GEO, Zenodo, Figshare, or official
package source are linked rather than redundantly archived. Our own external
archive is limited to:

1. canonical processed inputs actually used by the paper;
2. canonical method-result H5AD files;
3. reference or intermediate data without a stable substitute.

The root `artifacts.tsv` is the portable index for those external files. It
records a stable artifact ID, dataset or method, role, external URI, checksum,
and the H5AD keys required by the consuming script. It must not contain a
server-specific working path.

It is not intended to become a complete data lake. Raw server trees, Conda
installations, logs, exploratory outputs, obsolete result generations, and
checkpoints that are not required for a publication claim remain remote-only.

There is no shared `preprocessing/` package and no generic
`prepare_counts_inputs.py`. Data preparation that belongs to an upstream
method remains in that method's script, where its order relative to the
official tutorial is visible.

Each method that requires raw counts owns its count-bearing input selection,
dataset-specific alignment, finite/non-negative/integer checks, lossless dtype
conversion, and protection against repeated preprocessing. The affected
datasets and key mappings are documented beside that method rather than imposed
on methods that do not consume counts. Count-bearing inputs use the
`_counts.h5ad` suffix and remain external artifacts indexed by `artifacts.tsv`.

## Result-H5AD boundary

Method H5AD files retain their native layer and embedding names. The benchmark
does not rewrite every result into a synthetic universal schema. Each
evaluation module translates the heterogeneous keys at its own input boundary
through a small explicit mapping in `config.json`.

At minimum, a canonical result preserves cell and gene identifiers and stores
the method-specific quantities consumed by the corresponding evaluation.
Large canonical result files are distributed externally and indexed in the
root `artifacts.tsv`, including their URI, checksum, and required keys. Compact
numerical tables are committed under the `results/` directory of the
evaluation or ensemble module that owns them.

## Parameter-recovery module

The public repository does not use the internal name `Case 3`. Its final name
and layout are fixed as:

```text
evaluation/parameter_recovery/
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

`plots.py` consumes metric outputs and writes plot-ready CSV files to the
module's `results/` directory and descriptively named panels such as
`rpe1_trajectories.png` or `organoid_branch1_boxplot.png` to its `figures/`
directory. It does not write another manifest and does not claim to assemble an
exact manuscript Figure PDF unless it actually does so.

## Module-owned results and figures

There is no repository-wide `results/` or `figures/` directory. Each evaluation
or ensemble module owns the compact, canonical, machine-readable values and
generated panels for its scientific question. This keeps extraction,
processing, metric definitions, numerical outputs, and visualization together.

Full method H5ADs and large intermediate matrices remain external. Compact CSV
files and selected reproducible panels may be committed. Large exploratory
PNGs, redundant format variants, and manually assembled manuscript composites
are not required repository contents. Cross-module manuscript assembly, such
as Figure 5 panels and the Figure 6 summary, is recorded in
`docs/figure_map.md`; `evaluation/summary/` owns the code and compact outputs
for integrated benchmark summaries.

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
absence of private absolute paths, and the heterogeneous cell-, gene-, and
parameter-recovery extraction and numerical processing on synthetic data.

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

Status recorded on 2026-09-04:

| Order | Work package | Status |
| ---: | --- | --- |
| 0 | Establish the paper-aligned repository skeleton | Complete |
| 1 | Freeze and retain the private audit baseline | Complete |
| 2 | Curate `evaluation/parameter_recovery/` as the first runnable evaluation module | Implemented; final result promotion and release packaging remain |
| 3 | Use the same module ownership pattern for cell- and gene-level evaluation | Not started; existing audit metadata and compact tables are source material only |
| 4 | Curate direct method scripts and reviewed environments one method at a time | Not started |
| 5 | Freeze canonical processed inputs and result H5AD files externally | Not started; `artifacts.tsv` is an incomplete placeholder index |
| 6 | Commit compact module-owned results and selected reproducible panels | Not started outside the existing audited cell metric table |
| 7 | Create the clean-history public release and complete citation/license metadata | Not started; the current GitHub repository remains private curation space |

The next implementation targets are cell- and gene-level evaluation. Their
relative order can follow the paper-curation workflow. Status changes update
this table; they do not reopen the frozen top-level architecture.

No later implementation phase reopens the repository's top-level design.
