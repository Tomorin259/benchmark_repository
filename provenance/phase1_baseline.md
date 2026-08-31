# Phase 1 provenance baseline

Inventory frozen on **2026-08-31**.

## Audited sources

- Working project root: `/mnt/second19T/zhaozelin/velocity`
- Paper snapshot root: `/Users/tomorin/Downloads/test1_1`
- Paper source of record: `merged.tex`
- `merged.tex` SHA-256: `172bf2d84af52d8b7f005f588fbe81a44dc11f962bc4ccb5d6c35f7827b9f13a`
- Ensemble repository: `/mnt/second19T/zhaozelin/velocity/transition_ensemble`
- Ensemble branch: `main`
- Ensemble commit: `de47402ea8c5625f765c595b279d8d53a4594cec`
- Ensemble worktree at audit: clean

The project root is not a Git repository. The ensemble subproject is the only audited component with a self-contained Git history.

## Working-tree scale

The following sizes explain why the public repository must be curated rather than copied wholesale:

| Source directory | Approximate size |
|---|---:|
| `data/` | 18 GB |
| `indicators/` | 45 GB |
| `methods/` | 441 GB |
| `plot/` | 2.8 GB |
| `results/` | 247 GB |
| `transition_ensemble/` | 7.3 GB |

## Ensemble result baseline

- Canonical result directory: `transition_ensemble/results/row_geometric_median_all14/`
- Final method: row-wise geometric median, referred to in the paper as `RowGeometricMedian-4` or `GeoMedian-4` according to context.
- Internal ablation: `Mean-4`.
- Dataset count: 14.
- Materialized result commit in ensemble history: `35071e4`.
- Current numerical verification source: `d0e7e5e`.
- Stable audited repository commit: `de47402`.
- Previously verified score comparison: 28 CBDir/BATC values matched the canonical table with maximum absolute difference `1.11e-16`.

## Phase 1 mutation boundary

This phase created only a new repository scaffold, manifests, and provenance documentation. It did **not**:

- modify the working project or paper;
- copy raw or processed datasets;
- copy method installations, checkpoints, or bulk results;
- rerun scientific analyses;
- declare a license, public URL, or archival data deposit.
