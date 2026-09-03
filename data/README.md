# Dataset metadata

`datasets.tsv` records the paper dataset identity, role, dimensions, stable
upstream source, and release status. It contains no server working paths.

Large canonical inputs are not stored in Git. Their external source or archive
location is recorded in the root [`artifacts.tsv`](../artifacts.tsv). A method
script may select a method-specific artifact, such as a `_counts.h5ad` input,
when its official workflow requires raw counts.
