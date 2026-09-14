# Dataset metadata

`datasets.tsv` records the paper dataset identity, role, dimensions, stable
upstream source, and release status. It contains no server working paths.

Large canonical inputs are not stored in Git. The root
[`artifacts.tsv`](../artifacts.tsv) records their logical identities, expected
checksums, and intended external-location fields. Stable public URIs have not
yet been published, so the registry is not currently a download index. A
method script may select a method-specific artifact, such as a `_counts.h5ad`
input, when its official workflow requires raw counts.
