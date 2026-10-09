# Dataset metadata

`datasets.tsv` records the paper dataset identity, role, dimensions, stable
upstream source, and release status. It contains no server working paths.

Large canonical inputs are not stored in Git. The root
[`artifacts.tsv`](../artifacts.tsv) records their logical identities, expected
checksums, and intended external-location fields. Stable public URIs have not
yet been published, so the registry is not currently a download index. A
method script may select a method-specific artifact, such as a `_counts.h5ad`
input, when its official workflow requires raw counts.

## External release preparation

`release_map.tsv` covers every external artifact ID exactly once. It adds only
package/member layout and logical source locators; hashes and publication state
remain owned by `artifacts.tsv`. Locators recovered from frozen Cell, kinetics,
veloVI and simulation manifests are not treated as current-byte verification.
Unfrozen input identity and unrecovered locators remain explicit.

Check one package without creating a large archive:

```bash
python data/prepare_release.py --artifact-root /path/to/velocity-artifacts \
  --package reference-moments --output-root outputs/reference-release-check
```

Add `--build` to write a `.tar.gz` only after every selected source matches its
registered SHA-256. The archive preserves the logical velocity-root paths used
by the existing launchers/evaluators. Missing sources or changed bytes block
archive construction; no new canonical input or score is inferred. The script
never uploads or fills download URIs. Record real URIs only after an actual
publication in the selected artifact store. Public storage selection is pending.
