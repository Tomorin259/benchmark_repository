# Obtain the original kinetic-reference tables

The input distribution provides 25 mirrored scientific files. Obtain the three
additional kinetic-reference dependencies, Battich Tables S1/S3/S4, from the
original study using the steps below.

1. Visit [Battich et al. (2020)](https://doi.org/10.1126/science.aax3072) and open Supplementary Materials.
2. Obtain Tables S1, S3 and S4 through the original publisher/author access route.
3. Place the CSV files listed below at the paths in
   [upstream_references.tsv](upstream_references.tsv), using the layout required
   by your evaluator.
4. From this repository root, run `python data/check_references.py --artifact-root /path/to/artifacts`.
   In the HF version folder, run `python check_references.py --artifact-root /path/to/artifacts` instead.

| Table | Expected CSV | Contents |
| --- | --- | --- |
| S1 | `aax3072_table-s1.csv` | RPE1 regulatory strategies and normalized kinetic trajectories |
| S3 | `aax3072_table-s3.csv` | Organoid differential-expression results and the `selected` gene flag |
| S4 | `aax3072_table-s4.csv` | Organoid regulatory strategies and branch-specific normalized kinetic trajectories |

There are two documented layouts: velocity-root paths used by existing
launchers and `references/kinetics/` aliases used by the kinetics evaluator.
If both copies exist, both must match the recorded benchmark CSV identity.

## Expected format and checksums

The benchmark reads CSVs with a descriptive first line and a column-header
second line (`skiprows=1` in the evaluator). A CSV matching the recorded size
and SHA-256 can be used directly without format conversion. Preserve its
headers, gene identifiers, time-column order, numerical values and missing
entries.

`upstream_references.tsv` records their sizes and SHA-256 values.
`check_references.py` checks those file identities without rewriting the data.
For an Excel source or a different CSV serialization, compare its column
layout, gene identities, reference trajectories and missing values with the
benchmark CSV before evaluation. A different file hash can reflect serialization
differences; it does not by itself establish different table values.

Original citation: Battich et al., *Sequencing metabolically labeled
transcripts in single cells reveals mRNA turnover strategies*, Science 367,
1151–1156 (2020), doi:10.1126/science.aax3072.
