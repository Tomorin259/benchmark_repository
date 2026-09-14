# Paper figures

This directory contains the confirmed paper-level composite PDFs and the
self-contained paper panel assembly. Every available panel in `assembly/` is a
regular file whose bytes are recorded by SHA256 in
`assembly/assembly_manifest.tsv`; canonical source paths are repository-relative
and no assembly panel points outside the repository.

The manifest covers the current paper panels for Figures 1–6 and SuppFig.01–09.
It explicitly records five unresolved panel-level provenance/availability
items: independent panel files are missing for Fig.1A, Fig.1B, Fig.3A, and
Fig.6 whole, while the canonical source for Fig.4D is no longer present. No
substitute was guessed. The preserved top-level PDFs are `figTwo.pdf`,
`figFive.pdf`, and `figSix.pdf`.

Do not add intermediate panels, alternate result copies, or obsolete
historical composites to this directory.
