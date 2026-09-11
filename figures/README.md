# Paper figures

This directory contains the confirmed paper-level composite PDFs and the
self-contained paper panel assembly. Every available panel in `assembly/` is a
regular file whose bytes are recorded by SHA256 in
`assembly/assembly_manifest.tsv`; canonical source paths are repository-relative
and no assembly panel points outside the repository.

The manifest covers the current paper panels for Figures 2–6 and SuppFig.01–09.
It explicitly records four unresolved provenance/availability items: Fig.4D
canonical source is no longer present, and independent panel files are missing
for Fig.3A, Fig.6 whole, and SuppFig.09. No substitute was guessed. The
preserved top-level PDFs are `figTwo.pdf`, `figFive.pdf`, and `figSix.pdf`.

Do not add intermediate panels, alternate result copies, or old
