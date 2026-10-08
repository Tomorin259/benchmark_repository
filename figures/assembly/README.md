# Paper figure assembly sources

`assembly_manifest.tsv` is the authoritative mapping from each paper panel to
one retained source file and its SHA-256. `export_path` is a destination inside
an on-demand PowerPoint/WPS material bundle, not a second tracked copy.
`referenced_source` denotes a stored source; `retained_composite` denotes an
online-matched final PDF; `authoring_source` denotes a native manual source.
`missing_source` remains explicit when independent final panel extraction is
unavailable. Retaining a PDF/PPTX is not proof of a numerical renderer.

Verify sources without writing or exporting:

```bash
python figures/export_materials.py --verify-only
```

Create a fresh material bundle, preserving every file byte:

```bash
python figures/export_materials.py --figures fig4 fig5 fig6 --output-root outputs/paper-materials
```

The bundle records source hashes and source gaps. Its directory must be empty;
canonical sources are never overwritten. Module panels remain in their owning
evaluation/diagnostic directory. Native PPTX/manual sources are stored here.

Figure 1 retains the original authoring PPTX and exact online PDF. The user
intentionally keeps the three original PPT labels; it is not described as an
exact independently re-exported source of the corrected online PDF.
Figure 4 retains the native schematic/assembly source and verified author-export
lineage. Figure 5D's two overview images were recovered from native slide
placement; its previously misidentified legacy gene-summary bitmap is archived.
The validated Figure 6 source panels and author PDF retain their existing
scientific inputs and hashes. Figure 7 retains the editable native charts and
embedded workbooks, with manual WPS PDF export.

Supplementary Figures 9/11 retain exact online whole-PDF bytes; independent
editable/numerical whole-figure reconstruction remains a separate boundary.
Supplementary Figures 12/13 retain their approved composites and provenance.
