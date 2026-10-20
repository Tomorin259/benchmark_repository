# Figure assembly

[assembly_manifest.tsv](assembly_manifest.tsv) maps paper panels to retained
sources, export destinations and checksums.
`source_path` identifies a retained source; `export_path` is its destination in
an on-demand material bundle.

| Status | Meaning |
| --- | --- |
| `referenced_source` | Retained experiment-level source |
| `retained_composite` | Retained final PDF with the recorded verification scope |
| `authoring_source` | Editable author source |

Figure 1 is a manually drawn schematic registered as a whole-figure PDF.

Use [the figure map](../../docs/figure_map.md) for figure-specific source
availability, layout and manual PDF-export steps.

```bash
python figures/export_materials.py --verify-only
python figures/export_materials.py --figures fig4 fig5 fig6 --output-root outputs/paper-materials
```

Run from the repository root with a fresh output directory.
The exporter checks source identity and copies materials without changing
canonical panels. Manual author layout and PDF export remain required where
listed in the figure map.
