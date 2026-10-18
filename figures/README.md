# Paper figures

Use the [figure map](../docs/figure_map.md) to find each paper figure,
its numerical sources, renderers and editable assembly.

Final composites are retained here; experiment-level panels stay with their
producing modules. [assembly/](assembly/README.md) contains editable author
sources and the panel-to-source map.

## Reproduce materials

Run from the repository root:

```bash
python figures/export_materials.py --verify-only
python figures/export_materials.py --figures fig4 fig5 fig6 --output-root outputs/paper-materials
```

The output directory must be empty. Export copies the registered source files
for the assembly steps listed in the figure map. Composite manifests record
the numerical inputs, author exports and any
reconstruction gaps. Historical kinetics artwork under `archive/kinetics/`
uses the former Figure 5 numbering; current Figure 5 is gene-level evaluation.
