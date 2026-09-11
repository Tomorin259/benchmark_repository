# Kinetics plotting

This directory owns the plot-generation code for the kinetics-recovery
experiment. Figure 5 quantitative panels retain the historical plotting
workflow and use the explicitly mapped legacy CSV inputs; the DeepKINET gamma
CSV inputs point to the final seed-123 result.

The visible method labels use the final display names, including
`DeepVelo` and the kinetics-specific `Dynamo(label data)`. The global
20-method configuration retains `Dynamo(RNA_seq)`; the internal `rpe1` and `organoid`
configuration keys and output filenames are preserved; their manuscript-facing
dataset names are `rpe1` and `organoid`.

## Entries

- figure5.py generates the `rpe1` trajectory panels, `rpe1` heatmap,
  kinetic-class overlap panels, and `organoid` branch boxplots using the legacy layout and
  plotting parameters.
- figure5_config.json maps the historical method-result CSVs and the
  alpha/gamma trajectories used to recompute both kinetic-class Venn panels.
  The plotter writes count tables and per-gene membership tables under
  `plot_data/`.
- streams.py generates the Figure 5E `organoid` stream panels and the current
  DeepKINET `rpe1` stream used by Figure 5A using configured method-specific
  velocity-layer keys.
- stream_config.json records the organoid H5AD inputs, the current DeepKINET
  RPE1 result/reference pair, and the four remaining manually selected RPE1
  stream assets.

The final multi-panel Figure 5 and the integrated Figure 6 are still assembled
in the paper project. Their individual panel inputs and generators are
recorded in the figure assembly manifest under figures/assembly/.

Examples:

    python evaluation/kinetics_recovery/plotting/figure5.py --artifact-root /path/to/kinetics-artifacts --output /path/to/figures
    python evaluation/kinetics_recovery/plotting/streams.py --artifact-root /path/to/velocity --panel organoid
    python evaluation/kinetics_recovery/plotting/streams.py --artifact-root /path/to/velocity --panel rpe1_deepkinet
