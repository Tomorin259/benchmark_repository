# Supplementary Figure 11 — veloAE representation diagnostics

The online `merged.tex` defines two panels for Supplementary Figure 11:

**(A)** PC1 explained variance of the state matrix (`Ms`, left) and velocity
(right) for veloAE, scVelo stochastic, and TIvelo across the 14 real datasets.
Filled and open points denote raw PCA and PCA after per-dimension
standardization; lines connect observations from the same dataset, and black
bars denote cross-dataset means. PCA is fitted independently to each state and
velocity matrix.

**(B)** veloAE state coordinates and velocities projected onto the same
PC1–PC2 axes fitted to raw `Ms` for each dataset. State coordinates are centered,
whereas velocity vectors are projected without centering. Cells are colored by
dataset annotations, and arrows show spatial-bin averages of projected
velocities. Axis ranges are dataset-specific, and arrows are scaled and
length-capped for display; arrow lengths must not be compared across datasets.

The repository keeps the A and B source assets and their provenance. The final
two-panel composite uses manual author assembly.
