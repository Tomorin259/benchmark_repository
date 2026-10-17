# Manuscript captions

## Supplementary Figure S15

Phase-assignment constraints in RNA-only MultiVelo.

(A) Gene-level velocity fields for a representative branching gene in the normal + branching simulation, showing ground truth (left), default MultiVelo (middle), and a matched refit allowing repression-phase candidates (Gate open; right). Both conditions retain the gene-level on classification; the Gate open refit disables the subsequent exclusion of repression-phase candidates during fitting. The annotated CosSim values are means over all 2,000 cells for this gene. All three plots share spliced--unspliced coordinates and representative arrow samples. Colors indicate known simulated progression-bin order from early to late. Arrow lengths are scaled separately for display.

(B) Paired gene-level CosSim for the 20 evaluated branching genes in the normal + branching (left) and hybrid (right) simulations. Filled and open points denote default and modified conditions, respectively; lines connect the same gene, and black bars and adjacent numbers indicate arithmetic means. For scVelo_dyn, the modified condition retains fitted rates and scaling but uses the induction trajectory and its corresponding smoothed times. For MultiVelo, it is the Gate open refit shown in panel A. The hybrid comparison includes branching genes only.

## Supplementary Table S9

Post-hoc controls of gene-level velocity inference in scVelo_dyn and MultiVelo.

Values in parts A and B are arithmetic means over the evaluated genes in each subset; $n$ denotes genes, not independent simulation replicates. Part A summarizes branching-gene controls, part B summarizes MURK-gene controls, and part C reports paired decreases after allowing repression candidates, expressed as the number of genes with a lower value divided by the number of paired genes. The scVelo induction-only condition retains fitted rates and scaling while using the induction trajectory and its corresponding smoothed times. MultiVelo Gate open retains the on classification and changes only repression-candidate exclusion during refitting. No output smoothing removes final velocity smoothing/rescaling. Complete trajectory also changes gap handling and is a supporting sensitivity control, not a single-factor ablation. Gate-open refits were evaluated only in branching subsets. All directional scores are displayed to four decimal places; full-precision results and gene-wise SDs are available in the repository.
