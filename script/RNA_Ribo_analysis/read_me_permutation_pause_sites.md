Codon Usage + Permutation Testing

R pipeline for detecting codon usage enrichments and depletion at specific ribosome footprint/pausing sites.

Overview & Analysis Workflow

Fisher Test Input Files & Bed Annotations – Import and align raw site comparisons.

Filter Histones & Problematic Transcripts – Exclude multi-copy histones and known noise sources.

Annotate Relative Footprint Positions & Define Up Sites – Map coordinates to transcript regions.

Permutation Testing: 10,000 Label Shuffles – Build background distributions.

Compute empirical p-values & log2FC.

Generate plots – Compare % Codon up vs non-significant background.


Installation & Requirements
This pipeline requires R (>= 4.1.0) and the following CRAN / Bioconductor packages:# Install CRAN dependencies

install.packages(c(

  "tidyverse", "ggplot2", "stringr", 

  "ggthemes", "ggtext", "gridExtra", "ggrepel"

))

# Install Bioconductor dependencies

if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")

Pipeline Features

Gene universe restriction: Enforces strict gene-set intersection between significant and non-significant site comparisons to remove gene-expression bias.

Empirical permutation: Performs 10,000 label shuffles per dataset to establish expected background codon frequencies and calculate empirical two-sided p-values ($p_{up}$ / $p_{down}$).

Data Visualization:

Scatter plots: Displays observed codon frequency percentages in pause sites versus non-significant background sites.
Heatmap: Groups codons by amino acid and highlights significantly enriched or depleted codons across conditions with black dots.

Expected Outputs

Scatter Plots (output/figures/scatter_[condition].pdf): Per-condition scatter plots displaying codon frequency comparison between up-regulated and non-significant sites.
Heatmap (output/figures/heatmap_log2FC_all_conditions.pdf): Combined matrix heatmap displaying log2(Fold Change) values ordered by amino acid across monosome (mono) and disome (di) datasets.

##  License

This repository is shared for research and scientific data analysis.

