Heatmap, figure 3 script

An R pipeline combining Ribo-seq footprint counts and matched RNA-seq data to quantify Translation Efficiency (TE) via DESeq2, map positional ribosome pause sites from Fisher test results, and generate coordinated heatmaps of expression and stalling density.

Analysis Workflow

Load cDNA & Footprint Count Matrices
DESeq2 normalization (RNA + Ribo) -> calculate Translation Efficiency
Fishertest pause site processing -> Map positions relative to Start Codon
Prepare logarithmic violin plots of pause sites (Up vs Down Stalls)
Paired Expression and Pause-Count Heatmaps Ordered by Gene Expression

Dependencies & Installation

Requires R (>= 4.1.0) and the following CRAN and Bioconductor packages:# Install CRAN packages

install.packages(c(

  "patchwork", "ggplot2", "cowplot", "dplyr", 

  "ggthemes", "readxl", "gridExtra", "tidyverse", 

  "data.table", "pheatmap", "scales"

))

# Install Bioconductor packages

if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "rstatix"))

Key Pipeline Features
DESeq2 : Constructs a multi-factor model (~ compare + assay + compare:assay) combining RNA-seq and Ribo-seq datasets across High Dose (HD), Low Dose (LD), and UV treatment conditions.

Translation Efficiency (TE) Calculation: Computes gene-level TE.

Logarithmic pause site mapping: Evaluates ribosome stalling positions relative to translation START codons and generates log10-scaled violin density plots.

Coordinated Heatmaps: Constructs paired expression and pause-site density heatmaps sorted by decreasing expression.

Expected Outputs
Violin Plots
Expression Heatmaps
Pause Count Heatmaps
