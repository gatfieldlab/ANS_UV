RNA-Seq Analysis

An R pipeline using DESeq2 and ggplot2 to perform quality control, differential expression analysis, and principal component analysis (PCA) across multiple experimental conditions (High Dose, Low Dose, UV treatment).

Pipeline Overview

Raw cDNA Count Table -> Filter Long Reads & Collapse Isoforms to Max Expression
QC Metrics -> Generate Sample Expression QC Plot
Differential Expression -> Run DESeq2 Analysis per Condition (HD, LD, UV)
Transformation -> Extract Variance Stabilizing Transformation (VST) & Generate Volcano Plots
Dimensionality Reduction -> Combined PCA on Top 500 Variable Genes

Installation & Requirements
Ensure R (>= 4.1.0) is installed along with the required CRAN and Bioconductor packages:# Install CRAN dependencies

install.packages(c("ggplot2", "data.table", "tidyverse", "dplyr", "tibble", "matrixStats", "ggrepel"))

# Install Bioconductor dependencies

if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("DESeq2")

Key Steps
Isoform Collapsing: Retains the single highest-expressed transcript isoform per gene_id based on total expression signal.
QC Metrics: Calculates total raw counts per condition and outputs a summary barplot.
Differential Expression: Constructs individual DESeq2 datasets (dds_hd, dds_ld, dds_uv) and evaluates fold changes against controls.
Volcano Plots: Highlights significant up- and down-regulated genes (padj < 0.05, |log2FC| > 1).
Combined PCA Analysis: Combines variance-stabilized expression matrices across experiments to compute PCA coordinates using the top 500 most variable genes.

Expected Outputs
QC Barplot
Total raw count distribution across sample conditions
Volcano Plots
PCA Figures
