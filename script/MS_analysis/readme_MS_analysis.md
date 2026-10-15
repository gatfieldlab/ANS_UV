Mass Spectrometry Analysis Pipeline for Monosome & Disome Profiling

This repository contains the R analysis pipeline used for processing, normalizing, and analyzing quantitative mass spectrometry (LC-MS/MS) dataset across Monosome, Disome, and Total Input fractions under various stress/treatment conditions (Low Dose, High Dose, UV).

---
This script requires R (>= 4.1.0) and the following standard and Bioconductor packages:

# CRAN Packages
install.packages(c(
  "tidyverse", "readxl", "ggrepel", "ggh4x", 
  "FactoMineR", "factoextra", "gridExtra"
))

# Bioconductor Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "clusterProfiler", "org.Mm.eg.db", "GO.db", "biomaRt"))

---

Input files : 

Place the following required files in the data/ directory:

Input & Ribo MS Tables:

Table8_LD_input.txt / Table8_LD_ribo.txt #input and ribo data and found in the MS table (see sheets)
Table8_HD_input.txt / Table8_HD_ribo.txt #input and ribo data and found in the MS table (see sheets)
Table8_UV_input.txt / Table8_UV_ribo.txt #input and ribo data and found in the MS table (see sheets)

Ribosomal Protein Reference List:

mgi_name_of_ribo_prot (MGI symbols for RPS and RPL ribosomal proteins)
Protein Highlight Lists (Optional for customized volcano/heatmaps):
list_prot_to_lab.xlsx  #file that contains the proteins that has to be labeled on the volcano plots
list_for_figure.txt #file that contains the proteins that has to be labeled on the volcano plots

---

Key Analysis Steps

Resolves ambiguous gene assignments by retaining the peptide group with the maximum total intensity signal.

Sample Mapping & Filtering: Renames technical sample codes into experimental conditions (e.g., Di_LD, mono_HDMSO) and filters out unassigned fractions.

Two-Step Normalization:
Input Normalization: Subtracts Total Input log2 intensities from Monosome/Disome fractions to account for overall protein abundance changes.Ribosomal Normalization: Normalizes signal to the median expression of core ribosomal proteins ($Rps$/$Rpl$) across samples.

Differential Analysis: Performs pairwise Student's t-tests and $log_2(\text{Fold Change})$ calculations across experimental conditions vs controls.

Enrichment & Visualization: Automated generation of GO term enrichment (clusterProfiler), heatmaps with significance indicators (ggh4x), and condition-specific Volcano plots (ggrepel).

---

Outputs : 

Differential Expression Tables: t_test_results.csv, t_test_results_withvalues.csv

QC Plots: Ribosomal normalization boxplots, overall expression violins, PCA plots per fraction.

Functional Figures: Volcano plots targeting Ribosome Quality Control (RQC) & Translation factors, GO term dotplots.

