Ribo-seq Analysis Pipeline

R pipeline for processing, normalizing, filtering, and visualizing Ribosome Profiling (Ribo-seq) data. 

The script handles transcript quantification, histone gene family disambiguation, sequencing depth normalization (RPM), metagene profiling around translation START and STOP sites, and cumulative expression analysis.

Workflow Overview

QC & Filtering — Filter Histone Multi-Gene Families & Noisy Transcripts
Depth normalization — Merge replicates & RPM depth normalization
Absolute Metagene Profiling — Metagene Profiling: Raw Pileup & CDS Zoom Plots
Cumulative expression calculation

Installation & Dependencies
This pipeline requires R ($\ge$ 4.1.0) along with the following CRAN and Bioconductor packages:# Install CRAN dependencies

install.packages(c(

  "patchwork", "ggplot2", "cowplot", "dplyr", 

  "ggthemes", "readxl", "gridExtra", "tidyverse", 

  "data.table", "scales"

))

# Install Bioconductor dependencies

if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

BiocManager::install("rstatix")

Key Pipeline Steps

Most Expressed Transcript Quantification
Summarizes RSEM CDS raw counts to calculate condition means (nt vs. ans).
Histone Disambiguation & Noise Filtering
Identifies multi-gene histone families to retain only the most abundant representative gene per family, while excluding known high-noise transcripts (Lars2, Grcc10, Ftl1).
RPM Normalization
Aggregates positional A-site footprint counts per library (PD and PF series) and normalizes by Reads Per Million (RPM).
Metagene Profiling
Raw Pileup Metagene: Stacks raw normalized counts using library-tailored Y-axis scales, accompanied by zoomed CDS region views.
Cumulative Expression
Plots cumulative translation curves across CDS coordinates on a log10 scale.

Expected Outputs
normalized count table
Plots : metagene and cumulative plots

License
This repository is shared for research and scientific data analysis.
