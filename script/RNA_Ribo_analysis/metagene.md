Markdown# Ribo-seq Data Analysis Pipeline

A comprehensive R pipeline for processing, filtering, normalising, and plotting Ribosome Profiling (Ribo-seq) data. This script manages transcript quantification, multi-gene histone family filtering, sequencing depth normalisation (RPM), metagene profiling around START/STOP sites, and cumulative translation analysis.

---

## 🛠️ Overview & Workflow

```mermaid
flowchart TD
    A[Raw Footprint SGA Files & References] --> B[Load Annotations & CDS Counts]
    B --> C[Filter Histone Multi-Gene Families & Blacklisted Transcripts]
    C --> D[Merge Replicates & Depth Normalisation RPM]
    D --> E1[Metagene Analysis: Relative Counts per Gene]
    D --> E2[Metagene Analysis: Raw Pileup & Zoom Plots]
    D --> E3[Cumulative Translation Plots]
📁 Recommended Repository StructureTo ensure the script runs smoothly using generic relative paths, organize your directory as follows:Plaintext├── annotations/                # Reference annotation files (BED, Gene Lengths, Biomart)
│   ├── mart_export.txt
│   ├── gene_transcript_list.csv
│   ├── Mmusculus.GRCm39.111.startStop.APPRIS.bed
│   ├── Mmusculus.GRCm39.111.startStop.IDS.bed
│   ├── Mmusculus.GRCm39.111.GeneLength.txt
│   ├── Histones_VR.txt
│   └── LengthsExonsIntronsGenes.txt
├── data/                       # Input data folders
│   ├── rsem/                   # RSEM outputs (MostExpressed & SummarizeOverlaps_CDS_raw)
│   └── sga/                    # SGA A-site footprint count files
├── output/                     # Generated results (auto-created by script)
│   ├── tables/                 # Merged counts, RPM tables, and histone filtering reports
│   └── figures/                # Output PDF plots (Metagene, Pileups, Cumulative curves)
│       └── primary_alignment/
├── riboseq_pipeline.R          # Main R analysis script
└── README.md                   # Documentation
⚙️ Dependencies & InstallationThis pipeline requires R (>= 4.1.0) and the following CRAN and Bioconductor packages:R# Install CRAN packages
install.packages(c(
  "patchwork", "ggplot2", "cowplot", "dplyr", 
  "ggthemes", "readxl", "gridExtra", "tidyverse", 
  "data.table", "scales"
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("rstatix")
🚀 Key Pipeline StepsReference & Annotation Setup: Automatically creates required output directories and loads genome-wide transcript maps, APPRIS annotations, exon lengths, and histone lists.CDS Quantification Summary: Processes RSEM CDS raw counts for sample condition means and identifies expression quantiles.Histone Family Disambiguation: Identifies multi-gene histone families, selects the single most abundant representative gene per family, and strips out secondary histone copies along with blacklisted transcripts (Lars2, Grcc10, Ftl1) to prevent mapping artifacts.RPM Normalisation: Aggregates positional footprint counts per library (PD and PF series) and normalises to Reads Per Million (RPM).Metagene Profiling:Gene-Normalised Metagene: Computes relative footprint density per gene across relative nucleotide ($nt$) positions from START and STOP sites (-75 to +200 nt).Raw Pileup Metagene: Aggregates raw footprint density with condition-tailored Y-axis limits (PF vs. PD series) and zoomed CDS views (-1000 to +6000 nt).Cumulative Expression Analysis: Generates cumulative translation curves across CDS coordinates on a $\log_{10}$ scale.📊 Expected OutputsTables (output/tables/):[lib]_merge_gene_sga_DESeq2.csv: Positional count matrices per library.[lib]_Asite_merge_gene_norm_depth.csv: RPM-normalised positional data with relative START/STOP positions.final_histone_list_to_keep.csv / histones_to_remove.csv: Histone filtering logs.[lib]_meta_merge_final_nohistone_start.csv: Relative metagene positional averages.Figures (output/figures/):[lib]_start_stop_metagene.pdf: Gene-normalised START/STOP metagene plots.[lib]_start_stop_pileup_metagene.pdf: Stacked raw pileup metagene panels.[lib]_start_stop_pileup_metagene_zoom.pdf: CDS region zoom plots.[lib]_Asite_merge_gene_norm_depth_cumulplot.pdf: Cumulative footprint density plots.📝 UsageClone this repository to your local machine.Place your raw SGA and RSEM files into data/sga/ and data/rsem/.Place your reference annotations in annotations/.Open riboseq_pipeline.R in RStudio or execute from command line:BashRscript riboseq_pipeline.R
📜 LicenseThis project is open for scientific research and internal laboratory pipelines.