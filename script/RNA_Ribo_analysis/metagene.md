# Ribo-seq Data Analysis Pipeline

A comprehensive R pipeline for processing, filtering, normalising, and plotting Ribosome Profiling (Ribo-seq) data. This script manages transcript quantification, multi-gene histone family filtering, sequencing depth normalisation (RPM), metagene profiling around START/STOP sites, and cumulative translation analysis.

---

## 🛠️ Overview & Workflow

```mermaid
flowchart TD
    A[Raw Counts & SGA Files] --> B[Load Reference & APPRIS Annotations]
    B --> C[Select Top 20% Expressed Transcripts via RSEM]
    C --> D[Filter Histone Multi-Gene Families]
    D --> E[Merge Replicates & Depth Normalisation RPM]
    E --> F1[Metagene Analysis: Relative Counts per Gene]
    E --> F2[Metagene Analysis: Raw Pileup & Zoom Plots]
    E --> F3[Cumulative Translation Plots]