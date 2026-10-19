## Anisomycin-induced ribosome collisions follow a dose-dependent positional rule, distinct from the UV-induced collision landscape
(https://www.biorxiv.org/content/10.64898/2026.10.12.758323v1)

Sébastien Chamois, Julien Soudet, Virginie Ricci, Félix Naef, Cédric Gobet, David Gatfield

### <ins> Summary:</ins>
This project investigates how distinct ribosome-impairing stimuli (Anisomycin and UV irradiation) shape the distribution and molecular composition of ribosome collisions across mRNAs in Mouse Embryonic Stem Cells (mESCs).

### <ins> Experimental methods & data:</ins>
Transcriptome analysis and ribosome profiling: RNAseq, Monosome & Disome Profiling (Ribo-seq / Disome-seq)
Proteomics: Quantitative mass spectrometry of isolated collided-ribosome fractions
Model System: Mouse Embryonic Stem Cells (mESCs)

### <ins> Scripts:</ins>

### RNA-seq
[pipeline](https://github.com/gatfieldlab/Snakemake_pipeline): 
- Reads trimming and mapping to mouse genome

01_RNASeqCounts.R: 
- Read count in CDS and cDNA

### Ribo-seq
[pipeline](https://github.com/gatfieldlab/Snakemake_pipeline):
- Reads trimming and filtering
- Sequential mapping to mouse rRNA, human rRNA, and mouse tRNA
- Demultiplexing (barcodes) and mapping to mouse genome
- Deduplication of reads mapped to mouse genome

01_RiboSeqCounts_SubBam.R:
- Read counts in CDS
- Get the most expressed transcript per gene
- filter BAM

01_Ribo_Asite_Monosome.R:
- Get A-site position in Monosome (offset estimated using RiboWaltz)

02_Ribo_Asite_Disome.R:
- Get A-site position in Disome data (offset: +45)

04_Ribo_CountsInAA.R:
- For each transcript, get counts per amino acid

04_Ribo_FisherTest.R:
- For each transcript, perform FisherTest for each amino acid position




