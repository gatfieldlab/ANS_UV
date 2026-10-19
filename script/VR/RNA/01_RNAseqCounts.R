### LICENCE
# Copyright (c) 2026 Virginie Ricci

# This file is part of ANS_UV (https://github.com/gatfieldlab/ANS_UV/).

# ANS_UV includes free scripts: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, 
# either version 3 of the License, or (at your option) any later version.

# ANS_UV is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this script. 
# If not, see <https://www.gnu.org/licenses/>.
### LICENCE 

library(DEXSeq)
library(DESeq2)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(biomaRt)
library(stringr)
library(ggplot2)
library(tximport)
library(gridExtra)
library(ggrepel)
library(ggfortify)
library(ggpubr)
library(ggsignif)
library(cowplot)
library(ggh4x)
library(reshape2)
library(scales)
library(DESeq2)
library(Rsubread)
library(GenomicFeatures)
library(GenomicAlignments)
library(txdbmaker)
library(Rsamtools)
library(rtracklayer)
library(apeglm)
library(ashr)
library(ensembldb)
library(eisaR)
library(Repitools)
library(BiocParallel)
library(BSgenome)
library(ribosomeProfilingQC)
library(tibble)


path_home='/path_to/RNAseq/'
rsem_data=paste0(path_home, 'RSEM/')
rsem_res = '/path_to/RNAseq_Counts/'

system(paste0('mkdir -p ', rsem_res))

setwd(path_home)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
GENOME='Mmusculus.GRCm39.111'
DB='/data/databases/mouse/'
GENOME_DB=paste0(DB, 'star/', GENOME, '/')
DB_GTF=paste0(DB, 'gtf/')
DB_BED=paste0(DB, 'bed/')
DB_FASTA=paste0(DB, 'fasta/')
GTF_GENOME=paste0(DB_GTF, GENOME, '.sorted.gtf')
FASTA_CDNA=paste0(DB_FASTA, GENOME, '.cdna.ensembl.fa')

GTF = import(GTF_GENOME) # 1887003 ranges 
GTF
GTF_split = split(GTF, mcols(GTF)$gene_id)

tx_gene = unique(data.frame(mcols(GTF[GTF$type == 'transcript',][,c('gene_id', 'transcript_id', 'gene_name')])))
tx_gene$gene_tx = paste0(tx_gene$gene_id, '|', tx_gene$transcript_id)
head(tx_gene); dim(tx_gene)

GTF_gene = GTF[GTF$type == 'transcript', ]
GTF_gene
GTF_gene_biotype = unique(data.frame(ensembl_gene_id=GTF_gene$gene_id, ensembl_transcript_id=GTF_gene$transcript_id, gene_biotype=GTF_gene$gene_biotype))
head(GTF_gene_biotype); dim(GTF_gene_biotype)

ProteinCoding = GTF_gene[GTF_gene$gene_biotype == 'protein_coding',]
ProteinCoding
table(ProteinCoding$gene_biotype)
table(ProteinCoding$transcript_biotype)
GTF_gene_ProteinCoding = unique(data.frame(ensembl_gene_id = ProteinCoding$gene_id, mgi_symbol = ProteinCoding$gene_name))
head(GTF_gene_ProteinCoding); dim(GTF_gene_ProteinCoding)


###
CDS_len = read.csv(paste0(DB, 'Mmusculus.GRCm39.111.startStop.IDS.bed'), sep='\t', header=F)
names(CDS_len) = c('Gene_Tx', 'start', 'end', 'biotype')
CDS_len$gene_id = do.call(rbind, str_split(CDS_len$Gene_Tx, '[|]'))[,1]
CDS_len$transcript_id = do.call(rbind, str_split(CDS_len$Gene_Tx, '[|]'))[,2]
head(CDS_len); dim(CDS_len)
range(CDS_len$start)
range(CDS_len$end)


### cDNA fasta
cDNA = readDNAStringSet(FASTA_CDNA)
cDNA # 115911 sequences
head(names(cDNA))
cDNA_GR = GRanges(names(cDNA), IRanges(start=1, width=width(cDNA)), strand="*")
cDNA_GR # 115911
cDNA_GR$width = width(cDNA_GR)


GeneTx = str_split_fixed(names(cDNA), ' ', 2)[,1]
head(GeneTx); length(GeneTx)
length(GeneTx) == length(unique(GeneTx))

Genes = str_split_fixed(GeneTx, '[|]', 2)[,1]
head(Genes); length(Genes)
length(unique(Genes)) # 35939

Txs = str_split_fixed(GeneTx, '[|]', 2)[,2]
head(Txs); length(Txs)
length(unique(Txs)) # 115911



ScfTx_tmp = str_split_fixed(names(cDNA), ' ', 5)[,4]
ScfTx = str_split_fixed(ScfTx_tmp, ':', 6)[,3]
ScfTx_start = str_split_fixed(ScfTx_tmp, ':', 6)[,4]
ScfTx_end = str_split_fixed(ScfTx_tmp, ':', 6)[,5]
ScfTx_strand = str_split_fixed(ScfTx_tmp, ':', 6)[,6]
ScfTx_strand[ScfTx_strand == '1'] = '+'
ScfTx_strand[ScfTx_strand == '-1'] = '-'

Tx_tmp = str_split_fixed(names(cDNA), ' ', 5)[,5]
TxGene_tmp = str_split_fixed(Tx_tmp, ' ', 5)[,1]
TxGene_nb = str_split_fixed(TxGene_tmp, ':', 2)[,2]
TxGene = str_split_fixed(TxGene_nb, '[.]', 2)[,1]

TxGeneBiotype_tmp = str_split_fixed(Tx_tmp, ' ', 5)[,2]
TxGeneBiotype = str_split_fixed(TxGeneBiotype_tmp, ':', 2)[,2]

TxBiotype_tmp = str_split_fixed(Tx_tmp, ' ', 5)[,3]
TxBiotype = str_split_fixed(TxBiotype_tmp, ':', 2)[,2]

TxGeneSymbol_tmp = str_split_fixed(Tx_tmp, ' ', 5)[,4]
TxGeneSymbol = str_split_fixed(TxGeneSymbol_tmp, ':', 2)[,2]

TxDescr_tmp = str_split_fixed(Tx_tmp, ' ', 5)[,5]
TxDescr = str_split_fixed(TxDescr_tmp, ':', 2)[,2]

mcols(cDNA_GR)$gene_id = Genes
mcols(cDNA_GR)$transcript_id = Txs
seqlevels(cDNA_GR) = cDNA_GR$transcript_id
seqnames(cDNA_GR)
cDNA_GR
cDNA_df = data.frame(cDNA_GR)
cDNA_df$width.1 = NULL
cDNA_df$width = NULL
names(cDNA_df)[names(cDNA_df) == 'start'] = 'cDNA_start'
names(cDNA_df)[names(cDNA_df) == 'end'] = 'cDNA_end'
names(cDNA_df)[names(cDNA_df) == 'width'] = 'cDNA_width'
head(cDNA_df); dim(cDNA_df)

head(CDS_len); dim(CDS_len)
CDS_len_cDNA = left_join(CDS_len, cDNA_df)
head(CDS_len_cDNA); dim(CDS_len_cDNA)

CDS_len_cDNA_ProtCoding = CDS_len_cDNA[CDS_len_cDNA$biotype == 'protein_coding|protein_coding',]
CDS_len_cDNA_ProtCoding = CDS_len_cDNA_ProtCoding[CDS_len_cDNA_ProtCoding$start != 0,]
head(CDS_len_cDNA_ProtCoding); dim(CDS_len_cDNA_ProtCoding) # 50691
length(unique(CDS_len_cDNA_ProtCoding$gene_id)) # 20790
###

### CDS start and stop coordinates
GR_CDS = GRanges(CDS_len_cDNA_ProtCoding$transcript_id, IRanges(start=CDS_len_cDNA_ProtCoding$start, end=CDS_len_cDNA_ProtCoding$end))
names(GR_CDS) = seqnames(GR_CDS)
GR_CDS
GR_CDS_info = left_join(data.frame(transcript_id=names(GR_CDS)), tx_gene)
head(GR_CDS_info); dim(GR_CDS_info)

GR_CDS$gene_id = GR_CDS_info$gene_id
GR_CDS$transcript_id = GR_CDS_info$transcript_id
GR_CDS$gene_tx = paste0(GR_CDS_info$gene_id, '|', GR_CDS_info$transcript_id)
GR_CDS

### Full cDNA
Full_cDNA = GRanges(CDS_len_cDNA_ProtCoding$transcript_id, IRanges(start=CDS_len_cDNA_ProtCoding$cDNA_start, end=CDS_len_cDNA_ProtCoding$cDNA_end))
names(Full_cDNA) = seqnames(Full_cDNA)
Full_cDNA
Full_cDNA_info = left_join(data.frame(transcript_id=names(Full_cDNA)), tx_gene)
head(Full_cDNA_info); dim(Full_cDNA_info)

Full_cDNA$gene_id = Full_cDNA_info$gene_id
Full_cDNA$transcript_id = Full_cDNA_info$transcript_id
Full_cDNA$gene_tx = paste0(Full_cDNA_info$gene_id, '|', Full_cDNA_info$transcript_id)
Full_cDNA



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# .bam files - PROJECTION ON TRANSCRIPTOME OF THE GENOME MAPPING with RSEM
BAM=list.files(path=rsem_data, pattern='*/*.transcript.sorted.bam$', full.names=T, recursive=TRUE)
BAM = str_replace(BAM, '//', '/')
BAM_split = do.call(rbind,str_split(BAM, '/'))
Sample = BAM_split[,7]
Sample # e.g. PD63_AGCTA
BAM # e.g. /path_to/ANS_HD_UV_LD/RSEM/PD63/AGCTA/PD63_AGCTA.transcript.sorted.bam
BAM_df = data.frame(Sample = Sample, bam = BAM)
BAM_df

Bam_Files = BamFileList(BAM_df$bam)
names(Bam_Files) = BAM_df$Sample


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
### HTSeq-count in CDS of protein-coding transcripts

### CDS
Counts_CDS <- GenomicAlignments::summarizeOverlaps(GR_CDS, Bam_Files, mode='Union')
colnames(Counts_CDS) = BAM_df$Sample
colData(Counts_CDS)$Sample = BAM_df$Sample
colSums(assay(Counts_CDS))

Counts_CDS_df = data.frame(assay(Counts_CDS))
Counts_CDS_df = rownames_to_column(Counts_CDS_df, "transcript_id")
head(Counts_CDS_df); dim(Counts_CDS_df)



### cDNA
Counts_cDNA <- GenomicAlignments::summarizeOverlaps(Full_cDNA, Bam_Files, mode='Union')
colnames(Counts_cDNA) = BAM_df$Sample
colData(Counts_cDNA)$Sample = BAM_df$Sample
colSums(assay(Counts_cDNA))

Counts_cDNA_df = data.frame(assay(Counts_cDNA))
Counts_cDNA_df = rownames_to_column(Counts_cDNA_df, "transcript_id")
head(Counts_cDNA_df); dim(Counts_cDNA_df)



Counts_CDS_df = left_join(Counts_CDS_df, tx_gene[c('gene_id', 'transcript_id')])
Counts_CDS_df = Counts_CDS_df[c(ncol(Counts_CDS_df), 1:(ncol(Counts_CDS_df)-1))]
head(Counts_CDS_df)

Counts_cDNA_df = left_join(Counts_cDNA_df, tx_gene[c('gene_id', 'transcript_id')])
Counts_cDNA_df = Counts_cDNA_df[c(ncol(Counts_cDNA_df), 1:(ncol(Counts_cDNA_df)-1))]
head(Counts_cDNA_df)

write.table(Counts_CDS_df, file = paste0(rsem_res, 'Counts_CDS.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
write.table(Counts_cDNA_df, file = paste0(rsem_res, 'Counts_cDNA.txt'), quote = F, sep = '\t', row.names = F, col.names = T)


