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


path_home='/path_to/ANS_HD_UV_LD/'
rsem_data=paste0(path_home, 'RSEM/')
rsem_res = paste0(rsem_data, 'sub_bam/')

system(paste0('mkdir -p ', rsem_res))

setwd(path_home)


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
### Reference genome
GENOME='Mmusculus.GRCm39.111'
DB='/data/databases/mouse/'
DB_GTF=paste0(DB, 'gtf/')
GTF_GENOME=paste0(DB_GTF, GENOME, '.sorted.gtf')
IDS_GENOME=paste0(DB, GENOME, '.startStop.IDS.bed')

##### ##### ##### 
# GTF - annotation
GTF = import(GTF_GENOME) # 1887003 ranges 
GTF

GTF_gene = GTF[GTF$type == 'transcript', ]

table(GTF_gene$gene_biotype)
table(GTF_gene$transcript_biotype)

tx_gene = unique(data.frame(mcols(GTF_gene[,c('gene_id', 'transcript_id', 'gene_name', 'transcript_name', 'gene_biotype', 'transcript_biotype')])))
tx_gene$gene_tx = paste0(tx_gene$gene_id, '|', tx_gene$transcript_id)
head(tx_gene); dim(tx_gene) # 149076

# Protein-coding transcripts
ProteinCoding = tx_gene[tx_gene$transcript_biotype == 'protein_coding',]$transcript_id
head(ProteinCoding); length(ProteinCoding) # 58405

tx_gene_ProteinCoding = tx_gene[tx_gene$transcript_id %in% ProteinCoding, ]
head(tx_gene_ProteinCoding); dim(tx_gene_ProteinCoding) # 58405
length(unique(tx_gene_ProteinCoding$gene_id)) # 21469
unique(tx_gene_ProteinCoding$transcript_biotype)

# canonical transcripts
Tx_canonical = GTF_gene[grep('canonical', GTF_gene$tag),]$transcript_id
head(Tx_canonical); length(Tx_canonical) # 57126

# Protein-coding transcripts that are canonical
tx_gene_ProteinCoding_canonical = tx_gene_ProteinCoding[tx_gene_ProteinCoding$transcript_id %in% Tx_canonical,]
head(tx_gene_ProteinCoding_canonical); dim(tx_gene_ProteinCoding_canonical) # 21469


##### ##### ##### ##### ##### ##### 
### CDS start and stop positions
CDS_len = read.csv(IDS_GENOME, sep='\t', header=F)
names(CDS_len) = c('Gene_Tx', 'start', 'end', 'biotype')
CDS_len$gene_id = do.call(rbind, str_split(CDS_len$Gene_Tx, '[|]'))[,1]
CDS_len$transcript_id = do.call(rbind, str_split(CDS_len$Gene_Tx, '[|]'))[,2]
head(CDS_len); dim(CDS_len)
range(CDS_len$start)
range(CDS_len$end)

# CDS of protein-coding transcripts
CDS_len_ProteinCoding = CDS_len[CDS_len$transcript_id %in% tx_gene_ProteinCoding$transcript_id,]
head(CDS_len_ProteinCoding); dim(CDS_len_ProteinCoding) # 58405
range(CDS_len_ProteinCoding$start)
range(CDS_len_ProteinCoding$end)
unique(CDS_len_ProteinCoding$biotype)

dim(CDS_len_ProteinCoding[CDS_len_ProteinCoding$start == 0,]) # 7714
dim(CDS_len_ProteinCoding[CDS_len_ProteinCoding$end > 1000000000,]) # 4774

CDS_len_ProteinCoding_filtered = CDS_len_ProteinCoding[CDS_len_ProteinCoding$start !=0, ]
head(CDS_len_ProteinCoding_filtered); dim(CDS_len_ProteinCoding_filtered) # 50691

CDS_len_ProteinCoding_filtered_final = CDS_len_ProteinCoding_filtered[CDS_len_ProteinCoding_filtered$end < 1000000000,]
head(CDS_len_ProteinCoding_filtered_final); dim(CDS_len_ProteinCoding_filtered_final) # 50691

range(CDS_len_ProteinCoding_filtered_final$start)
range(CDS_len_ProteinCoding_filtered_final$end)

# CDS of protein-coding transcripts that are canonical (start=0 and end=10^6 removed)
CDS_len_ProteinCoding_canonical = CDS_len_ProteinCoding_filtered_final[CDS_len_ProteinCoding_filtered_final$transcript_id %in% tx_gene_ProteinCoding_canonical$transcript_id,]
head(CDS_len_ProteinCoding_canonical); dim(CDS_len_ProteinCoding_canonical) # 20265
range(CDS_len_ProteinCoding_canonical$start)
range(CDS_len_ProteinCoding_canonical$end)


##### ##### ##### ##### ##### ##### 
# GRANGES of CDS start and stop coordinates
GR_CDS = GRanges(CDS_len_ProteinCoding_filtered_final$transcript_id, IRanges(start=CDS_len_ProteinCoding_filtered_final$start+1, end=CDS_len_ProteinCoding_filtered_final$end+1))
names(GR_CDS) = seqnames(GR_CDS)
GR_CDS
GR_CDS_info = left_join(data.frame(transcript_id=names(GR_CDS)), tx_gene)
head(GR_CDS_info); dim(GR_CDS_info)

GR_CDS$gene_id = GR_CDS_info$gene_id
GR_CDS$transcript_id = GR_CDS_info$transcript_id
GR_CDS$gene_tx = paste0(GR_CDS_info$gene_id, '|', GR_CDS_info$transcript_id)
GR_CDS

range(start(GR_CDS))
range(end(GR_CDS))


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# .bam files - PROJECTION ON TRANSCRIPTOME OF THE GENOME MAPPING with RSEM
BAM=list.files(path=rsem_data, pattern='*/*.transcript.sorted.bam$', full.names=T, recursive=TRUE)
BAM = str_replace(BAM, '//', '/')
BAM_split = do.call(rbind,str_split(BAM, '/'))
Sample = paste0(BAM_split[,6], '_', BAM_split[,7])
Sample # e.g. PD63_AGCTA
BAM # e.g. /path_to/ANS_HD_UV_LD/RSEM/PD63/AGCTA/PD63_AGCTA.transcript.sorted.bam
BAM_df = data.frame(Sample = Sample, bam = BAM)
BAM_df

Bam_Files = BamFileList(BAM_df$bam)
names(Bam_Files) = BAM_df$Sample


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
### HTSeq-count in CDS of protein-coding transcripts
Counts_CDS <- GenomicAlignments::summarizeOverlaps(GR_CDS, Bam_Files, mode='IntersectionNotEmpty')
colnames(Counts_CDS) = BAM_df$Sample
colData(Counts_CDS)$Sample = BAM_df$Sample
colSums(assay(Counts_CDS))

Counts_CDS_df = data.frame(assay(Counts_CDS))
Counts_CDS_df = rownames_to_column(Counts_CDS_df, "transcript_id")
head(Counts_CDS_df); dim(Counts_CDS_df)



########################
### select the most expressed transcript per gene
Counts_CDS_df_copy = Counts_CDS_df
head(Counts_CDS_df_copy); dim(Counts_CDS_df_copy) # 50691

Counts_CDS_df_copy = left_join(Counts_CDS_df_copy, tx_gene[c('gene_id', 'transcript_id')])
head(Counts_CDS_df_copy); dim(Counts_CDS_df_copy)
Counts_CDS_df_copy = Counts_CDS_df_copy[c(ncol(Counts_CDS_df_copy), 1:(ncol(Counts_CDS_df_copy)-1))]
head(Counts_CDS_df_copy); dim(Counts_CDS_df_copy)

Counts_CDS_df_copy$geneMeans = rowMeans(Counts_CDS_df_copy[,-c(1,2)])
Counts_CDS_df_copy = Counts_CDS_df_copy[Counts_CDS_df_copy[,-c(1,2)]$geneMeans > 0,]
head(Counts_CDS_df_copy); dim(Counts_CDS_df_copy) # 42801

write.table(Counts_CDS_df_copy, file = paste0(rsem_res, 'SummarizeOverlaps_CDS_raw.txt'), quote = F, sep = '\t', row.names = F, col.names = T)

# most expressed tx per gene
Counts_CDS_df_copy_melt = reshape2::melt(Counts_CDS_df_copy[c('gene_id', 'transcript_id', 'geneMeans')])
head(Counts_CDS_df_copy_melt); dim(Counts_CDS_df_copy_melt)
Counts_CDS_df_copy_melt_filter = Counts_CDS_df_copy_melt %>% dplyr::group_by(gene_id) %>% dplyr::slice(which.max(value))
head(Counts_CDS_df_copy_melt_filter); dim(Counts_CDS_df_copy_melt_filter)
length(unique(Counts_CDS_df_copy_melt_filter$gene_id)) # 16955
selected_tx = Counts_CDS_df_copy_melt_filter$transcript_id
head(selected_tx); length(selected_tx) # 16955


write.table(data.frame(transcript_id=selected_tx), file = paste0(rsem_res, 'MostExpressedTranscripts.txt'), quote = F, sep = '\t', row.names = F, col.names = T)


# selected_tx >> the most expressed protein-coding transcripts per gene (start=0 and end=10^6 removed)
# CDS_len_ProteinCoding_canonical >> protein-coding transcripts that are canonical (start=0 and end=10^6 removed)
nb_selected_tx_that_are_canonical = length(intersect(selected_tx, CDS_len_ProteinCoding_canonical$transcript_id))
(nb_selected_tx_that_are_canonical / length(selected_tx))*100 # 78.71424% of selected tx that are canonical


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# scanBamWhat()
params = ScanBamParam(what=scanBamWhat(), tag=c('AS', 'XS', 'XN', 'XM', 'XO', 'XG', 'NM', 'MD', 'YT'))

bam_files = BAM_df$bam
bam_f = bam_files[28]
bam_f

stats_all = data.frame(matrix(ncol=0, nrow=0))
stats_all
for (bam_f in bam_files){
    print('------------------------------')
    print(bam_f)
    print(Sys.time())

    sampleID = BAM_df[BAM_df$bam == bam_f,]$Sample

    bam_n = do.call(rbind, str_split(bam_f, '/'))[,8]

    final_bam = paste0(rsem_res, str_replace(bam_n, 'bam', 'filtered.bam'))
    final_bam

    b_S <- readGAlignments(bam_f, use.names=TRUE, param=params)

    init_lib_depth = length(b_S)
    init_lib_depth

    read_tx = paste0(names(b_S), '_', seqnames(b_S))
    head(read_tx); length(read_tx)
    length(which(duplicated(read_tx)))
    mcols(b_S)$read_tx = read_tx

    # filter reads mapped multiple times to a given transcript - keep only one read
    b_S_filtered = b_S[-which(duplicated(mcols(b_S)$read_tx)),]
    b_S_filtered
    length(b_S_filtered)

    b_S_filtered_selected_tx = b_S_filtered[seqnames(b_S_filtered) %in% selected_tx]
    length(b_S_filtered_selected_tx)

    export(b_S_filtered_selected_tx, BamFile(final_bam))

    stats_df = data.frame(Sample=sampleID, Init_BAM=init_lib_depth, ReadsFiltered_BAM=length(b_S_filtered), FinalSelectedTx_BAM=length(b_S_filtered_selected_tx))
    stats_df

    if (nrow(stats_all) == 0){
        stats_all = stats_df
    }else{
        stats_all = rbind(stats_all, stats_df)
    }

    print(Sys.time())
    print('------------------------------')
}

# Init_BAM >> lib depth of BAM
# ReadsFiltered_BAM >> after removing reads present more than once in a transcript (only first occurrence of the read kept)
# FinalSelectedTx_BAM >> keeping only reads mapped to the selected tx
write.table(stats_all, file = paste0(rsem_res, 'Stats_BAM.txt'), quote=F, sep='\t', row.names=F, col.names=T)

