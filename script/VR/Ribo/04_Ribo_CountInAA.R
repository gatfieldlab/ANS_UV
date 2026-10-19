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


######################################################################
### Find pausing sites in genes using a method similar to the one
### published in Stein et al. Nature 2022
### PART ONE

library(tidyverse)
library(DESeq2)
library(org.Mm.eg.db)
library(cowplot)
library(Logolas)
library(stringr)
library(dplyr)
library(devtools)
library(data.table)
library(GenomicRanges)
library(biomaRt)



########################
### REFERENCE GENOME ### 
########################
GENOME='Mmusculus.GRCm39.111'
DB='/data/databases/mouse/'
DB_GTF=paste0(DB, 'gtf/')
GTF_GENOME=paste0(DB_GTF, GENOME, '.sorted.gtf')

GTF = rtracklayer::import(GTF_GENOME)
GTF_gene = GTF[GTF$type == 'gene_id', ]
GTF_gene

GTF_tx = GTF[GTF$type == 'transcript', ]
GTF_gene_tx = data.frame(gene_id = GTF_tx$gene_id, transcript_id = GTF_tx$transcript_id, mgi_symbol=GTF_tx$gene_name)
head(GTF_gene_tx); dim(GTF_gene_tx)

mart = useEnsembl('ensembl', dataset = 'mmusculus_gene_ensembl', version = '111')


### Histones list by Sebastien
histones = read.csv(paste0('/path_to/ANS/', 'final_histone_list_to_keep.csv'), sep=',', header=T)
histones$X=NULL
head(histones); dim(histones) # 19

histones_to_remove = read.csv(paste0('/path_to/ANS/', 'histones_to_remove.csv'), sep=',', header=T)
histones_to_remove$X=NULL
head(histones_to_remove); dim(histones_to_remove) # 108
### Histones list by Sebastien



CDS_coo = read.csv(paste0(DB, 'Mmusculus.GRCm39.111.startStop.length.txt'), sep='\t', header=FALSE)
names(CDS_coo) = c('transcript_id', 'length_cDNA', 'start', 'end')
CDS_coo$start = CDS_coo$start + 1
CDS_coo$end = CDS_coo$end + 1
head(CDS_coo); dim(CDS_coo)
range(CDS_coo$start)
range(CDS_coo$end)

########################
### SAMPLE INFO ###
########################
path_home='/data/vricci/projects/Sebastien/ANS/'
SampleInfo = read.csv(paste0(path_home, 'Samples_info_final.txt'), sep=' ', header=F)
names(SampleInfo) = c('Library', 'Condition', 'Treatment', 'Barcode')
SampleInfo$name = paste0(SampleInfo$Library, '_', SampleInfo$Barcode)
SampleInfo$ID = paste0(SampleInfo$Condition, '_', SampleInfo$Treatment)
SampleInfo$Type = 'Mono'
SampleInfo[grep('PD', SampleInfo$Library),]$Type = 'Di'
head(SampleInfo); dim(SampleInfo)


######################################################################
######################################################################
######################################################################
### SummarizeOverlaps(CDS)
path_bam = '/path_to/ANS_HD_UV_LD/RSEM/sub_bam/'
Counts_CDS_df = read.csv(paste0(path_bam, 'SummarizeOverlaps_CDS_raw.txt'), sep='\t', header=T)
head(Counts_CDS_df); dim(Counts_CDS_df)
Counts_CDS = Counts_CDS_df[-c(1:2, ncol(Counts_CDS_df))]
rownames(Counts_CDS) = Counts_CDS_df$transcript_id
head(Counts_CDS); dim(Counts_CDS)


LOST = setdiff(SampleInfo$name, colnames(Counts_CDS))
LOST
SampleInfo = SampleInfo[!SampleInfo$name %in% LOST,]
rownames(SampleInfo) = NULL
head(SampleInfo); dim(SampleInfo)


######################################################################
######################################################################
######################################################################
### codon counts (1 codon window)
# first combine counts for each codon (e.g. pos 142, 143, 144 >> pos 142 == codon '142')
SGA_codons_all = fread('/path_to/ANS_UV_project/chip_cor/Mmusculus.GRCm39.111.transcriptCodonsNorm.sga', sep='\t', header=F)
head(SGA_codons_all)
length(unique(SGA_codons_all$V1)) # 53071

# counts are not used - only transcript ID, AA and position
SGA_codons_all_df = data.frame(SGA_codons_all)
names(SGA_codons_all_df) = c('transcript_id', 'codon', 'pos', 'strand', 'nb')
SGA_codons_all_df = SGA_codons_all_df[c('transcript_id', 'pos', 'codon')]
head(SGA_codons_all_df); dim(SGA_codons_all_df)

SGA_codons_all_GR = GRanges(
    seqnames = SGA_codons_all$V1,
    ranges = IRanges(SGA_codons_all$V3, width=3),
    codon = SGA_codons_all$V2,
    strand= '+'
    )
SGA_codons_all_GR


### SGA files A-site clean (one tx per gene (the most expressed), reads occur only once per tx)
path_sga = '/path_to/ANS_HD_UV_LD/RSEM/SGA/' # Disome - Asite offset +45
path_sga_new = '/path_to/ANS/SGA/' # Monosome - Asite offset using RiboWaltz
import_sga <- c(list.files(path=path_sga_new, pattern="PF.*Asite.sga$", full.names=T), list.files(path=path_sga, pattern="PD.*.Asite.sga$", full.names=T))
head(import_sga); length(import_sga)

SGA_codon = SGA_codons_all_df
head(SGA_codon); dim(SGA_codon)
f_sga = import_sga[1]
for (f_sga in import_sga){

    tmp = do.call(rbind, str_split(f_sga, '/'))[,ncol(do.call(rbind, str_split(f_sga, '/')))]
    ID = do.call(rbind, str_split(tmp, '[.]'))[,1]
    print(ID)
    
    SGA = fread(f_sga)
    head(SGA); dim(SGA)
    SGA = SGA %>% dplyr::filter(!V1 %in% histones_to_remove$ensembl_transcript_id)
    head(SGA); dim(SGA)
    length(unique(SGA$V1))

    SGA_GR = GRanges(
        seqnames = SGA$V1,
        ranges = IRanges(SGA$V3, width=1),
        strand= '+',
        count=SGA$V5
        )
    SGA_GR

    overlaps = findOverlaps(SGA_codons_all_GR, SGA_GR)
    overlaps

    SGA_codons_merge = SGA_codons_all_GR[queryHits(overlaps)]
    SGA_codons_merge$pos = start(SGA_GR[subjectHits(overlaps)])
    SGA_codons_merge

    SGA_merge = SGA_GR[subjectHits(overlaps)]
    SGA_merge$codon = SGA_codons_merge$codon
    SGA_merge$pos = start(SGA_codons_merge)
    SGA_merge

    SGA_merge_df = data.frame(transcript_id=seqnames(SGA_merge), pos=SGA_merge$pos, codon=SGA_merge$codon, count=SGA_merge$count)
    head(SGA_merge_df); dim(SGA_merge_df)

    SGA_merge_df_aggr = aggregate(SGA_merge_df$count, by=list(SGA_merge_df$transcript_id, SGA_merge_df$pos, SGA_merge_df$codon), FUN=sum)
    names(SGA_merge_df_aggr) = names(SGA_merge_df)
    SGA_merge_df_aggr = SGA_merge_df_aggr[order(SGA_merge_df_aggr$transcript_id, SGA_merge_df_aggr$pos),]
    names(SGA_merge_df_aggr)[4] = ID
    head(SGA_merge_df_aggr); dim(SGA_merge_df_aggr)

    SGA_codon = full_join(SGA_codon, SGA_merge_df_aggr)
    SGA_codon[is.na(SGA_codon)] = 0

}
head(SGA_codon); dim(SGA_codon) # 25489181

SGA_codon_info = left_join(SGA_codon, GTF_gene_tx)
head(SGA_codon_info); dim(SGA_codon_info)
ncol(SGA_codon_info)
SGA_codon_info = SGA_codon_info[c(1:3, (ncol(SGA_codon_info)-1), ncol(SGA_codon_info), 4:(ncol(SGA_codon_info)-2))]
head(SGA_codon_info); dim(SGA_codon_info)
SGA_codon_info_filtered = SGA_codon_info[which(rowSums(SGA_codon_info[-c(1:5)]) != 0),]
head(SGA_codon_info_filtered); dim(SGA_codon_info_filtered) # 5019632
# 20645072 sites with rowSums == 0

path_output='/data/vricci/projects/Sebastien/ANS/CodonCounts/'
system(paste0('mkdir -p ', path_output))
fwrite(SGA_codon_info_filtered, paste0(path_output, 'CodonCounts.txt'), sep ='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)

# w StartCodon
head(SGA_codon_info_filtered); dim(SGA_codon_info_filtered)
SGA_codon_info_filtered_freq = aggregate(SGA_codon_info_filtered[-c(1:5)], by=list(SGA_codon_info_filtered$codon), FUN=sum)
names(SGA_codon_info_filtered_freq)[1] = 'codon'
head(SGA_codon_info_filtered_freq); dim(SGA_codon_info_filtered_freq) # 62 (codons)
fwrite(SGA_codon_info_filtered_freq, paste0(path_output, 'CodonCountsFreq.txt'), sep ='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)

#
head(SGA_codon_info_filtered); dim(SGA_codon_info_filtered) # 5019632
head(CDS_coo); dim(CDS_coo) # 149076
CDS_coo_copy = CDS_coo[c('transcript_id', 'start')]
names(CDS_coo_copy)[2] = 'pos'
head(CDS_coo_copy); dim(CDS_coo_copy)


# startCodon only
SGA_codon_info_filtered_StartCodonOnly = semi_join(SGA_codon_info_filtered, CDS_coo_copy, by=c('transcript_id', 'pos'))
head(SGA_codon_info_filtered_StartCodonOnly); dim(SGA_codon_info_filtered_StartCodonOnly) # 9549
unique(SGA_codon_info_filtered_StartCodonOnly$codon)

SGA_codon_info_filtered_StartCodonOnly_freq = aggregate(SGA_codon_info_filtered_StartCodonOnly[-c(1:5)], by=list(SGA_codon_info_filtered_StartCodonOnly$codon), FUN=sum)
names(SGA_codon_info_filtered_StartCodonOnly_freq)[1] = 'codon'
head(SGA_codon_info_filtered_StartCodonOnly_freq); dim(SGA_codon_info_filtered_StartCodonOnly_freq) # 3
fwrite(SGA_codon_info_filtered_StartCodonOnly_freq, paste0(path_output, 'CodonCountsFreq_StartCodonOnly.txt'), sep ='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)


## wo startCodon
SGA_codon_info_filtered_woStartCodon = anti_join(SGA_codon_info_filtered, CDS_coo_copy, by=c('transcript_id', 'pos'))
head(SGA_codon_info_filtered_woStartCodon); dim(SGA_codon_info_filtered_woStartCodon) # 5010083

SGA_codon_info_filtered_woStartCodon_freq = aggregate(SGA_codon_info_filtered_woStartCodon[-c(1:5)], by=list(SGA_codon_info_filtered_woStartCodon$codon), FUN=sum)
names(SGA_codon_info_filtered_woStartCodon_freq)[1] = 'codon'
head(SGA_codon_info_filtered_woStartCodon_freq); dim(SGA_codon_info_filtered_woStartCodon_freq) # 62
fwrite(SGA_codon_info_filtered_woStartCodon_freq, paste0(path_output, 'CodonCountsFreq_woStartCodon.txt'), sep ='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)
