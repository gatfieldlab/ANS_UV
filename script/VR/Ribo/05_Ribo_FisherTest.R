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
### PART TWO

print(.libPaths())
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
library(tibble)


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

CDS_coo = read.csv(paste0(DB, 'Mmusculus.GRCm39.111.startStop.length.txt'), sep='\t', header=FALSE)
names(CDS_coo) = c('transcript_id', 'length_cDNA', 'start', 'end')
CDS_coo$start = CDS_coo$start + 1
CDS_coo$end = CDS_coo$end + 1
head(CDS_coo); dim(CDS_coo)


########################
### SAMPLE INFO ###
########################
path_home='/path_to/ANS/'
SampleInfo = read.csv(paste0(path_home, 'Samples_info_final.txt'), sep=' ', header=F)
names(SampleInfo) = c('Library', 'Condition', 'Treatment', 'Barcode')
SampleInfo$name = paste0(SampleInfo$Library, '_', SampleInfo$Barcode)
SampleInfo$ID = paste0(SampleInfo$Condition, '_', SampleInfo$Treatment)
SampleInfo$Type = 'Mono'
SampleInfo[grep('PD', SampleInfo$Library),]$Type = 'Di'
head(SampleInfo); dim(SampleInfo)

LOST = setdiff(SampleInfo$name, colnames(Counts_CDS)); LOST # "PF69_TAGAC" "PD63_CGTAA"
SampleInfo = SampleInfo[!SampleInfo$name %in% LOST,]
rownames(SampleInfo) = NULL
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
# !!! last column >> geneMeans !!!


Counts_CDS_mono = Counts_CDS[grep('PF', colnames(Counts_CDS))]
head(Counts_CDS_mono); dim(Counts_CDS_mono)

Counts_CDS_di = Counts_CDS[grep('PD', colnames(Counts_CDS))]
head(Counts_CDS_di); dim(Counts_CDS_di)


SampleInfo_mono = SampleInfo[grep('PF', SampleInfo$Library),]
head(SampleInfo_mono); dim(SampleInfo_mono)

SampleInfo_di = SampleInfo[grep('PD', SampleInfo$Library),]
head(SampleInfo_di); dim(SampleInfo_di)


dds_mono = DESeqDataSetFromMatrix(countData = Counts_CDS_mono, colData = SampleInfo_mono, design = ~ ID)
dds_mono <- DESeq(dds_mono)
resultsNames(dds_mono)

dds_di = DESeqDataSetFromMatrix(countData = Counts_CDS_di, colData = SampleInfo_di, design = ~ ID)
dds_di <- DESeq(dds_di)
resultsNames(dds_di)


### Most expressed tx
MostExprTx = read.csv(paste0(path_bam, 'MostExpressedTranscripts.txt'), sep='\t', header=T)
head(MostExprTx); dim(MostExprTx)

Counts_CDS_df_MostExprTx = Counts_CDS_df[Counts_CDS_df$transcript_id %in% MostExprTx$transcript_id,]
head(Counts_CDS_df_MostExprTx); dim(Counts_CDS_df_MostExprTx)
summary(Counts_CDS_df_MostExprTx$geneMeans)


### OUTPUTs of PART ONE (FindingPausesites_Part1.R)
path_output='/path_to/ANS/CodonCounts/'
### CodonCounts.txt
SGA_codon_info_filtered = data.frame(fread(paste0(path_output, 'CodonCounts.txt'), sep='\t', header=TRUE))
head(SGA_codon_info_filtered); dim(SGA_codon_info_filtered)


# w StartCodon - CodonCountsFreq.txt
SGA_codon_info_filtered_freq = data.frame(fread(paste0(path_output, 'CodonCountsFreq.txt'), sep='\t', header=TRUE))
head(SGA_codon_info_filtered_freq); dim(SGA_codon_info_filtered_freq) # 64 (codons)


# startCodon only -  CodonCountsFreq_StartCodonOnly.txt
SGA_codon_info_filtered_StartCodonOnly_freq = data.frame(fread(paste0(path_output, 'CodonCountsFreq_StartCodonOnly.txt'), sep='\t', header=TRUE))
head(SGA_codon_info_filtered_StartCodonOnly_freq); dim(SGA_codon_info_filtered_StartCodonOnly_freq) # 3


## wo startCodon - CodonCountsFreq_woStartCodon.txt
SGA_codon_info_filtered_woStartCodon_freq = data.frame(fread(paste0(path_output, 'CodonCountsFreq_woStartCodon.txt'), sep='\t', header=TRUE))
head(SGA_codon_info_filtered_woStartCodon_freq); dim(SGA_codon_info_filtered_woStartCodon_freq) # 62
### OUTPUTs of PART ONE (FindingPausesites_Part1.R)


# stop codon
SGA_codon_info_filtered[SGA_codon_info_filtered$codon == 'TGA',]
SGA_codon_info_filtered[SGA_codon_info_filtered$codon == 'TAA',]
SGA_codon_info_filtered[SGA_codon_info_filtered$codon == 'TAG',]



### Compare codons - differential analysis
contrasts_list = list(c('ID', 'ANS_LD_treat', 'ANS_LD_ctrl'),
                    c('ID', 'ANS_HD_treat', 'ANS_HD_ctrl'),
                    c('ID', 'UV_long_treat', 'UV_long_ctrl'),
                    c('ID', 'UV_short_treat', 'UV_short_ctrl'))
contrasts_list

cmin = 0
codonData = SGA_codon_info_filtered
rownames(codonData) = paste0(codonData$transcript_id, '_', codonData$pos)
head(codonData)

type='mono'
contr=1
for (type in c('mono', 'di')){
    if (type == 'mono'){
        cdsData = Counts_CDS_mono
        sampleAnno = SampleInfo_mono
    }else{
        cdsData = Counts_CDS_di
        sampleAnno = SampleInfo_di
        
    }
    rownames(sampleAnno) = sampleAnno$name
    head(sampleAnno)

    for (contr in seq(1:4)){
        contrasts = contrasts_list[[contr]]

        print(type)
        print(contrasts)

        id1 <- which(sampleAnno[,contrasts[1]] == contrasts[2]); id1; sampleAnno[sampleAnno[,contrasts[1]] == contrasts[2],]
        id2 <- which(sampleAnno[,contrasts[1]] == contrasts[3]); id2; sampleAnno[sampleAnno[,contrasts[1]] == contrasts[3],]
        head(cdsData[,rownames(sampleAnno)[id1]])
        head(cdsData[,rownames(sampleAnno)[id2]])
        cdsCounts1 <- round(rowMeans(cdsData[,rownames(sampleAnno)[id1]])); head(cdsCounts1)
        cdsCounts2 <- round(rowMeans(cdsData[,rownames(sampleAnno)[id2]]))
        head(codonData[,rownames(sampleAnno)[id1]])
        head(codonData[,rownames(sampleAnno)[id2]])
        codonCounts1 <- round(rowMeans(codonData[,rownames(sampleAnno)[id1]])); head(codonCounts1)
        codonCounts2 <- round(rowMeans(codonData[,rownames(sampleAnno)[id2]]))
        codCounts.mean <- data.frame(codonData[,1:5],
                                        Sample_treat = codonCounts1,
                                        Sample_ctrl = codonCounts2,
                                        Sample_treat_tot = cdsCounts1[codonData[,1]],
                                        Sample_ctrl_tot = cdsCounts2[codonData[,1]]
                                        )
        head(codCounts.mean); dim(codCounts.mean) # 4844109
        codCounts.kept <- codCounts.mean[unique(which(codCounts.mean$Sample_ctrl != 0), which(codCounts.mean$Sample_treat != 0)), ]
        head(codCounts.kept); dim(codCounts.kept) # 1218829

        codCounts.kept$mini_tot = codCounts.kept$Sample_ctrl + codCounts.kept$Sample_treat
        head(codCounts.kept); dim(codCounts.kept)
        summary(codCounts.kept$mini_tot)

        codCounts.kept_filtered = codCounts.kept[codCounts.kept$mini_tot >=3,]
        head(codCounts.kept_filtered); dim(codCounts.kept_filtered) # 319844
        summary(codCounts.kept_filtered$mini_tot)


        oddRatio <- vector()
        pValues <- vector()
        for(I in 1:dim(codCounts.kept_filtered)[1]){
            if (I %in% seq(1, dim(codCounts.kept_filtered)[1], by=10000)){
                print(I)
            }
            x1 <- matrix(ncol = 2, nrow = 2)
            x1[1,1] <- codCounts.kept_filtered[I,6]
            x1[2,1] <- codCounts.kept_filtered[I,7]
            x1[1,2] <- codCounts.kept_filtered[I,8] - codCounts.kept_filtered[I,6]
            x1[2,2] <- codCounts.kept_filtered[I,9] - codCounts.kept_filtered[I,7]
            ft <- fisher.test(x=x1)
            # tt <- t.test(x=x1, )
            pValues <- c(pValues, ft$p.value)
            oddRatio <- c(oddRatio, ft$estimate)
        }
        padj <- p.adjust(pValues)
        codCounts.kept_filtered <- cbind(codCounts.kept_filtered,
                                oddRatios = log2(oddRatio),
                                pValues = pValues, padj = padj)
        colnames(codCounts.kept_filtered) <- c(colnames(codonData)[1:5],
                                        contrasts[2:3],
                                        paste(contrasts[2:3], "cds", sep="."),
                                        'mini_tot', "oddRatios", "pValue", "padj")
        codCounts.kept_filtered$oddRatios[which(abs(codCounts.kept_filtered$oddRatios) == Inf)] <- NA
        head(codCounts.kept_filtered)

        fwrite(codCounts.kept_filtered, paste0(path_output, 'FisherTest_', paste0(contrasts[2], '_vs_', contrasts[3]), '_', type,  '.txt'), sep ='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)

    }
}


