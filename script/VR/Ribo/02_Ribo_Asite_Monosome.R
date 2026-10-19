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
library(devtools, lib.loc='/home/vricci/miniconda3/lib/R/library')
# install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
library(riboWaltz)


path_home='/path_to/ANS_HD_UV_LD/'
rsem_data=paste0(path_home, 'RSEM/')
rsem_res = paste0(rsem_data, 'sub_bam/')
sga_res = paste0('/path_to/ANS/', 'SGA/')
system(paste0('mkdir -p ', sga_res))


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
########################
### SAMPLE INFO ###
########################
path_home='/path_to/ANS/'
SampleInfo = read.csv(paste0(path_home, 'Samples_info_final.txt'), sep=' ', header=F)
names(SampleInfo) = c('Library', 'Condition', 'Treatment', 'Barcode')
SampleInfo$sample = paste0(SampleInfo$Library, '_', SampleInfo$Barcode)
SampleInfo$ID = paste0(SampleInfo$Condition, '_', SampleInfo$Treatment)
SampleInfo$Type = 'Mono'
SampleInfo[grep('PD', SampleInfo$Library),]$Type = 'Di'
head(SampleInfo); dim(SampleInfo)


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# .bam files - PROJECTION ON TRANSCRIPTOME OF THE GENOME MAPPING with RSEM
BAM=list.files(path=rsem_res, pattern='*.transcript.sorted.filtered.bam$', full.names=T, recursive=TRUE)
BAM = str_replace(BAM, '//', '/')
BAM_split = do.call(rbind,str_split(BAM, '/'))
Sample = do.call(rbind, str_split(BAM_split[,7], '[.]'))[,1] 
Sample # e.g. PD63_AGCTA
BAM # e.g. /path_to/ANS_HD_UV_LD/RSEM/PD63/AGCTA/PD63_AGCTA.transcript.sorted.bam
BAM_df = data.frame(Sample = Sample, bam = BAM)
BAM_df

Bam_Files = BamFileList(BAM_df$bam)
names(Bam_Files) = BAM_df$Sample
###



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



### RiboWaltz
GTF_anno = create_annotation(GTF_GENOME)
GTF_anno
# transcript, length_transcript, length_utr5, length_cds, length_utr3


reads_list <- bamtolist(bamfolder = rsem_res, annotation = GTF_anno)
length(reads_list)
head(reads_list[1])


psite_offset <- psite(reads_list, flanking = 6, extremity = "5end")
psite_offset

psite_offset_copy = psite_offset
psite_offset_copy$sample = str_replace(psite_offset_copy$sample, '.transcript.sorted.filtered', '')
psite_offset_copy

psite_offset_info = left_join(data.frame(psite_offset_copy), SampleInfo[c('sample', 'ID')])
psite_offset_info = psite_offset_info[c((ncol(psite_offset_info)-1), ncol(psite_offset_info), 1:(ncol(psite_offset_info)-2))]
head(psite_offset_info); dim(psite_offset_info)

write.table(psite_offset_info, file = paste0(sga_res, 'Psite_output.txt'), quote = F, sep = '\t', row.names = F, col.names = T)

psite_offset_info_final = psite_offset_info %>% dplyr::filter(ID %in% c('ANS_HD_ctrl', 'ANS_HD_treat', 
                                            'ANS_LD_ctrl', 'ANS_LD_treat',
                                            'UV_short_ctrl', 'UV_short_treat'))

psite_offset_info_final = psite_offset_info_final %>% dplyr::filter(grepl('PF', sample))

head(psite_offset_info_final); dim(psite_offset_info_final)
write.table(psite_offset_info_final, file = paste0(sga_res, 'Psite_output_final.txt'), quote = F, sep = '\t', row.names = F, col.names = T)


reads_psite_list <- psite_info(reads_list, psite_offset)
length(reads_psite_list)
reads_psite_list[1]

names(reads_psite_list)
length(reads_psite_list)
x=1
for (x in seq(1, length(reads_psite_list))){
    print(x)
    sample = do.call(rbind, str_split(names(reads_psite_list[x]), '[.]'))[,1]
    reads_psite_list[[x]]$asite = reads_psite_list[[x]]$psite + 3

    sga_x = reads_psite_list[[x]][,c('transcript', 'asite')]
    sga_x

    sga_x_aggr = sga_x %>% count(transcript, asite)
    sga_x_aggr
    sga_x_aggr_order = sga_x_aggr[str_order(sga_x_aggr$transcript),]
    sga_x_aggr_order

    sga_x_final = data.frame(transcript_id=sga_x_aggr_order$transcript, cDNA='cDNA', Asite=sga_x_aggr_order$asite, strand='+', count=sga_x_aggr_order$n)
    head(sga_x_final); dim(sga_x_final)

    out_f=paste0(sample, '.mouse_cDNA.clean.Asite.sga')
    write.table(sga_x_final, file = paste0(sga_res, out_f), quote = F, sep = '\t', row.names = F, col.names = F)

}

