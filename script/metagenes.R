
library(patchwork)
library(ggplot2)
library(rlist)  
library(cowplot)
library(dplyr)
library(ggthemes)
library(readxl)
library(stats)
library(gridExtra)
library(grid)
library(tidyverse)
library(data.table)



#################
#######import .sga file
#################

setwd("~/path/to/sga_files")


mart_export <- read.csv("C:/Users/schamois/Downloads/mart_export.txt") #file generated on ensembl biomart to have gene names and gene ID
colnames(mart_export)[1]="gene"
colnames(mart_export)[3]="gene_name"


gene_transcript_list <- read.csv("~/gene_transcript_list.csv", header=FALSE) #gene+transcript IDs
colnames(gene_transcript_list)[1]="gene"
colnames(gene_transcript_list)[2]="transcripts"

#get start and stop info from APPRIS
Mmusculus.GRCm39.111.startStop.APPRIS <- read.delim("~/annotations/Mmusculus.GRCm39.111.startStop.APPRIS.bed", header=FALSE)
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[2]="start"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[3]="stop"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[4]="type"
anno_file <- as.data.frame(Mmusculus.GRCm39.111.startStop.APPRIS)

#build function to check if a col exist in a df (ctrl for missing samples)
existing <- function(cols, df) {
  cols[cols %in% names(df)]
}

#import count files 
files_tx <- list.files(pattern = "Most")
MostExpressedTranscripts <- as.data.frame(fread(files_tx), sep='\t', header=T)

files_RSEM <- list.files(pattern = "SummarizeOverlaps_CDS_raw")
SummarizeOverlaps_CDS_raw <- as.data.frame(fread(files_RSEM), sep='\t', header=T)
SummarizeOverlaps_CDS_raw <- SummarizeOverlaps_CDS_raw[SummarizeOverlaps_CDS_raw$transcript_id %in% MostExpressedTranscripts$transcript_id, ]
SummarizeOverlaps_CDS_raw <- dplyr::select(SummarizeOverlaps_CDS_raw, -geneMeans)

lib <- unique(str_split_fixed(names(SummarizeOverlaps_CDS_raw)[3:ncol(SummarizeOverlaps_CDS_raw)], "_", 2)[,1])
i="PF80"
for (i in lib) {
  
  lib <- i
  nt_cols <- c(paste0(lib, "_AGCTA"), paste0(lib, "_CTAGA"), paste0(lib, "_CGTAA"), paste0(lib, "_ATCGT"))
  ans_cols <- c(paste0(lib, "_GATCA"), paste0(lib, "_GCATA"), paste0(lib, "_TAGAC"), paste0(lib,"_TCTAG"))
  
  
  tmp_file <- select(SummarizeOverlaps_CDS_raw, transcript_id, starts_with(lib))
  
  mutate_exprs <- list()
  
  nt_cols_pres <- existing(nt_cols, tmp_file)
  ans_cols_pres <- existing(ans_cols, tmp_file)
 
  if (length(nt_cols_pres) > 0) {
    mutate_exprs$mean_nt <- rowMeans(tmp_file[, nt_cols_pres, drop = FALSE], na.rm = TRUE)
  }
  
  if (length(ans_cols_pres) > 0) {
    mutate_exprs$mean_ans <- rowMeans(tmp_file[, ans_cols_pres, drop = FALSE], na.rm = TRUE)
  }
  
  tmp_file <- tmp_file %>%
    mutate(!!!mutate_exprs)
  
  colnames(tmp_file)[ncol(tmp_file)-1]=paste0("mean_nt_", lib)
  colnames(tmp_file)[ncol(tmp_file)]=paste0("mean_ans_", lib)
  
  tmp_file_nt <- select(tmp_file, transcript_id, starts_with("mean_nt"))
  
  tmp_file_nt <- tmp_file_nt %>%
    dplyr::select(transcript_id, all_of(paste0("mean_nt_", lib))) %>%
    arrange(desc(.data[[paste0("mean_nt_", lib)]])) %>%
    slice_head(prop = 0.2) %>%
    dplyr::select(transcript_id)
  
  tmp_file_ans <- dplyr::select(tmp_file, transcript_id, starts_with("mean_ans"))
  
  tmp_file_ans <- tmp_file_ans %>%
    dplyr::select(transcript_id, all_of(paste0("mean_ans_", lib))) %>%
    arrange(desc(.data[[paste0("mean_ans_", lib)]])) %>%
    slice_head(prop = 0.2) %>%
    dplyr::select(., transcript_id)

  assign(paste0(lib, "_CDS_counts"), tmp_file)
  assign(paste0(lib, "_top_20perc_exp_nt"), tmp_file_nt)
  assign(paste0(lib, "_top_20perc_exp_treat"), tmp_file_ans)
  
  rm(tmp_file)
}


##filter multi-gene family histones

files_sga <- list.files("~/path/to/sga/", pattern = "Asite")

list_histone <-read.table("~/annotations/Histones_VR.txt", sep="\t", header = TRUE)

colnames(list_histone)[1] <- "gene_id"

Mmusculus.GRCm39.111.GeneLength <- read.delim("~/Clock uORF/RNAseq/Metacycle/HTseqcounts/Mmusculus.GRCm39.111.GeneLength.txt") #file for gene length info for RPKM
Mmusculus.GRCm39.111.GeneLength_mod<- dplyr::select(Mmusculus.GRCm39.111.GeneLength,gene_id,gene_name) 

list_histone <-list(Mmusculus.GRCm39.111.GeneLength_mod, list_histone) %>% purrr::reduce(inner_join, by="gene_id")

Mmusculus.GRCm39.111.startStop.APPRIS <- read.delim("C:/Users/schamois/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/schamois/ANS_UV/tools_tables/Mmusculus.GRCm39.111.startStop.APPRIS.bed", header=FALSE)
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[1]="merge"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[2]="start"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[3]="stop"
Mmusculus.GRCm39.111.startStop.APPRIS<-Mmusculus.GRCm39.111.startStop.APPRIS[-4]
split_merge <-as.data.frame(str_split_fixed(Mmusculus.GRCm39.111.startStop.APPRIS$merge, "[|]", 2))
Mmusculus.GRCm39.111.startStop.APPRIS<-cbind(Mmusculus.GRCm39.111.startStop.APPRIS, split_merge)
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[4]="gene_id"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[5]="transcript_id"


for (i in files_sga){
  lib<-substr(i,1,4)
  adapter<-substr(i,6,10)
  tmp_file <-read.table(paste0("~/BTF3 study/20250624_lib_ANSHD_UV30/sga/",i), header = FALSE)
  names(tmp_file) <- c("gene_id", "cDNA", "position", "strand", "counts")
  tmp_file <- dplyr::select(tmp_file, "gene_id", "position", "counts")
  tmp_file <- subset(tmp_file,gene_id != "ENSMUST00000038863") #REMOVE LARS2 explained why in methods
  tmp_file <- subset(tmp_file,gene_id != "ENSMUST00000004389") #remove GRCC10
  tmp_file <- subset(tmp_file,gene_id != "ENSMUST00000094434") #remove Ftl1
  assign(paste0(lib,"_",adapter,".Asite.sga"),tmp_file)
}

for(i in ls(pattern = "Asite")){
  
  tmp_file <- get(i)

  lib <- substr(i, 1, 10)
  
  list_histone_tmp <- list_histone
  names(list_histone_tmp) <- c("ENSEMBL_ID", "gene_name","gene_id", "biotype", "SYMBOL", "description", "canonical")
  
  tmp_file <- subset(tmp_file, gene_id %in% list_histone_tmp$gene_id)
  
  tmp_file <- list(tmp_file, list_histone_tmp)  %>% purrr::reduce(inner_join, by="gene_id")%>%
    group_by(SYMBOL, gene_id, ENSEMBL_ID)%>%
    summarise(sum_counts = sum(counts))%>%
    ungroup()%>%
    mutate(family = sub("[0-9]+$", "", SYMBOL))%>%   # extract family (remove trailing numbers)
    group_by(family) %>%
    slice_max(order_by = sum_counts , n = 1, with_ties = FALSE) %>% # keep top gene per family
    ungroup()
  
  assign(paste0(lib, "_histones"),tmp_file )
  
  
}

df_list <- lapply(ls(pattern = "(A|C|G|T)_histone"), get) #merge them in a df

histone_df <- df_list %>%  purrr::reduce(full_join, by="SYMBOL")

histone_to_keep <- data.frame(SYMBOL=histone_df$SYMBOL) %>% mutate(family = sub("[0-9]+$", "", SYMBOL))   # extract family (remove trailing numbers)

#select one representative per family based on alphabetical order of SYMBOL
final_histone_list_to_keep <- histone_to_keep %>%
  arrange(family, SYMBOL) %>%
  group_by(family) %>%
  slice(1) %>%
  ungroup()

histones_to_remove <- subset(list_histone, !(mgi_symbol %in% final_histone_list_to_keep$SYMBOL))


write.csv(final_histone_list_to_keep, "~/BTF3 study/20250624_lib_ANSHD_UV30/final_histone_list_to_keep.csv")
write.csv(histones_to_remove, "~/BTF3 study/20250624_lib_ANSHD_UV30/histones_to_remove.csv")

df_names <- ls(pattern = "Asite")   
df_list  <- mget(df_names)                  

for(i in df_names){
  tmp_file <- get(i)
  tmp_file <- subset(tmp_file, !(gene_id %in% histones_to_remove$ensembl_transcript_id))
  assign(i,tmp_file )
  
  
}


lib <- unique(str_split_fixed(names(SummarizeOverlaps_CDS_raw)[3:ncol(SummarizeOverlaps_CDS_raw)], "_", 2)[,1])

sample_list<-ls(pattern = "Asite.sga")

for (i in sample_list){   #file precessing -> remove col and merge transcript/pos
  tmp_file <- get(i)
  
  adapter<-substr(i,6,10)
  lib<-substr(i,1,4)
  names(tmp_file)<-c("transcripts",
                     paste0("position_", adapter),
                     paste0("counts_",adapter))
  
  tmp_file$transcripts <- paste(tmp_file$transcripts, sep="_" , tmp_file$position)
  
  
  tmp_file<-dplyr::select(tmp_file,
                          transcripts,
                          paste0("counts_",adapter))
  
  assign(i, tmp_file)
  
}

libs <- c("^PD80_.*.Asite.sga","^PF80_.*.Asite.sga",
          "^PD79_.*.Asite.sga","^PF79_.*.Asite.sga", 
          "^PD63_.*.Asite.sga","^PF69_.*.Asite.sga")

for (i in libs){
  
  lib<-substr(i,2,5)
  list_file_merge <- lapply(ls(pattern = i), get) #merge them in a df
  
  merge_gene_sga_DESeq2 <- list_file_merge %>%
    purrr::reduce(full_join, by = "transcripts") %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  
  assign(paste0(lib, "_merge_gene_sga_DESeq2"),merge_gene_sga_DESeq2)
  
  write.csv(merge_gene_sga_DESeq2, file=paste0("~/BTF3 study/20240318_library ANS and UV/tables/", lib, "_merge_gene_sga_DESeq2.csv"))
  
  rm(merge_gene_sga_DESeq2)

}


Mmusculus.GRCm39.111.GeneLength <- read.delim("~/path/to/Mmusculus.GRCm39.111.GeneLength.txt") #file for tx and gene length info 
Mmusculus.GRCm39.111.GeneLength_mod<- dplyr::select(Mmusculus.GRCm39.111.GeneLength,gene_id,gene_name) 

file_list <- ls(pattern = ("(PD63|PF69|PD64|PF70|PF79|PF80|PD79|PD80)_merge_gene_sga_DESeq2" ))


file_merged <- ls(pattern = ("(PD63|PF69|PD64|PF70|PF79|PF80|PD79|PD80)_merge_gene_sga_DESeq2" ))
for (i in file_merged){
  lib<-substr(i,1,4)
  
  
  tmp_file <- get(i)
  
  tmp_file <- column_to_rownames(tmp_file,"transcripts")
  depth_file <- apply(tmp_file, 2 , sum)
  depth_file_gene <- sweep(tmp_file, 2, depth_file, `/`)*1000000
  
  depth_file_gene <- tibble::rownames_to_column(depth_file_gene, "transcript_position") #remplace le row.names par une colonne que l'on appelle "gene_id"
  
  
  split_trasncript_name <-   as.data.frame(str_split_fixed(depth_file_gene$transcript_position, "_",2))
  colnames(split_trasncript_name)[1]="transcripts"
  colnames(split_trasncript_name)[2]="position" 
  
  
  lib_Asite_merge <- cbind(split_trasncript_name,depth_file_gene )
  lib_Asite_merge_gene <- list(lib_Asite_merge, gene_transcript_list) %>% purrr::reduce(inner_join, by="transcripts")
  colnames(lib_Asite_merge_gene)[ncol(lib_Asite_merge_gene)]="gene_id"
  
  lib_Asite_merge_gene <- list(lib_Asite_merge_gene, Mmusculus.GRCm39.111.GeneLength_mod) %>% purrr::reduce(inner_join, by="gene_id")
  
  assign(paste0(lib, "_Asite_merge_gene_norm_depth"), lib_Asite_merge_gene)
  rm(lib_Asite_merge_gene)
  rm(tmp_file)
  rm(depth_file)
  rm(depth_file_gene)
  rm(split_trasncript_name)
  rm(lib_Asite_merge)
  rm(i)
}

#define the column groups
nt_cols <- c("counts_AGCTA", "counts_CTAGA", "counts_CGTAA", "counts_ATCGT")
ans_cols <- c("counts_GATCA", "counts_GCATA", "counts_TAGAC","counts_TCTAG" )

list_file <- ls(pattern = "_Asite_merge_gene_norm_depth")
for (i in list_file) {
  tmp_file <- get(i)
  
  mutate_exprs <- list()
  
  
  nt_cols_pres <- existing(nt_cols, tmp_file)
  ans_cols_pres <- existing(ans_cols, tmp_file)

  if (length(nt_cols_pres) > 0) {
    mutate_exprs$mean_nt <- rowMeans(tmp_file[, nt_cols_pres, drop = FALSE], na.rm = TRUE)
  }
  
  if (length(ans_cols_pres) > 0) {
    mutate_exprs$mean_ans <- rowMeans(tmp_file[, ans_cols_pres, drop = FALSE], na.rm = TRUE)
  }
 
  tmp_file <- tmp_file %>%
    mutate(!!!mutate_exprs)
  
  assign(i, tmp_file)
  rm(tmp_file)
}



Mmusculus.GRCm39.111.startStop.APPRIS <- read.delim("~/annotations/Mmusculus.GRCm39.111.startStop.IDS.bed", header=FALSE)
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[1]="merge"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[2]="start"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[3]="stop"
Mmusculus.GRCm39.111.startStop.APPRIS<-Mmusculus.GRCm39.111.startStop.APPRIS[-4]
split_merge <-as.data.frame(str_split_fixed(Mmusculus.GRCm39.111.startStop.APPRIS$merge, "[|]", 2))

Mmusculus.GRCm39.111.startStop.APPRIS<-cbind(Mmusculus.GRCm39.111.startStop.APPRIS, split_merge)

Mmusculus.GRCm39.111.startStop.APPRIS <- dplyr::select(Mmusculus.GRCm39.111.startStop.APPRIS, 
                                                       V2,
                                                       start,
                                                       stop)
Mmusculus.GRCm39.111.startStop.APPRIS$start <- Mmusculus.GRCm39.111.startStop.APPRIS$start+1
Mmusculus.GRCm39.111.startStop.APPRIS$stop <- Mmusculus.GRCm39.111.startStop.APPRIS$stop+1


colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[1]="transcripts"


list_file <- ls(pattern = "_Asite_merge_gene_norm_depth")
for(i in list_file){
  
  tmp_file <- get(i)
  
  tmp_file <- list(tmp_file, Mmusculus.GRCm39.111.startStop.APPRIS) %>% purrr::reduce(inner_join, by="transcripts")
  
  
  tmp_file$position<-as.numeric(tmp_file$position)
  tmp_file$relative_to_start <- tmp_file$position - tmp_file$start
  
  assign(i, tmp_file)
  
  rm(tmp_file)
}

GenesExonsLengths <- read.delim("~/annotations/LengthsExonsIntronsGenes.txt")
Genelength <- rstatix::select(GenesExonsLengths, GeneID, ExonLength) #genelength for rpkm -> exon based because the counts are based on exons
colnames(Genelength)[1] = 'gene_id'
colnames(Genelength)[2] = 'length'


list_file <- ls(pattern = ("(PD63|PF69|PD64|PF70|PF79|PF80|PD79|PD80)_Asite_merge_gene_norm_depth" ))
for(i in list_file){
  
  lib<-substr(i,1,4)
  
  tmp_file <- get(i)
  tmp_file <-   list(tmp_file,Genelength) %>% purrr::reduce(inner_join, by="gene_id")
  tmp_file$relative_to_stop= (tmp_file$relative_to_start+tmp_file$start)-tmp_file$stop
  
  assign(i,tmp_file)
  write.csv(i, file=paste0("~/BTF3 study/20240318_library ANS and UV/tables/",i,".csv"))
  
  rm(tmp_file)
}


file_list <- ls(pattern = "(PF69|PF70|PF79|PF80)_Asite_merge_gene_norm_depth$")

for(i in file_list){
  
  tmp_file <- get(i)
  lib<-substr(i,1,4)
  
  start_df <- tmp_file %>%
    dplyr::group_by(relative_to_start) %>%
    summarise(sum_counts_nt=sum(mean_nt),
              sum_counts_ans=sum(mean_ans))
  
  stop_df <- tmp_file %>%
    group_by(relative_to_stop) %>%
    summarise(sum_counts_nt=sum(mean_nt),
              sum_counts_ans=sum(mean_ans))
  
  metagene_plot_start_top <- ggplot(start_df) + 
    geom_line( aes(relative_to_start,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_start,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-75,175)+
    # ylim(0,37000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")+
    coord_cartesian(ylim = c(15000, 37000))
  
  
  metagene_plot_start_bot<- ggplot(start_df) + 
    geom_line( aes(relative_to_start,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_start,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-75,175)+
    # ylim(0,37000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")+
    coord_cartesian(ylim = c(0, 5000))
  # scale_y_break(c(5000, 15000), scales=2) 
  
  # ggsave(metagene_plot_start, file=paste0("~/BTF3 study/20240318_library ANS and UV/figures/primary_alignment/",lib,"_start_pileup_metagene.pdf"),width = 5, height = 5, units = "in")
  
  metagene_plot_start <- grid.arrange(metagene_plot_start_top,metagene_plot_start_bot,nrow=2 )
  
  metagene_plot_stop_top <- ggplot(stop_df) + 
    geom_line( aes(relative_to_stop,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_stop,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-175,75)+
    # ylim(0,37000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")+
    coord_cartesian(ylim = c(15000, 37000))
  
  metagene_plot_stop_bot <- ggplot(stop_df) + 
    geom_line( aes(relative_to_stop,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_stop,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-175,75)+
    # ylim(0,37000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")+
    coord_cartesian(ylim = c(0, 5000))
  
  # metagene_plot_merge <- grid.arrange(metagene_plot_start_top,
  #                                     metagene_plot_stop_top,
  #                                     metagene_plot_start_bot,
  #                                     metagene_plot_stop_bot,
  #                                    nrow=2 )
  # 
  metagene_plot_merge <- plot_grid(
    metagene_plot_start_top,
    metagene_plot_stop_top,
    metagene_plot_start_bot,
    metagene_plot_stop_bot,
    nrow = 2,
    align = "hv"
  )
  
  
  # ggsave(metagene_plot_stop, file=paste0("~/BTF3 study/20240318_library ANS and UV/figures/primary_alignment/",lib,"_stop_pileup_metagene.pdf"),width = 5, height = 5, units = "in")
  
  # metagene_plot_merge <- grid.arrange(metagene_plot_stop_top,metagene_plot_stop_bot,nrow=2 )
  
  # metagene_plot_start+metagene_plot_stop
  
  
  ggsave(metagene_plot_merge, file=paste0("~/BTF3 study/20250624_lib_ANSHD_UV30/figure/",lib,"_start_stop_pileup_metagene.pdf"),width = 15, height = 10, units = "in")
  
  
  
  
  metagene_plot_start_zoom <- ggplot(start_df) + 
    geom_line( aes(relative_to_start,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_start,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-1000,6000)+
    ylim(0,10000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")
  # coord_cartesian(ylim = c(0, 6000))
  
  
  # ggsave(metagene_plot_start, file=paste0("~/BTF3 study/20240318_library ANS and UV/figures/primary_alignment/",lib,"_start_pileup_metagene.pdf"),width = 5, height = 5, units = "in")
  
  
  metagene_plot_stop_zoom <- ggplot(stop_df) + 
    geom_line( aes(relative_to_stop,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_stop,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-6000,1000)+
    ylim(0,10000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")
  
  # ggsave(metagene_plot_stop, file=paste0("~/BTF3 study/20240318_library ANS and UV/figures/primary_alignment/",lib,"_stop_pileup_metagene.pdf"),width = 5, height = 5, units = "in")
  
  metagene_plot_merge_zoom <- metagene_plot_start_zoom+metagene_plot_stop_zoom
  
  
  ggsave(metagene_plot_merge_zoom, file=paste0("~/BTF3 study/20250624_lib_ANSHD_UV30/figure/",lib,"_start_stop_pileup_metagene_zoom.pdf"),width = 15, height = 10, units = "in")
  
  rm(tmp_file)
}


file_list <- ls(pattern = "(PD63|PD64|PD79|PD80)_Asite_merge_gene_norm_depth$")
i="PD63_Asite_merge_gene_norm_depth"

for(i in file_list){
  
  tmp_file <- get(i)
  lib<-substr(i,1,4)
  
  start_df <- tmp_file %>%
    group_by(relative_to_start) %>%
    summarise(sum_counts_nt=sum(mean_nt),
              sum_counts_ans=sum(mean_ans))
  
  stop_df <- tmp_file %>%
    group_by(relative_to_stop) %>%
    summarise(sum_counts_nt=sum(mean_nt),
              sum_counts_ans=sum(mean_ans))
  
  metagene_plot_start_bot <- ggplot(start_df) + 
    geom_line( aes(relative_to_start,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_start,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-75,175)+
    # ylim(0,10000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")+
    coord_cartesian(ylim = c(0, 5500))+
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  
  
  metagene_plot_start_top <- ggplot(start_df) + 
    geom_line( aes(relative_to_start,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_start,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-75,175)+
    # ylim(0,10000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")+
    coord_cartesian(ylim = c(8000, 10000))+
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  
  
  
  # ggsave(metagene_plot_start, file=paste0("~/BTF3 study/20240318_library ANS and UV/figures/primary_alignment/",lib,"_start_pileup_metagene.pdf"),width = 5, height = 5, units = "in")
  
  
  metagene_plot_stop_bot <- ggplot(stop_df) + 
    geom_line( aes(relative_to_stop,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_stop,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-175,75)+
    # ylim(0,10000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")+
    coord_cartesian(ylim = c(0, 5500))+
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  
  
  metagene_plot_stop_top <- ggplot(stop_df) + 
    geom_line( aes(relative_to_stop,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_stop,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-175,75)+
    # ylim(0,10000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")+
    coord_cartesian(ylim = c(8000, 10000))+
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  

  
  
  metagene_plot_merge <- plot_grid(
    metagene_plot_start_top,
    metagene_plot_stop_top,
    metagene_plot_start_bot,
    metagene_plot_stop_bot,
    nrow = 2,
    align = "hv"
  )
  
  
  
  ggsave(metagene_plot_merge, file=paste0("~/BTF3 study/20250624_lib_ANSHD_UV30/figure/",lib,"_start_stop_pileup_metagene.pdf"),width = 15, height = 10, units = "in")
  
  
  
  
  metagene_plot_start_zoom <- ggplot(start_df) + 
    geom_line( aes(relative_to_start,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_start,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-1000,6000)+
    ylim(0,10000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")
  
  
  # ggsave(metagene_plot_start, file=paste0("~/BTF3 study/20240318_library ANS and UV/figures/primary_alignment/",lib,"_start_pileup_metagene.pdf"),width = 5, height = 5, units = "in")
  
  
  metagene_plot_stop_zoom <- ggplot(stop_df) + 
    geom_line( aes(relative_to_stop,sum_counts_nt ,color="nt")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_CGTAA_merge,color="nt")) +
    geom_line( aes(relative_to_stop,sum_counts_ans,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GATCA_merge,color="ans")) +
    # geom_line( aes(corrected_to_start,mean_counts_PD63_274_GCATA_merge,color="ans")) +
    scale_color_manual(values=c("seagreen","black"))+
    xlim(-6000,1000)+
    ylim(0,10000)+
    theme_classic()+
    ylab("Average of relative counts gene and per position")+
    xlab("position relative to START (nt)")+
    ggtitle(paste0(lib, "_metagene plot on all genes"))+
    theme(legend.position="none")
  
  # ggsave(metagene_plot_stop, file=paste0("~/BTF3 study/20240318_library ANS and UV/figures/primary_alignment/",lib,"_stop_pileup_metagene.pdf"),width = 5, height = 5, units = "in")
  
  metagene_plot_merge_zoom <- metagene_plot_start_zoom+metagene_plot_stop_zoom
  
  
  ggsave(metagene_plot_merge_zoom, file=paste0("~/BTF3 study/20250624_lib_ANSHD_UV30/figure/",lib,"_start_stop_pileup_metagene_zoom.pdf"),width = 20, height = 10, units = "in")
  
}

#cumulative expression plot 
file_list <- ls(pattern = "(PD63|PF69|PD64|PF70|PF79|PF80|PD79|PD80)_Asite_merge_gene_norm_depth$")
i="PD79_Asite_merge_gene_norm_depth"
for (i in file_list) {
  tmp_file <- get(i)
  
  merge_gene <- tmp_file %>%
    group_by(relative_to_start) %>%
    summarise(sum_counts_nt = sum(mean_nt),
              sum_counts_ans = sum(mean_ans))
  
  df_sorted_dmso <- merge_gene %>%
    arrange(relative_to_start) %>%
    mutate(cum_expr = cumsum(sum_counts_nt))
  
  df_sorted_ans <- merge_gene %>%
    arrange(relative_to_start) %>%
    mutate(cum_expr = cumsum(sum_counts_ans))
  
  df_sorted_dmso <- subset(df_sorted_dmso, relative_to_start >= 0)
  df_sorted_ans <- subset(df_sorted_ans, relative_to_start >= 0)
  
  # Compute range of x values
  x_min <- min(df_sorted_dmso$relative_to_start + 1, df_sorted_ans$relative_to_start + 1)
  x_max <- max(df_sorted_dmso$relative_to_start + 1, df_sorted_ans$relative_to_start + 1)
  
  # Generate breaks: 10 per log unit
  log_range <- floor(log10(x_min)):ceiling(log10(x_max))
  custom_breaks <- unlist(lapply(log_range, function(p) seq(10^p, 10^(p + 1), length.out = 10)))
  
  if (grepl("PF69", i)) {
    color_nt <- "#909090"
    color_treated <- "#F1811E"
  } else if (grepl("PF79", i)) {
    color_nt <- "#909090"
    color_treated <- "#D4007C"
  } else if (grepl("PF80", i)) {
    color_nt <- "#909090"
    color_treated <- "#00ADEB"
  } else if (grepl("PD63", i)) {
    color_nt <- "#000000"
    color_treated <- "#FF000D"
  } else if (grepl("PD79", i)) {
    color_nt <- "#000000"
    color_treated <- "#F000FF"
  } else if (grepl("PD80", i)) {
    color_nt <- "#000000"
    color_treated <- "#0037FF"
  }
  
  
  plot <- ggplot() +
    geom_point(data = df_sorted_dmso, aes(x = (relative_to_start + 1), y = cum_expr),
               color = color_nt) +
    geom_point(data = df_sorted_ans, aes(x = (relative_to_start + 1), y = cum_expr),
               color = color_treated) +
    scale_x_log10(breaks = custom_breaks, labels = scales::label_number(accuracy = 1)) +
    ggtitle(paste0(i, "_cumul_plot")) +
    theme_classic() +
    theme(legend.position = "none")
  
  assign(paste0(i, "_plot"), plot)
  ggsave(plot, file = paste0("~/BTF3 study/20250624_lib_ANSHD_UV30/figure/", i, "_cumulplot.pdf"),
         width = 4, height = 10, units = "in")
}

