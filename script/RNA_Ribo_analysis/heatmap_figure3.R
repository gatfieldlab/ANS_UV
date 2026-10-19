##################################
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
library(DESeq2)
library(ggplot2)
library(rstatix)
library(data.table)


#################
#######import .sga file 
#################

setwd("~/path/to/file.sga")



mart_export <- read.csv("C:/path/mart_export.txt") #file generated on ensembl biomart to have gene names and gene ID
colnames(mart_export)[1]="gene"
colnames(mart_export)[3]="gene_name"


gene_transcript_list <- read.csv("~/path/to/gene_transcript.csv", header=FALSE) #gene+transcript IDs table
colnames(gene_transcript_list)[1]="gene"
colnames(gene_transcript_list)[2]="transcripts"
# 
histones_to_remove <- read.csv("~/path/histones_to_remove.csv") #histone list to remove to avoid family multimapping


Mmusculus.GRCm39.111.startStop.APPRIS <- read.delim("~/path/Mmusculus.GRCm39.111.startStop.APPRIS.bed", header=FALSE) #APPRIS file GRCm39.111
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[2]="start"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[3]="stop"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[4]="type"
anno_file <- as.data.frame(Mmusculus.GRCm39.111.startStop.APPRIS)

setwd("~/path/to/wd")


existing <- function(cols, df) { #function to check if columns are presents in a df 
  cols[cols %in% names(df)]
}

files_tx <- list.files(pattern = "Most")
MostExpressedTranscripts <- as.data.frame(fread(files_tx), sep='\t', header=T) #table of most expressed transcripts

files_RSEM <- list.files(pattern = "SummarizeOverlaps_CDS_raw")  #import counts

SummarizeOverlaps_CDS_raw <- as.data.frame(fread(files_RSEM), sep='\t', header=T)
SummarizeOverlaps_CDS_raw <- SummarizeOverlaps_CDS_raw[SummarizeOverlaps_CDS_raw$transcript_id %in% MostExpressedTranscripts$transcript_id, ]
SummarizeOverlaps_CDS_raw <- subset(SummarizeOverlaps_CDS_raw, !(transcript_id %in% histones_to_remove$ensembl_transcript_id))

ribo_counts <- select(SummarizeOverlaps_CDS_raw, -geneMeans)
rownames(ribo_counts) <- NULL
ribo_counts <- select(ribo_counts, transcript_id, gene_id, contains("PF")) %>%
  column_to_rownames(., "gene_id")


rna_counts <- as.data.frame(fread("~/path/Counts_cDNA.txt", sep='\t', header=T)) #import RNA counts

rna_counts <- subset(rna_counts, transcript_id  %in% ribo_counts$transcript_id) 
names(rna_counts) <- gsub("UV", "u", names(rna_counts))


rownames(rna_counts) <- NULL


rna_counts <- rna_counts%>%column_to_rownames(., "gene_id") %>%
  select(-transcript_id)

ribo_counts <- ribo_counts%>%
  select(-transcript_id)
names(ribo_counts) <- gsub("PF", "F", names(ribo_counts))


count_matrix <- cbind(
  RNA  = rna_counts,
  Ribo = ribo_counts
)

count_matrix <- select(count_matrix, !contains("F70"))
count_matrix <- select(count_matrix, !contains("long"))

names(count_matrix) <- gsub("([0-9])([A-Za-z])", "\\1_\\2", names(count_matrix))
names(count_matrix) <- gsub("[.]", "_", names(count_matrix))


coldata <- data.frame(
  sample = colnames(count_matrix)
) %>%
  mutate(assay=str_split_fixed(sample,"_", 3)[,1])

tmp_coldata <-coldata


tmp_coldata <- coldata[coldata$assay=="Ribo",]
infos_tmp <- str_split_fixed(tmp_coldata$sample, "_",3)
tmp_coldata$condition <- infos_tmp[,2]
tmp_coldata$rep <- infos_tmp[,3]
tmp_coldata$rep <- infos_tmp[,3]


tmp_coldata[tmp_coldata$rep=="ATCGT",]$rep <- 1
tmp_coldata[tmp_coldata$rep=="AGCTA",]$rep <- 2
tmp_coldata[tmp_coldata$rep=="CGTAA",]$rep <- 3
tmp_coldata[tmp_coldata$rep=="CTAGA",]$rep <- 4
tmp_coldata[tmp_coldata$rep=="GATCA",]$rep <- 5
tmp_coldata[tmp_coldata$rep=="GCATA",]$rep <- 6
tmp_coldata[tmp_coldata$rep=="TAGAC",]$rep <- 7
tmp_coldata[tmp_coldata$rep=="TCTAG",]$rep <- 8

tmp_coldata[tmp_coldata$condition=="F69",]$condition <- "LD"
tmp_coldata[tmp_coldata$condition=="F79",]$condition <- "HD"
tmp_coldata[tmp_coldata$condition=="F80",]$condition <- "UV"

tmp_coldata$compare <- infos_tmp[,3]

tmp_coldata[tmp_coldata$compare %in% c("ATCGT","AGCTA","CGTAA","CTAGA") & tmp_coldata$condition=="LD" ,]$compare <- 1
tmp_coldata[tmp_coldata$compare %in% c("GATCA","GCATA","TAGAC","TCTAG") & tmp_coldata$condition=="LD" ,]$compare <- 2
tmp_coldata[tmp_coldata$compare %in% c("ATCGT","AGCTA","CGTAA","CTAGA") & tmp_coldata$condition=="HD" ,]$compare <- 3
tmp_coldata[tmp_coldata$compare %in% c("GATCA","GCATA","TAGAC","TCTAG") & tmp_coldata$condition=="HD" ,]$compare <- 4
tmp_coldata[tmp_coldata$compare %in% c("ATCGT","AGCTA","CGTAA","CTAGA") & tmp_coldata$condition=="UV" ,]$compare <- 5
tmp_coldata[tmp_coldata$compare %in% c("GATCA","GCATA","TAGAC","TCTAG") & tmp_coldata$condition=="UV" ,]$compare <- 6

ribo_counts <- tmp_coldata

tmp_coldata <- coldata[coldata$assay=="RNA",]
infos_tmp <- str_split_fixed(tmp_coldata$sample, "_",3)
tmp_coldata$condition <- infos_tmp[,3]
tmp_coldata[tmp_coldata$condition=="short",]$condition <- "UV"

tmp_coldata$compare <- str_split_fixed(tmp_coldata$sample, "_",2)[,2]

tmp_coldata[grepl("^D\\d+_LD$", tmp_coldata$compare), "compare"] <- 1
tmp_coldata[grepl("^A\\d+_LD$", tmp_coldata$compare), "compare"] <- 2
tmp_coldata[grepl("^D\\d+_HD$", tmp_coldata$compare), "compare"] <- 3
tmp_coldata[grepl("^A\\d+_HD$", tmp_coldata$compare), "compare"] <- 4
tmp_coldata[grepl("^N\\d+_short$", tmp_coldata$compare), "compare"] <- 5
tmp_coldata[grepl("^u\\d+_short$", tmp_coldata$compare), "compare"] <- 6


tmp_coldata$rep <- substr(infos_tmp[,2],2,2)

rna_counts <- tmp_coldata

coldata <- rbind(ribo_counts, rna_counts)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = ~ compare + assay + compare:assay
)

dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)



mean_counts <- norm_counts %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(
    -gene,
    names_to = "sample",
    values_to = "norm_count"
  ) %>%
  left_join(
    coldata,
    by = "sample"
  ) %>%
  mutate(info=paste0(assay,"_", condition, "_", compare))%>%
  group_by(gene, info) %>%
  summarise(
    mean_norm_count = mean(norm_count),
    .groups = "drop"
  )%>%
  pivot_wider(
    names_from = info,
    values_from = mean_norm_count
  )



df_long <- mean_counts %>%
  pivot_longer(
    -gene,
    names_to = "sample",
    values_to = "value"
  )

df_long$value <- log2(df_long$value+1)

df2 <- df_long %>%
  mutate(
    assay     = ifelse(str_detect(sample, "^RNA"), "RNA", "Ribo"),
    dose      = str_extract(sample, "LD|HD|UV"),
    replicate = str_extract(sample, "\\d+$")
  )

df_wide <- df2 %>%
  pivot_wider(
    id_cols    = c(gene, dose, replicate),
    names_from = assay,
    values_from = value
  )

df_diff <- df_wide %>%
  mutate(value = Ribo - RNA, 
         sample = paste0("TE_", dose,"_", replicate)) %>%
  select(gene,sample, value)

df_long <- rbind(df_long,df_diff)



df <- df_long %>% pivot_wider(values_from = value,names_from = sample)

cols <- setdiff(colnames(df), "gene")

library(pheatmap)

##prepare the files of genes with up sites

import_file <- list.files("~/path/to/fisher_test_results/", pattern = "vs")

link <- "~/path/to/fisher_test_results"
FisherTest_tx_Codoncontent <- read.table(paste0(link, "FisherTest_tx_Codoncontent.txt"), sep="\t", header = TRUE)


Mmusculus.GRCm39.111.startStop.APPRIS <- read.delim("~/path/Mmusculus.GRCm39.111.startStop.IDS.bed", header=FALSE)
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[1]="merge"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[2]="start"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[3]="stop"
Mmusculus.GRCm39.111.startStop.APPRIS<-Mmusculus.GRCm39.111.startStop.APPRIS[-4]
split_merge <-as.data.frame(str_split_fixed(Mmusculus.GRCm39.111.startStop.APPRIS$merge, "[|]", 2))

Mmusculus.GRCm39.111.startStop.APPRIS<-cbind(Mmusculus.GRCm39.111.startStop.APPRIS, split_merge)
names(Mmusculus.GRCm39.111.startStop.APPRIS) <- c("merge", "start", "stop", "gene_id", "transcript_id")

for (i in import_file){   #generates position/sites violin plot figure
  tmp_file <-read.csv(paste0(link,i), sep="\t", header = TRUE)
  name_file <- str_split_fixed(i, "[.]", 2)[1]
  
  tmp_file <- subset(tmp_file, !(transcript_id %in% list_histone_to_remove$ensembl_transcript_id))
  
  
  tmp_file <- list(tmp_file,Mmusculus.GRCm39.111.startStop.APPRIS) %>% purrr::reduce(left_join, by="transcript_id")%>%
    dplyr::select(-gene_id.y) %>%
    dplyr::rename(gene_id = gene_id.x)%>%
    mutate(relative_to_start = pos-start)
  
  tmp_file <- subset(tmp_file,transcript_id  != "ENSMUST00000038863") #REMOVE LARS2 :as described in methods
  tmp_file <- subset(tmp_file,transcript_id  != "ENSMUST00000004389") #remove GRCC10 :as described in methods
  tmp_file <- subset(tmp_file,transcript_id  != "ENSMUST00000094434") #remove Ftl1 :as described in methods
  
  n_gene <- length(unique(tmp_file$mgi_symbol))
  n_site <- length(tmp_file$pos)
  vec_all <- c(n_gene,n_site )
  
  
  sites_oddratios_sig <- tmp_file[tmp_file$padj<=0.05,]
  sites_oddratios_sig$oddRatios[is.na(sites_oddratios_sig$oddRatios)] <- -10
  
  sites_oddratios_ns_sig <- tmp_file[tmp_file$padj>0.05,]
  sites_oddratios_ns_sig$oddRatios[is.na(sites_oddratios_ns_sig$oddRatios)] <- -10
  
  n_gene_sig <- length(unique(sites_oddratios_sig$mgi_symbol))
  n_site_sig <- length(sites_oddratios_sig$pos)
  
  n_gene_ns_sig <- length(unique(sites_oddratios_ns_sig$mgi_symbol))
  n_site_ns_sig <- length(sites_oddratios_ns_sig$pos)
  
  vec_sig <- c(n_gene_sig,n_site_sig )
  vec_ns_sig <- c(n_gene_ns_sig,n_site_ns_sig )
  
  
  
  sites_oddratios_sig$group <- ifelse(sites_oddratios_sig$oddRatios > 0, "up", "down")
  sites_oddratios_ns_sig$group <- ifelse(sites_oddratios_ns_sig$oddRatios > 0, "up", "down")
  
  
  assign(paste0(i, "_sig_data"), sites_oddratios_sig)
  assign(paste0(i, "ns_sig_data"), sites_oddratios_sig)

  
  # Compute range of x values
  x_min <- min(sites_oddratios_sig$relative_to_start + 1, sites_oddratios_sig$relative_to_start + 1)
  x_max <- max(sites_oddratios_sig$relative_to_start + 1, sites_oddratios_sig$relative_to_start + 1)
  log_range <- floor(log10(x_min)):ceiling(log10(x_max))
  custom_breaks <- unlist(lapply(log_range, function(p) seq(10^p, 10^(p + 1), length.out = 10)))
 
  
  violin_plot_log <- ggplot(sites_oddratios_sig, aes(x = (relative_to_start) , y = group, fill = group)) +
    geom_violin(stat = "ydensity") +
    theme_minimal() +
    labs(title = "Pause Sites Relative to Start Codon",
         y = "") +
    scale_fill_manual(values = c("up" = "#1f77b4", "down" = "#ff7f0e")) +
    theme_clean() +
    scale_x_log10(
      limits = c(1, 14000),  # Force range from 1 to 14000
      breaks = custom_breaks,
      labels = scales::label_number(accuracy = 1)
    ) +
    ggtitle(paste0(
      "n_up= ", length(subset(sites_oddratios_sig, oddRatios > 0)$transcript_id),
      "   n_down= ", length(subset(sites_oddratios_sig, oddRatios < 0)$transcript_id)
    ))
  
  ggsave(violin_plot_log, file=paste0("~/path/figure.pdf"),width = 6, height = 3, units = "in")
  
  
  tmp_df <- data.frame(all=vec_all,
                       sig=vec_sig)
  row.names(tmp_df) <- c("ngene", "nsite")
  
  assign(paste0(name_file, "_stats"), tmp_df) 
  
  assign(name_file, tmp_file) 
  assign(paste0(name_file, "_sig"), sites_oddratios_sig) 
  assign(paste0(name_file, "_ns_sig"), sites_oddratios_ns_sig) 
  
  
}


for (i in ls(pattern = "^Fisher.*(mono|di)_sig$")) {
  
  tmp_file <- get(i)
  tmp_file <- tmp_file[tmp_file$group=="up", ]
  
  tmp_file <- tmp_file %>% 
    group_by(gene_id) %>%
    dplyr::count(.)
  
  assign(paste0(i, "_count_table"), tmp_file)
  
}

list_count_tables <- ls(pattern = "count_table")
list_count_tables <- list_count_tables[!str_detect(list_count_tables, "long")]


cols <- setdiff(colnames(df), "gene")
# Select numeric matrix only once
mat_all <- as.matrix(select(df[, -1],!contains("TE") ))

# Global range
global_range <- range(mat_all, na.rm = TRUE)
# global_range <- c(5,16)

# Global breaks
global_breaks <- seq(global_range[1], global_range[2], length.out = 100)

cols <- colnames(df)[-1]

rm(FisherTest_UV_long_treat_vs_UV_long_ctrl_mono_sig_count_table)
color_palette <- colorRampPalette(c("white", "red"))(n = 100)

for (col in cols) {
  
  ID <- str_split_fixed(col, "_", 3)[, 2]
  
  # Get the count_table for this ID
  
  # Order df by the current column
  tmp <- df[order(-df[[col]]), c("gene", col)]
  # Create the main heatmap (p)
  mat <- as.matrix(tmp[[col]])
  rownames(mat) <- tmp$gene
  colnames(mat) <- col
  
  p <- pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    show_rownames = FALSE,
    show_colnames =  FALSE,
    main = col,
    # breaks = global_breaks
  )
  
  # Save the gene order from tmp
  gene_order <- tmp$gene
  
  # Reorder count_table to match the gene order of p
  count_table <- get(ls(pattern = paste0(ID, ".*mono_sig_count_table")))
  mat_df <- data.frame(mat)
  mat_df <- mat_df %>% rownames_to_column("gene_id")
  
  count_table <- list(mat_df,count_table) %>% 
    purrr::reduce(left_join, by="gene_id")%>% 
    replace(is.na(.), 0) %>%
    select(gene_id, n) %>%
    column_to_rownames("gene_id")
  
  count_table_ordered <- count_table[match(gene_order, rownames(count_table)), ]
  
  # Plot the count_table heatmap with the same gene order
  color_palette <- colorRampPalette(c("white", "red"))(n = 100)
  
  # Plot the heatmap with the custom color palette
  p_count <- pheatmap(
    count_table_ordered,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    show_rownames = FALSE,
    show_colnames = FALSE,
    
    main = paste0(col, "_count_table"),
    color = color_palette  # Apply the custom color palette,
    
  )
  
  # Assign both heatmaps
  assign(paste0(col, "_heatmap"), p)
  assign(paste0(col, "_count_heatmap"), p_count)
  
  ggsave(p, file=paste0("~/save/path/",col, "_figure.pdf"), width = 2, height =5, units = "in")
  ggsave(p_count, file=paste0("~/save/path/",col, "_figure.pdf"), width = 2, height =5, units = "in")
  
}

