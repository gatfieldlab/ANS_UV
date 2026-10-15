library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggthemes)
library(readxl)
library(stats)
library(gridExtra)
library(grid)
library(tidyverse)
library(data.table)
library(rstatix)
library(scales)

execute <- "Yes"  # if metagene plots are to be prepared (takes a lot of time) otherwise put anything

# ------------------------------------------------------------------------------
# PATH CONFIGURATION (GENERIC GITHUB)
# ------------------------------------------------------------------------------
# Project root directory
base_dir <- getwd()

# Subfolder definitions
data_dir        <- file.path(base_dir, "data")
sga_dir         <- file.path(data_dir, "sga")
rsem_dir        <- file.path(data_dir, "rsem")
annot_dir       <- file.path(base_dir, "annotations")
output_dir      <- file.path(base_dir, "output")
tables_dir      <- file.path(output_dir, "tables")
figures_dir     <- file.path(output_dir, "figures")
primary_fig_dir <- file.path(figures_dir, "primary_alignment")

# Automatic directory creation if they do not exist
dirs_to_create <- c(data_dir, sga_dir, rsem_dir, annot_dir, output_dir, tables_dir, figures_dir, primary_fig_dir)
for (d in dirs_to_create) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 1. SINGLE-TIME LOADING OF ANNOTATIONS AND REFERENCE FILES
# ------------------------------------------------------------------------------
mart_export <- read.csv(file.path(annot_dir, "mart_export.txt"))
colnames(mart_export)[1] <- "gene"
colnames(mart_export)[3] <- "gene_name"

gene_transcript_list <- read.csv(file.path(annot_dir, "gene_transcript_list.csv"), header = FALSE)
colnames(gene_transcript_list)[1] <- "gene"
colnames(gene_transcript_list)[2] <- "transcripts"

# APPRIS Annotations BED (single load)
appris_bed_path <- file.path(annot_dir, "Mmusculus.GRCm39.111.startStop.APPRIS.bed")
Mmusculus.GRCm39.111.startStop.APPRIS <- read.delim(appris_bed_path, header = FALSE)
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[1] <- "merge"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[2] <- "start"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[3] <- "stop"
Mmusculus.GRCm39.111.startStop.APPRIS <- Mmusculus.GRCm39.111.startStop.APPRIS[-4]
split_merge <- as.data.frame(str_split_fixed(Mmusculus.GRCm39.111.startStop.APPRIS$merge, "[|]", 2))
Mmusculus.GRCm39.111.startStop.APPRIS <- cbind(Mmusculus.GRCm39.111.startStop.APPRIS, split_merge)
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[4] <- "gene_id"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[5] <- "transcript_id"

# APPRIS IDS Annotations (single load)
appris_ids_path <- file.path(annot_dir, "Mmusculus.GRCm39.111.startStop.IDS.bed")
APPRIS_IDS <- read.delim(appris_ids_path, header = FALSE)
colnames(APPRIS_IDS)[1] <- "merge"
colnames(APPRIS_IDS)[2] <- "start"
colnames(APPRIS_IDS)[3] <- "stop"
APPRIS_IDS <- APPRIS_IDS[-4]
split_merge_ids <- as.data.frame(str_split_fixed(APPRIS_IDS$merge, "[|]", 2))
APPRIS_IDS <- cbind(APPRIS_IDS, split_merge_ids)
APPRIS_IDS <- rstatix::select(APPRIS_IDS, V2, start, stop)
APPRIS_IDS$start <- APPRIS_IDS$start + 1
APPRIS_IDS$stop  <- APPRIS_IDS$stop + 1
colnames(APPRIS_IDS)[1] <- "transcripts"

# Gene Lengths (single load)
Mmusculus.GRCm39.111.GeneLength <- read.delim(file.path(annot_dir, "Mmusculus.GRCm39.111.GeneLength.txt"))
Mmusculus.GRCm39.111.GeneLength_mod <- rstatix::select(Mmusculus.GRCm39.111.GeneLength, gene_id, gene_name)

# Histone Annotations (single load)
list_histone <- read.table(file.path(annot_dir, "Histones_VR.txt"), sep = "\t", header = TRUE)
colnames(list_histone)[1] <- "gene_id"
list_histone <- list(Mmusculus.GRCm39.111.GeneLength_mod, list_histone) %>% purrr::reduce(inner_join, by = "gene_id")

# Exon Lengths (single load)
GenesExonsLengths <- read.delim(file.path(annot_dir, "LengthsExonsIntronsGenes.txt"))
Genelength <- rstatix::select(GenesExonsLengths, GeneID, ExonLength)
colnames(Genelength)[1] <- 'gene_id'
colnames(Genelength)[2] <- 'length'

# Helper function
existing <- function(cols, df) {
  cols[cols %in% names(df)]
}

# ------------------------------------------------------------------------------
# 2. SELECTION OF TOP EXPRESSED TRANSCRIPTS (RSEM)
# ------------------------------------------------------------------------------
files_tx <- list.files(rsem_dir, pattern = "Most", full.names = TRUE)
MostExpressedTranscripts <- as.data.frame(fread(files_tx), sep = '\t', header = TRUE)

files_RSEM <- list.files(rsem_dir, pattern = "SummarizeOverlaps_CDS_raw", full.names = TRUE)
SummarizeOverlaps_CDS_raw <- as.data.frame(fread(files_RSEM), sep = '\t', header = TRUE)
SummarizeOverlaps_CDS_raw <- SummarizeOverlaps_CDS_raw[SummarizeOverlaps_CDS_raw$transcript_id \%in\% MostExpressedTranscripts$transcript_id, ]
SummarizeOverlaps_CDS_raw <- select(SummarizeOverlaps_CDS_raw, -geneMeans)

lib <- unique(str_split_fixed(names(SummarizeOverlaps_CDS_raw)[3:ncol(SummarizeOverlaps_CDS_raw)], "_", 2)[, 1])

for (i in lib) {
  nt_cols  <- c(paste0(i, "_AGCTA"), paste0(i, "_CTAGA"), paste0(i, "_CGTAA"), paste0(i, "_ATCGT"))
  ans_cols <- c(paste0(i, "_GATCA"), paste0(i, "_GCATA"), paste0(i, "_TAGAC"), paste0(i, "_TCTAG"))
  
  tmp_file <- select(SummarizeOverlaps_CDS_raw, transcript_id, starts_with(i))
  
  mutate_exprs <- list()
  nt_cols_pres  <- existing(nt_cols, tmp_file)
  ans_cols_pres <- existing(ans_cols, tmp_file)
  
  if (length(nt_cols_pres) > 0) {
    mutate_exprs$mean_nt <- rowMeans(tmp_file[, nt_cols_pres, drop = FALSE], na.rm = TRUE)
  }
  if (length(ans_cols_pres) > 0) {
    mutate_exprs$mean_ans <- rowMeans(tmp_file[, ans_cols_pres, drop = FALSE], na.rm = TRUE)
  }
  
  tmp_file <- tmp_file %>% mutate(!!!mutate_exprs)
  
  colnames(tmp_file)[ncol(tmp_file)-1] <- paste0("mean_nt_", i)
  colnames(tmp_file)[ncol(tmp_file)]   <- paste0("mean_ans_", i)
  
  tmp_file_nt <- tmp_file %>%
    select(transcript_id, all_of(paste0("mean_nt_", i))) %>%
    arrange(desc(.data[[paste0("mean_nt_", i)]])) %>%
    slice_head(prop = 0.2) %>%
    select(transcript_id)
  
  tmp_file_ans <- tmp_file %>%
    select(transcript_id, all_of(paste0("mean_ans_", i))) %>%
    arrange(desc(.data[[paste0("mean_ans_", i)]])) %>%
    slice_head(prop = 0.2) %>%
    select(transcript_id)
  
  assign(paste0(i, "_CDS_counts"), tmp_file)
  assign(paste0(i, "_top_20perc_exp_nt"), tmp_file_nt)
  assign(paste0(i, "_top_20perc_exp_treat"), tmp_file_ans)
  
  rm(tmp_file)
}

# ------------------------------------------------------------------------------
# 3. SGA FILE PROCESSING AND FILTERING
# ------------------------------------------------------------------------------
files_sga <- list.files(sga_dir, pattern = "Asite")

for (i in files_sga) {
  lib <- substr(i, 1, 4)
  adapter <- substr(i, 6, 10)
  tmp_file <- read.table(file.path(sga_dir, i), header = FALSE)
  names(tmp_file) <- c("gene_id", "cDNA", "position", "strand", "counts")
  tmp_file <- select(tmp_file, "gene_id", "position", "counts")
  tmp_file <- subset(tmp_file, gene_id != "ENSMUST00000038863") # REMOVE LARS2
  tmp_file <- subset(tmp_file, gene_id != "ENSMUST00000004389") # remove GRCC10
  tmp_file <- subset(tmp_file, gene_id != "ENSMUST00000094434") # remove Ftl1
  
  assign(paste0(lib, "_", adapter, ".Asite.sga"), tmp_file)
}

# Histone filtering
for (i in ls(pattern = "Asite")) {
  tmp_file <- get(i)
  lib <- substr(i, 1, 10)
  
  list_histone_tmp <- list_histone
  names(list_histone_tmp) <- c("ENSEMBL_ID", "gene_name", "gene_id", "biotype", "SYMBOL", "description", "canonical")
  
  tmp_file <- subset(tmp_file, gene_id %in% list_histone_tmp$gene_id)
  
  tmp_file <- list(tmp_file, list_histone_tmp) %>% purrr::reduce(inner_join, by = "gene_id") %>%
    group_by(SYMBOL, gene_id, ENSEMBL_ID) %>%
    summarise(sum_counts = sum(counts), .groups = "drop") %>%
    mutate(family = sub("[0-9]+$", "", SYMBOL)) %>%
    group_by(family) %>%
    slice_max(order_by = sum_counts, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  assign(paste0(lib, "_histones"), tmp_file)
}

df_list <- lapply(ls(pattern = "(A|C|G|T)_histone"), get)
histone_df <- df_list %>% purrr::reduce(full_join, by = "SYMBOL")

histone_to_keep <- data.frame(SYMBOL = histone_df$SYMBOL) \%>\% mutate(family = sub("[0-9]+$", "", SYMBOL))

final_histone_list_to_keep <- histone_to_keep %>%
  arrange(family, SYMBOL) %>%
  group_by(family) %>%
  slice(1) %>%
  ungroup()

histones_to_remove <- subset(list_histone, !(mgi_symbol %in% final_histone_list_to_keep$SYMBOL))

write.csv(final_histone_list_to_keep, file.path(tables_dir, "final_histone_list_to_keep.csv"), row.names = FALSE)
write.csv(histones_to_remove, file.path(tables_dir, "histones_to_remove.csv"), row.names = FALSE)

df_names <- ls(pattern = "Asite")
for (i in df_names) {
  tmp_file <- get(i)
  tmp_file <- subset(tmp_file, !(gene_id %in% histones_to_remove$ensembl_transcript_id))
  assign(i, tmp_file)
}

# Column formatting and merging by position
sample_list <- ls(pattern = "Asite.sga")
for (i in sample_list) {
  tmp_file <- get(i)
  adapter <- substr(i, 6, 10)
  lib     <- substr(i, 1, 4)
  
  names(tmp_file) <- c("transcripts", paste0("position_", adapter), paste0("counts_", adapter))
  tmp_file$transcripts <- paste(tmp_file$transcripts, sep = "_", tmp_file[[paste0("position_", adapter)]])
  
  tmp_file <- select(tmp_file, transcripts, paste0("counts_", adapter))
  assign(i, tmp_file)
}

libs <- c("PD80_.*.Asite.sga", "PF80_.*.Asite.sga",
          "PD79_.*.Asite.sga", "PF79_.*.Asite.sga", 
          "PD63_.*.Asite.sga", "PF69_.*.Asite.sga",
          "PD64_.*.Asite.sga", "PF70_.*.Asite.sga")

for (i in libs) {
  lib <- substr(i, 1, 4)
  list_file_merge <- lapply(ls(pattern = i), get)
  
  merge_gene_sga_DESeq2 <- list_file_merge %>%
    purrr::reduce(full_join, by = "transcripts") %>%
    replace(is.na(.), 0)
  
  assign(paste0(lib, "_merge_gene_sga_DESeq2"), merge_gene_sga_DESeq2)
  write.csv(merge_gene_sga_DESeq2, file = file.path(tables_dir, paste0(lib, "_merge_gene_sga_DESeq2.csv")), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 4. SEQUENCING DEPTH NORMALIZATION (RPM)
# ------------------------------------------------------------------------------
file_merged <- ls(pattern = "(PD63|PF69|PD64|PF70|PF79|PF80|PD79|PD80)_merge_gene_sga_DESeq2")

for (i in file_merged) {
  lib <- substr(i, 1, 4)
  tmp_file <- get(i)
  
  tmp_file <- column_to_rownames(tmp_file, "transcripts")
  depth_file <- apply(tmp_file, 2, sum)
  depth_file_gene <- sweep(tmp_file, 2, depth_file, `/`) * 1000000
  depth_file_gene <- tibble::rownames_to_column(depth_file_gene, "transcript_position")
  
  split_trasncript_name <- as.data.frame(str_split_fixed(depth_file_gene$transcript_position, "_", 2))
  colnames(split_trasncript_name)[1] <- "transcripts"
  colnames(split_trasncript_name)[2] <- "position" 
  
  lib_Asite_merge <- cbind(split_trasncript_name, depth_file_gene)
  lib_Asite_merge_gene <- list(lib_Asite_merge, gene_transcript_list) %>% purrr::reduce(inner_join, by = "transcripts")
  colnames(lib_Asite_merge_gene)[ncol(lib_Asite_merge_gene)] <- "gene_id"
  
  lib_Asite_merge_gene <- list(lib_Asite_merge_gene, Mmusculus.GRCm39.111.GeneLength_mod) %>% purrr::reduce(inner_join, by = "gene_id")
  
  assign(paste0(lib, "_Asite_merge_gene_norm_depth"), lib_Asite_merge_gene)
}

# Mean calculation per condition (nt vs ans)
nt_cols  <- c("counts_AGCTA", "counts_CTAGA", "counts_CGTAA", "counts_ATCGT")
ans_cols <- c("counts_GATCA", "counts_GCATA", "counts_TAGAC", "counts_TCTAG")

list_file <- ls(pattern = "_Asite_merge_gene_norm_depth")

for (i in list_file) {
  tmp_file <- get(i)
  mutate_exprs <- list()
  
  nt_cols_pres  <- existing(nt_cols, tmp_file)
  ans_cols_pres <- existing(ans_cols, tmp_file)
  
  if (length(nt_cols_pres) > 0) {
    mutate_exprs$mean_nt <- rowMeans(tmp_file[, nt_cols_pres, drop = FALSE], na.rm = TRUE)
  }
  if (length(ans_cols_pres) > 0) {
    mutate_exprs$mean_ans <- rowMeans(tmp_file[, ans_cols_pres, drop = FALSE], na.rm = TRUE)
  }
  
  tmp_file <- tmp_file %>% mutate(!!!mutate_exprs)
  assign(i, tmp_file)
}

# Addition of relative positions to START and STOP sites
for (i in list_file) {
  tmp_file <- get(i)
  
  tmp_file <- list(tmp_file, APPRIS_IDS) %>% purrr::reduce(inner_join, by = "transcripts")
  tmp_file$position <- as.numeric(tmp_file$position)
  tmp_file$relative_to_start <- tmp_file$position - tmp_file$start
  
  tmp_file <- list(tmp_file, Genelength) %>% purrr::reduce(inner_join, by = "gene_id")
  tmp_file$relative_to_stop  <- (tmp_file$relative_to_start + tmp_file$start) - tmp_file$stop
  
  assign(i, tmp_file)
  write.csv(tmp_file, file = file.path(tables_dir, paste0(i, ".csv")), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 5. GENE-NORMALIZED AVERAGE METAGENE (0.015 MAX)
# ------------------------------------------------------------------------------
file_list <- ls(pattern = "(PD63|PF69|PD64|PF70|PF79|PF80|PD79|PD80)_Asite_merge_gene_norm_depth$")

if (execute == "yes") {
  for (i in file_list) {
    lib <- substr(i, 1, 4)
    tmp_file <- get(i)
    
    meta_SUB_count <- rstatix::select(tmp_file, gene_id, transcripts, gene_name, transcript_position, start, stop, relative_to_start, relative_to_stop, mean_nt, mean_ans, starts_with("counts"))
    
    av1 <- aggregate(mean_nt ~ gene_id, data = meta_SUB_count, FUN = sum)
    colnames(av1)[2] <- "sum_nt"
    av2 <- aggregate(mean_ans ~ gene_id, data = meta_SUB_count, FUN = sum)
    colnames(av2)[2] <- "sum_ans"
    
    meta_merge <- list(av1, av2) %>% purrr::reduce(full_join, by = "gene_id") 
    meta_merge_final <- list(meta_merge, meta_SUB_count) %>% purrr::reduce(full_join, by = "gene_id")
    
    meta_merge_final$nt_relative  <- meta_merge_final$mean_nt / meta_merge_final$sum_nt
    meta_merge_final$ans_relative <- meta_merge_final$mean_ans / meta_merge_final$sum_ans
    meta_merge_final <- meta_merge_final %>% drop_na()
    
    meta_merge_final_start <- meta_merge_final %>%
      group_by(relative_to_start) %>%
      summarise(mean_counts_nt  = sum(nt_relative) / n_distinct(meta_merge_final$gene_id),
                mean_counts_ans = sum(ans_relative) / n_distinct(meta_merge_final$gene_id))
    
    write.csv(meta_merge_final_start, file.path(tables_dir, paste0(lib, "_meta_merge_final_nohistone_start.csv")), row.names = FALSE)
    
    meta_merge_final_stop <- meta_merge_final %>%
      group_by(relative_to_stop) %>%
      summarise(mean_counts_nt  = sum(nt_relative) / n_distinct(meta_merge_final$gene_id),
                mean_counts_ans = sum(ans_relative) / n_distinct(meta_merge_final$gene_id))
    
    write.csv(meta_merge_final_stop, file.path(tables_dir, paste0(lib, "_meta_merge_final_stop.csv")), row.names = FALSE)
    
    metagene_plot_start <- ggplot(meta_merge_final_start) + 
      geom_line(aes(relative_to_start, mean_counts_nt, color = "nt")) +
      geom_line(aes(relative_to_start, mean_counts_ans, color = "ans")) +
      scale_color_manual(values = c("seagreen", "black")) +
      xlim(-75, 200) + ylim(0, 0.015) + theme_classic() +
      ylab("Average of relative counts gene and per position") +
      xlab("position relative to START (nt)") +
      ggtitle(paste0(lib, "_metagene plot on all genes"))
    
    metagene_plot_stop <- ggplot(meta_merge_final_stop) + 
      geom_line(aes(relative_to_stop, mean_counts_nt, color = "nt")) +
      geom_line(aes(relative_to_stop, mean_counts_ans, color = "ans")) +
      scale_color_manual(values = c("seagreen", "black")) +
      xlim(-100, 100) + ylim(0, 0.015) + theme_classic() +
      ylab("Average of relative counts gene and per position") +
      xlab("position relative to STOP (nt)") +
      ggtitle(paste0(lib, "_metagene plot on all genes"))
    
    start_stop_metagene <- metagene_plot_start + metagene_plot_stop
    ggsave(start_stop_metagene, file = file.path(primary_fig_dir, paste0(lib, "_start_stop_metagene.pdf")), width = 10, height = 5, units = "in")
  }
}

# ------------------------------------------------------------------------------
# 6. RAW PILEUP METAGENE (PF = 37000/5000 | PD = 10000/5500) + ZOOM
# ------------------------------------------------------------------------------
for (i in file_list) {
  tmp_file <- get(i)
  lib <- substr(i, 1, 4)
  
  start_df <- tmp_file %>%
    dplyr::group_by(relative_to_start) %>%
    summarise(sum_counts_nt = sum(mean_nt), sum_counts_ans = sum(mean_ans))
  
  stop_df <- tmp_file %>%
    group_by(relative_to_stop) %>%
    summarise(sum_counts_nt = sum(mean_nt), sum_counts_ans = sum(mean_ans))
  
  is_pf <- grepl("PF", lib)
  ylim_top <- if (is_pf) c(15000, 37000) else c(8000, 10000)
  ylim_bot <- if (is_pf) c(0, 5000) else c(0, 5500)
  
  metagene_plot_start_top <- ggplot(start_df) + 
    geom_line(aes(relative_to_start, sum_counts_nt, color = "nt")) +
    geom_line(aes(relative_to_start, sum_counts_ans, color = "ans")) +
    scale_color_manual(values = c("seagreen", "black")) + xlim(-75, 175) + theme_classic() +
    ylab("Average of relative counts gene and per position") + xlab("position relative to START (nt)") +
    ggtitle(paste0(lib, "_metagene plot on all genes")) + theme(legend.position = "none") +
    coord_cartesian(ylim = ylim_top) + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  
  metagene_plot_start_bot <- ggplot(start_df) + 
    geom_line(aes(relative_to_start, sum_counts_nt, color = "nt")) +
    geom_line(aes(relative_to_start, sum_counts_ans, color = "ans")) +
    scale_color_manual(values = c("seagreen", "black")) + xlim(-75, 175) + theme_classic() +
    ylab("Average of relative counts gene and per position") + xlab("position relative to START (nt)") +
    ggtitle(paste0(lib, "_metagene plot on all genes")) + theme(legend.position = "none") +
    coord_cartesian(ylim = ylim_bot) + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  
  metagene_plot_stop_top <- ggplot(stop_df) + 
    geom_line(aes(relative_to_stop, sum_counts_nt, color = "nt")) +
    geom_line(aes(relative_to_stop, sum_counts_ans, color = "ans")) +
    scale_color_manual(values = c("seagreen", "black")) + xlim(-175, 75) + theme_classic() +
    ylab("Average of relative counts gene and per position") + xlab("position relative to START (nt)") +
    ggtitle(paste0(lib, "_metagene plot on all genes")) + theme(legend.position = "none") +
    coord_cartesian(ylim = ylim_top) + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  
  metagene_plot_stop_bot <- ggplot(stop_df) + 
    geom_line(aes(relative_to_stop, sum_counts_nt, color = "nt")) +
    geom_line(aes(relative_to_stop, sum_counts_ans, color = "ans")) +
    scale_color_manual(values = c("seagreen", "black")) + xlim(-175, 75) + theme_classic() +
    ylab("Average of relative counts gene and per position") + xlab("position relative to START (nt)") +
    ggtitle(paste0(lib, "_metagene plot on all genes")) + theme(legend.position = "none") +
    coord_cartesian(ylim = ylim_bot) + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  
  metagene_plot_merge <- plot_grid(metagene_plot_start_top, metagene_plot_stop_top, metagene_plot_start_bot, metagene_plot_stop_bot, nrow = 2, align = "hv")
  ggsave(metagene_plot_merge, file = file.path(figures_dir, paste0(lib, "_start_stop_pileup_metagene.pdf")), width = 15, height = 10, units = "in")
  
  # Zoom Metagene Plots
  metagene_plot_start_zoom <- ggplot(start_df) + 
    geom_line(aes(relative_to_start, sum_counts_nt, color = "nt")) +
    geom_line(aes(relative_to_start, sum_counts_ans, color = "ans")) +
    scale_color_manual(values = c("seagreen", "black")) + xlim(-1000, 6000) + ylim(0, 10000) + theme_classic() +
    ylab("Average of relative counts gene and per position") + xlab("position relative to START (nt)") +
    ggtitle(paste0(lib, "_metagene plot on all genes")) + theme(legend.position = "none")
  
  metagene_plot_stop_zoom <- ggplot(stop_df) + 
    geom_line(aes(relative_to_stop, sum_counts_nt, color = "nt")) +
    geom_line(aes(relative_to_stop, sum_counts_ans, color = "ans")) +
    scale_color_manual(values = c("seagreen", "black")) + xlim(-6000, 1000) + ylim(0, 10000) + theme_classic() +
    ylab("Average of relative counts gene and per position") + xlab("position relative to START (nt)") +
    ggtitle(paste0(lib, "_metagene plot on all genes")) + theme(legend.position = "none")
  
  metagene_plot_merge_zoom <- metagene_plot_start_zoom + metagene_plot_stop_zoom
  zoom_width <- if (is_pf) 15 else 20
  ggsave(metagene_plot_merge_zoom, file = file.path(figures_dir, paste0(lib, "_start_stop_pileup_metagene_zoom.pdf")), width = zoom_width, height = 10, units = "in")
}

# ------------------------------------------------------------------------------
# 7. CUMULATIVE EXPRESSION PLOT (LOG10 SCALE)
# ------------------------------------------------------------------------------
for (i in file_list) {
  tmp_file <- get(i)
  
  merge_gene <- tmp_file %>%
    group_by(relative_to_start) %>%
    summarise(sum_counts_nt  = sum(mean_nt),
              sum_counts_ans = sum(mean_ans))
  
  df_sorted_dmso <- merge_gene %>% arrange(relative_to_start) %>% mutate(cum_expr = cumsum(sum_counts_nt)) %>% subset(relative_to_start >= 0)
  df_sorted_ans  <- merge_gene %>% arrange(relative_to_start) %>% mutate(cum_expr = cumsum(sum_counts_ans)) %>% subset(relative_to_start >= 0)
  
  x_min <- min(df_sorted_dmso$relative_to_start + 1, df_sorted_ans$relative_to_start + 1)
  x_max <- max(df_sorted_dmso$relative_to_start + 1, df_sorted_ans$relative_to_start + 1)
  
  log_range <- floor(log10(x_min)):ceiling(log10(x_max))
  custom_breaks <- unlist(lapply(log_range, function(p) seq(10^p, 10^(p + 1), length.out = 10)))
  
  if (grepl("PF69", i)) {
    color_nt <- "#909090"; color_treated <- "#F1811E"
  } else if (grepl("PF79", i)) {
    color_nt <- "#909090"; color_treated <- "#D4007C"
  } else if (grepl("PF80", i)) {
    color_nt <- "#909090"; color_treated <- "#00ADEB"
  } else if (grepl("PD63", i)) {
    color_nt <- "#000000"; color_treated <- "#FF000D"
  } else if (grepl("PD79", i)) {
    color_nt <- "#000000"; color_treated <- "#F000FF"
  } else if (grepl("PD80", i)) {
    color_nt <- "#000000"; color_treated <- "#0037FF"
  } else {
    color_nt <- "#000000"; color_treated <- "#FF0000"
  }
  
  plot <- ggplot() +
    geom_point(data = df_sorted_dmso, aes(x = (relative_to_start + 1), y = cum_expr), color = color_nt) +
    geom_point(data = df_sorted_ans,  aes(x = (relative_to_start + 1), y = cum_expr), color = color_treated) +
    scale_x_log10(breaks = custom_breaks, labels = scales::label_number(accuracy = 1)) +
    ggtitle(paste0(i, "_cumul_plot")) +
    theme_classic() +
    theme(legend.position = "none")
  
  assign(paste0(i, "_plot"), plot)
  ggsave(plot, file = file.path(figures_dir, paste0(i, "_cumulplot.pdf")), width = 4, height = 10, units = "in")
}