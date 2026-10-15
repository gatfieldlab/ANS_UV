library(gridExtra)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(ympes)
library(readxl)
library(ggrepel)
library(ympes)
library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)
library(tidyverse)
library(GO.db)
library(ggh4x)
library(FactoMineR)
library(factoextra)


save_path="~/Path/to/wd/"
file_path="~/path/to/file/"
setwd("~/path/to/wd/")



mapping <- c(
  "X18723" = "LD1",
  "X18724" = "LD2",
  "X18725" = "LD3",
  "X18726" = "LDMSO1",
  "X18727" = "LDMSO2",
  "X18728" = "LDMSO3",
  "X20128" = "HD1",
  "X20129" = "HD2"	,
  "X20130" = "HD3",
  "X20125" = "HDMSO1",
  "X20126" = "HDMSO2",
  "X20127" = "HDMSO3",
  "X20119" = "NTUV1",
  "X20120" = "NTUV2",
  "X20121" = "NTUV3",
  "X20122" = "UV1",
  "X20123" = "UV2",
  "X20124" = "UV3",

  "X15709" = "exclude",  #"Disome"
  "X15721" = "exclude",  #"Disome"
  "X15733" = "exclude" ,  #"Disome"
  "X15715" = "exclude",  #"Disome-ANS"
  "X15727" = "exclude",  #"Disome-ANS"
  "X15739" = "exclude",  #"Disome-ANS"
  "X15718" = "Di_LD1",  #"Disome-ANS-Mnase"
  "X15730" = "Di_LD2",
  "X15742" = "Di_LD3" ,
  "X15712" = "Di_LDMSO1",  #"Disome-Mnase"
  "X15724" = "Di_LDMSO2",
  "X15736" = "Di_LDMSO3",
  "X15716" = "mono_LD1",  #"Mono-ANS-Mnase"
  "X15728" = "mono_LD2",
  "X15740" = "mono_LD3" ,
  "X15707" = "exclude",  #"Monosome"
  "X15719" = "exclude",
  "X15731" = "exclude",
  "X15713" = "exclude",  #"Monosome-ANS"
  "X15725" = "exclude",
  "X15737" = "exclude",
  "X15710" = "mono_LDMSO1",  #"Monosome-Mnase"
  "X15722" = "mono_LDMSO2",
  "X15734" = "mono_LDMSO3",

  "X20152" = "Di_HD1",  #"ANS_HD_disome"
  "X20153" = "Di_HD2",
  "X20154"	= "Di_HD3",
  "X20140"	= "mono_HD1",  #"ANS_HD_mono"
  "X20141"	= "mono_HD2",
  "X20142"	= "mono_HD3",
  "X20149"	= "Di_HDMSO1",  #DMSO_disome
  "X20150"	= "Di_HDMSO2",
  "X20151"	= "Di_HDMSO3",
  "X20137"	= "mono_HDMSO1",  #DMSO_mono
  "X20138"	= "mono_HDMSO2",
  "X20139"	= "mono_HDMSO3",
  "X20143"	= "Di_NTUV1",  #NT_disome
  "X20144"	= "Di_NTUV2",
  "X20145"	= "Di_NTUV3",
  "X20131"	= "mono_NTUV1",  #NT_mono
  "X20132"	= "mono_NTUV2",
  "X20133"	= "mono_NTUV3",
  "X20146"	= "Di_UV1",  #UV_disome
  "X20147"	= "Di_UV2",
  "X20148"	= "Di_UV3",
  "X20134"	= "mono_UV1",  #UV_mono
  "X20135"	= "mono_UV2",
  "X20136" = "mono_UV3"
)

input_df_LD <- read.delim(paste0(file_path, "Table8_LD_input.txt"))
input_df_HD <- read.delim(paste0(file_path, "Table8_HD_input.txt"))
input_df_UV <- read.delim(paste0(file_path, "Table8_UV_input.txt"))

ribo_df_LD <- read.delim(paste0(file_path, "Table8_LD_ribo.txt"))
ribo_df_HD <- read.delim(paste0(file_path, "Table8_HD_ribo.txt"))
ribo_df_UV <- read.delim(paste0(file_path, "Table8_UV_ribo.txt"))


#keep one PG.genes per gene based on the max quantity log2
for(i in ls(pattern = "_df_")){
  
  tmp_file <- get(i)
  
  tmp_file <- tmp_file %>%
    rowwise() %>%
    mutate(total_signal = sum(c_across(ends_with(".Quantity.log2")), na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(PG.Genes) %>%
    slice_max(total_signal, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  assign(i, tmp_file)
}

#rename columns, optional - if column is renamed exclude -> exclude
for(i in ls(pattern = "_df_")) {
  tmp_file <- get(i)
  for (old in names(mapping)) {
    colnames(tmp_file) <- sub(paste0("^", old), mapping[old], colnames(tmp_file))

  }
  
  tmp_file <- tmp_file %>% dplyr::select(., !starts_with(c("exclude", "X")))
 
  #keep columns and conditions of interest
  tmp_file <- tmp_file %>% dplyr::select(., !contains(c("NrOfPrecursorsIdentified", 
                                                        "OrganismId",
                                                        "Organisms",
                                                        "Database",
                                                        "FastaFiles",
                                                        "IBAQ",
                                                        "NrOfPrecursorsUsedForQuantification",
                                                        "MolecularWeight",
                                                        "UniProtIds",
                                                        "ProteinGroups"
  )))
  
  assign(i, tmp_file)
  
}


#keep only the genes in input+ribo/disome (to avoid inf or NA when norm to input)
for(i in c("UV", "HD", "LD")){  
  
  input_df <- get(ls(pattern = paste0("input_df_",i)))
  input_df <-   dplyr::select(input_df, !total_signal)
  ribo_df <-  get(ls(pattern = paste0("ribo_df_",i)))
  ribo_df <-   dplyr::select(ribo_df, !total_signal)
  
  intersect_genes <- intersect(input_df$PG.Genes,ribo_df$PG.Genes ) #extract the proteins in common
  
  
  ribo_df <- subset(ribo_df, PG.Genes %in% intersect_genes)  #keep protein present in every datasets -> avoid NA of inf
  input_df <- subset(input_df, PG.Genes %in% intersect_genes) #keep prot present in input & ribo/disomes
  names(input_df) <- sub("\\.Quantity\\.log2$", "", names(input_df))
  
  assign(paste0("input_df_",i), input_df)
  assign(paste0("ribo_df_",i), ribo_df)
  
}



#build the function to normalize to input
normalize_by_input <- function(main_df, input_df) {
  for (col in names(main_df)) {
    #keep columns containing values only
    if (grepl("^PG\\.", col) || col == "prot_name") next
    
    #extract group id
    group_id <- sub(".*(LD\\d+|HD\\d+|HDMSO\\d+|LDMSO\\d+|UV\\d+|NTUV\\d+).*", "\\1", col)
    
    if (group_id %in% names(input_df)) {
      main_df[[paste0(col, "_norm")]] <- as.numeric(main_df[[col]]) - as.numeric(input_df[[group_id]])
    }
  }
  return(main_df)
}

#apply first normalization to input
for(i in c("UV", "HD", "LD")){
  
  input_df <- get(ls(pattern = paste0("input_df_",i)))
  ribo_df <-  get(ls(pattern = paste0("ribo_df_",i)))
  
  
  norm_ribo_df <- normalize_by_input(ribo_df, input_df) %>%
    dplyr::select(.,contains(c("PG.Genes", "norm"))) %>%
    mutate(prot_name=str_split_fixed(PG.Genes, ";", 2)[,1]) %>%
    dplyr::select(.,contains(c("PG.Genes", "prot_name", "norm")))
  
  assign(paste0("norm_ribo_df_",i), norm_ribo_df)
  
  
  norm_ribo_df_nolog <- data.frame(prot_name=norm_ribo_df$prot_name, 
                                   2^norm_ribo_df[,3:length(names(norm_ribo_df))])
  
  assign(paste0("norm_ribo_df_nolog_",i), norm_ribo_df_nolog)
  
  rm(norm_ribo_df)
  rm(norm_ribo_df_nolog)
  
  
}


#import list of RPS/L (ribo prot list) to norm to RPS/L
rpsl_list <- read.csv(paste0(file_path, "mgi_name_of_ribo_prot."))

#norm to ribo and plot some QC
for( i in c("UV", "HD", "LD")){
  
  norm_ribo_df_nolog <- get(ls(pattern=paste0("norm_ribo_df_nolog_",i)))
  
  norm_ribo_df_log <- data.frame(prot_name=norm_ribo_df_nolog$prot_name, 
                                 log2(norm_ribo_df_nolog[,2:length(norm_ribo_df_nolog)]))
  
  boxplot_nonorm <- norm_ribo_df_log %>%  
    pivot_longer(cols = contains("Quantity"), values_to = "counts", names_to = "sample")
  
  
  boxplot_notnorm_plot <- ggplot()+
    geom_boxplot(data=boxplot_nonorm, aes(x=sample, y=counts))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  ggsave(paste0(save_path, i,"_boxplot_notnorm_plot.pdf"), boxplot_notnorm_plot, ,width = 30, height = 16, units = "cm")
  
  
  rpsl_df <- subset(norm_ribo_df_nolog, prot_name %in% rpsl_list) 
  
  rpsl_df_for_plot <- subset(norm_ribo_df_nolog, prot_name %in% rpsl_list) %>%
    pivot_longer(., col=contains("norm"), values_to = "quantity", names_to = "sample")
  
  
  rpls_plot <- ggplot()+
    geom_point(data=rpsl_df_for_plot, aes(x=sample, y=log2(quantity), color=prot_name))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  ggsave(paste0(save_path, i,"_rpls_plot.pdf"), rpls_plot, ,width = 30, height = 16, units = "cm")
  
  median_rpsl_df <- apply(dplyr::select(rpsl_df, contains("norm")), 2, median) #norm to rpsl median
  
  median_for_plot <- median_rpsl_df %>% data.frame(.) %>% rownames_to_column("cdt")
  # pivot_longer(., col=contains("norm"), values_to = "quantity", names_to = "sample")
  
  
  median_plot <- ggplot()+
    geom_point(data=median_for_plot, aes(x=cdt, y=log2(.)))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  
  norm_ribo_df_nolog <- data.frame(prot_name=norm_ribo_df_nolog$prot_name, 
                                   sweep(dplyr::select(norm_ribo_df_nolog, contains("norm")), 2, median_rpsl_df, FUN = "/"))
  
  
  
  rpsl_df <- subset(norm_ribo_df_nolog, prot_name %in% rpsl_list) 
  
  
  rpsl_df_for_plot <- subset(norm_ribo_df_nolog, prot_name %in% rpsl_list) %>%
    pivot_longer(., col=contains("norm"), values_to = "quantity", names_to = "sample")
  
  
  rpsl_df_plot <- rpsl_df_for_plot
  rpsl_df_plot$group <- NA
  rpsl_df_plot[grep("Rpl",rpsl_df_plot$prot_name),]$group <- "Rpsl"
  # rpsl_df_plot[grep("Rpl",rpsl_df_plot$prot_name),]$group <- "Rpl"
  
  rpsl_df_plot[grep("Rps",rpsl_df_plot$prot_name),]$group <- "Rpsl"
  # rpsl_df_plot[grep("Rps",rpsl_df_plot$prot_name),]$group <- "Rps"
  
  rpls_plot_postnorm <- ggplot() +
    geom_boxplot(data = rpsl_df_plot, aes(x = sub("\\..*", "", sample), y = log2(quantity), color = group)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "red"))
  
  print(rpls_plot_postnorm)
  
  assign(paste0("rpsl_df_plot",i), rpsl_df_plot)
  ggsave(paste0(save_path, i,"_rpls_plot_postnorm.pdf"), rpls_plot_postnorm, ,width = 30, height = 16, units = "cm")

  norm_ribo_df_log2 <- data.frame(prot_name=norm_ribo_df_nolog$prot_name, 
                                  log2(dplyr::select(norm_ribo_df_nolog, contains("norm")))) 
  
  
  assign(paste0("norm_ribo_df_log2_",i), norm_ribo_df_log2)
  
  
  boxplot_norm <- norm_ribo_df_log2 %>%  
    pivot_longer(cols = contains("Quantity"), values_to = "counts", names_to = "sample")
  
  
  
  boxplot_norm_plot <- ggplot()+
    geom_boxplot(data=boxplot_norm, aes(x=sample, y=counts))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  ggsave(paste0(save_path, i,"boxplot_norm_plot_in_percent.pdf"), boxplot_norm_plot, ,width = 30, height = 16, units = "cm")
  
  rm(rpsl_df_plot)
}

ls(pattern = "norm_ribo_df_log2_")
ls(pattern = "rpsl_df_plot")

rpsl_df_list <- bind_rows(mget(ls(pattern = "rpsl_df_plot")))


mean_expression <- rpsl_df_list %>%
  mutate(condition = str_remove(sample, "\\d?\\..*")) %>%
  group_by(prot_name, condition, group) %>%
  summarise(mean_sample = mean(quantity))

#additional QC plots
rpls_plot_postnorm <- ggplot() +
  geom_boxplot(data = mean_expression, aes(x = condition, y = log2(mean_sample), color = group)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("darkblue", "red"))+
  scale_x_discrete(limits=c("mono_LDMSO",
                            "mono_HDMSO", 
                            "mono_NTUV",
                            "mono_LD",
                            "mono_HD", 
                            "mono_UV",
                            "Di_LDMSO",
                            "Di_HDMSO", 
                            "Di_NTUV",
                            "Di_LD",
                            "Di_HD", 
                            "Di_UV"))

print(rpls_plot_postnorm)

ggsave(paste0(save_path, "rpls_plot_postnorm.pdf"), rpls_plot_postnorm, ,width = 20, height = 10, units = "cm")


#boxplot RPS/L prot quantity
rpls_plot_postnorm_individuals <- ggplot() +
  geom_boxplot(data = rpsl_df_list, aes(x = sample, y = log2(quantity), color = group)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("darkblue", "red"))

ggsave(paste0(save_path, "rpls_plot_postnorm_individuals.pdf"), rpls_plot_postnorm_individuals, ,width = 20, height = 10, units = "cm")


ls(pattern = "norm_ribo_df_log2_")

#prepare table for plot
for(i in ls(pattern = "norm_ribo_df_log2_")){
  
  tmp_file <- get(i)
  tmp_file <- tmp_file %>%
    pivot_longer(
      cols = -prot_name, 
      names_to = "sample", 
      values_to = "quantity"
    ) %>%
    mutate(condition = str_remove(sample, "\\d?\\..*")) %>%
    group_by(prot_name, condition) %>%
    summarise(mean_sample = mean(quantity, na.rm = TRUE), .groups = "drop")
  
  assign(paste0(i, "long_table"), tmp_file)
  
  
}

all_prot_df <- bind_rows(mget(ls(pattern = "long_table$")))

#violin plot all prot quantity (+median)
violin_all_prot <- ggplot() +
  geom_violin(data = all_prot_df, aes(x = condition, y = mean_sample), draw_quantiles = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  # scale_color_manual(values = c("darkblue", "red"))+
  scale_x_discrete(limits=c("mono_LDMSO",
                            "mono_HDMSO", 
                            "mono_NTUV",
                            "mono_LD",
                            "mono_HD", 
                            "mono_UV",
                            "Di_LDMSO",
                            "Di_HDMSO", 
                            "Di_NTUV",
                            "Di_LD",
                            "Di_HD", 
                            "Di_UV"))

ggsave(paste0(save_path, "violin_all_prot.pdf"), violin_all_prot, ,width = 20, height = 10, units = "cm")

condition_limits <- c("mono_LDMSO", "mono_HDMSO", "mono_NTUV", 
                      "mono_LD", "mono_HD", "mono_UV",
                      "Di_LDMSO", "Di_HDMSO", "Di_NTUV", 
                      "Di_LD", "Di_HD", "Di_UV")

#figure XXX
merged_plot <- ggplot() +
  # 1. Violins en fond
  geom_violin(data = all_prot_df, aes(x = condition, y = mean_sample), draw_quantiles = 0.5) +
  
  # 2. Boxplot au premier plan (largeur réduite)
  geom_boxplot(data = mean_expression, aes(x = condition, y = log2(mean_sample), color = group)) +
  
  # Thème et échelles
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("darkblue", "red")) +
  scale_x_discrete(limits = condition_limits) +
  labs(
    y = "log2(mean_sample)",
    x = "Condition",
    color = "Groupe"
  )

ggsave(paste0(save_path, "violin_all_prot_merged.pdf"), merged_plot, ,width = 20, height = 10, units = "cm")


#prepare table for t-tests
for(i in c("HD", "UV", "LD")){
  
  norm_ribo_df_log2 <- get(ls(pattern = paste0("norm_ribo_df_log2_", i, "$")))
  
  long_df <- norm_ribo_df_log2 %>%
    pivot_longer( cols = contains("log2_norm"),
                  # !c("prot_name"),
                  values_to = "norm_counts", names_to = "sample") %>%
    mutate(
      condition = sub("[0-9]$", "", str_split_fixed(sample, "[.]", 2)[, 1]),
      replicate = str_extract(str_split_fixed(sample, "[.]", 2)[, 1], "[0-9]$")
    )
  
  assign(paste0("long_df_",i), long_df)
  rm(long_df)
}


long_df <- rbind(long_df_HD, long_df_LD, long_df_UV) 

#t-tests and log2FC per protein
t_test_results <- long_df %>%
  group_by(prot_name) %>%
  summarise(
    # Di_LD vs Di_LDMSO
    ttest_Di_LD_vs_DMSO = {
      vals1 <- norm_counts[condition == "Di_LD"]
      vals2 <- norm_counts[condition == "Di_LDMSO"]
      if (length(vals1) > 1 && length(vals2) > 1) t.test(vals1, vals2)$p.value else NA_real_
    },
    log2FC_Di_LD_vs_DMSO = {
      vals1 <- norm_counts[condition == "Di_LD"]
      vals2 <- norm_counts[condition == "Di_LDMSO"]
      if (length(vals1) > 0 && length(vals2) > 0) mean(vals1) - mean(vals2) else NA_real_
    },
    
    # mono_LD vs mono_LDMSO
    ttest_mono_LD_vs_DMSO = {
      vals1 <- norm_counts[condition == "mono_LD"]
      vals2 <- norm_counts[condition == "mono_LDMSO"]
      if (length(vals1) > 1 && length(vals2) > 1) t.test(vals1, vals2)$p.value else NA_real_
    },
    log2FC_mono_LD_vs_DMSO = {
      vals1 <- norm_counts[condition == "mono_LD"]
      vals2 <- norm_counts[condition == "mono_LDMSO"]
      if (length(vals1) > 0 && length(vals2) > 0) mean(vals1) - mean(vals2) else NA_real_
    },
    
    # mono_HD vs mono_HDMSO
    ttest_mono_HD_vs_DMSO = {
      vals1 <- norm_counts[condition == "mono_HD"]
      vals2 <- norm_counts[condition == "mono_HDMSO"]
      if (length(vals1) > 1 && length(vals2) > 1) t.test(vals1, vals2)$p.value else NA_real_
    },
    log2FC_mono_HD_vs_DMSO = {
      vals1 <- norm_counts[condition == "mono_HD"]
      vals2 <- norm_counts[condition == "mono_HDMSO"]
      if (length(vals1) > 0 && length(vals2) > 0) mean(vals1) - mean(vals2) else NA_real_
    },
    
    # Di_HD vs Di_HDMSO
    ttest_Di_HD_vs_DMSO = {
      vals1 <- norm_counts[condition == "Di_HD"]
      vals2 <- norm_counts[condition == "Di_HDMSO"]
      if (length(vals1) > 1 && length(vals2) > 1) t.test(vals1, vals2)$p.value else NA_real_
    },
    log2FC_Di_HD_vs_DMSO = {
      vals1 <- norm_counts[condition == "Di_HD"]
      vals2 <- norm_counts[condition == "Di_HDMSO"]
      if (length(vals1) > 0 && length(vals2) > 0) mean(vals1) - mean(vals2) else NA_real_
    },
    
    # mono_UV vs mono_NTUV
    ttest_mono_UV_vs_NTUV = {
      vals1 <- norm_counts[condition == "mono_UV"]
      vals2 <- norm_counts[condition == "mono_NTUV"]
      if (length(vals1) > 1 && length(vals2) > 1) t.test(vals1, vals2)$p.value else NA_real_
    },
    log2FC_mono_UV_vs_NTUV = {
      vals1 <- norm_counts[condition == "mono_UV"]
      vals2 <- norm_counts[condition == "mono_NTUV"]
      if (length(vals1) > 0 && length(vals2) > 0) mean(vals1) - mean(vals2) else NA_real_
    },
    
    # Di_UV vs Di_NTUV
    ttest_Di_UV_vs_NTUV = {
      vals1 <- norm_counts[condition == "Di_UV"]
      vals2 <- norm_counts[condition == "Di_NTUV"]
      if (length(vals1) > 1 && length(vals2) > 1) t.test(vals1, vals2)$p.value else NA_real_
    },
    log2FC_Di_UV_vs_NTUV = {
      vals1 <- norm_counts[condition == "Di_UV"]
      vals2 <- norm_counts[condition == "Di_NTUV"]
      if (length(vals1) > 0 && length(vals2) > 0) mean(vals1) - mean(vals2) else NA_real_
    },
    
    .groups = "drop"
  )


# join to replicate-level table
final_df <- long_df %>% dplyr::select(., -c(condition, replicate)) %>%
  pivot_wider(names_from = "sample", values_from = "norm_counts")%>%
  left_join(t_test_results, by = "prot_name")

write.csv(t_test_results, file=paste0(save_path, "t_test_results.csv"))
write.csv(final_df, file=paste0(save_path, "t_test_results_withvalues.csv"))

ttest_cols <- grep("^ttest_", names(final_df), value = TRUE)

#filter rows where any ttest column <= 0.05
res <- final_df$prot_name[
  apply(final_df[ , ttest_cols], 1, function(x) any(x <= 0.05, na.rm = TRUE))
]

# unique set of protein names
unique_proteins <- unique(res) #extract all prot sign at least once

write.csv(unique_proteins, file=paste0(save_path, "unique_proteins.csv"))

#prepare figure XX 

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

sub_df <- final_df
sub_df <- subset(sub_df, prot_name %in% unique_proteins)

all_logfc_cols <- grep("^log2FC_", names(sub_df), value = TRUE)

global_min <- min(unlist(sub_df[all_logfc_cols]), na.rm = TRUE)
global_max <- max(unlist(sub_df[all_logfc_cols]), na.rm = TRUE)

#set symmetric scale around 0
max_abs <- max(abs(c(global_min, global_max)))
global_min <- -max_abs
global_max <-  max_abs

#import the categories to label : table XXX
list_prot_to_lab <- read_excel("~/path/to/list_prot_to_lab.xlsx", 
                               col_names = FALSE)

names(list_prot_to_lab) <- c("prot", "type")

head(list_prot_to_lab)

#loop to plot per category (~GO category)
for(i in unique(list_prot_to_lab$type)){

  y="mono_di"
  
  lab_table <- list_prot_to_lab[list_prot_to_lab$type==i,]
  
  tmp_df <- final_df
  tmp_df <- subset(tmp_df, prot_name %in% lab_table$prot)

  tmp_df <- dplyr::select(tmp_df, prot_name, contains("ttest"), contains("log2FC"))

  
  #parameters to define
  alpha <- 0.05        #significance
  na_to_nonsig <- TRUE # treat NA p-values as non-significant
 
  
  df <- tmp_df
  
  # 1) grep t-test and log2FC columns
  ttest_cols  <- grep("^ttest_",  names(df), value = TRUE)
  logfc_cols  <- grep("^log2FC_", names(df), value = TRUE)
  
  contrast_key <- function(x) sub("^[^_]*_", "", x)
  
  #build data frame of matched pairs
  pairs <- tibble(
    ttest_col = ttest_cols,
    tkey      = contrast_key(ttest_cols)
  ) %>%
    inner_join(
      tibble(logfc_col = logfc_cols, lkey = contrast_key(logfc_cols)),
      by = c("tkey" = "lkey")
    ) %>%
    rename(contrast = tkey)
  
  #prepare df to long df (protein x contrast)
  long <- map_dfr(seq_len(nrow(pairs)), function(i) {
    tcol <- pairs$ttest_col[i]
    lcol <- pairs$logfc_col[i]
    con  <- pairs$contrast[i]
    
    tibble(
      prot_name = df$prot_name,
      contrast  = con,
      pval      = df[[tcol]],
      log2FC    = df[[lcol]]
    )
  })
  
  #NA p-values=non-significant
  if (na_to_nonsig) {
    long <- long %>% mutate(pval = if_else(is.na(pval), 1, pval))
  }
  
  #assign significance and categories
  long <- long %>%
    mutate(
      significant = !is.na(pval) & (pval < alpha),
      category = case_when(
        significant & !is.na(log2FC) & log2FC < 0 ~ "sig_down",
        significant & !is.na(log2FC) & log2FC > 0 ~ "sig_up",
        TRUE                                      ~ "nonsig"
      ),
      category = factor(category, levels = c("sig_down", "nonsig", "sig_up"))
    )
  
  #set colors and plot parameters
  grad_colors <- scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(global_min, global_max),  # <-- fixed across all GO terms
    oob = scales::squish,                # squish out-of-bounds to the edge colors
    na.value = "#BDBDBD",
    name = "log2FC"
  )
  
  dot_size <- 3
  
  p2 <- ggplot(long, aes(x = contrast, y = prot_name)) +
    geom_tile(aes(fill = log2FC), color = "white", size = 0.2) +
    grad_colors +
    geom_point(
      data = subset(long, significant),
      aes(x = contrast, y = prot_name),
      color = "black", size = dot_size
    ) +
    labs(
      x = NULL, y = NULL,
      title = paste0("log2FC Heatmap + Significance Dots: ", i),
      subtitle = sprintf("Dot = p < %.3f | Scale = [%.2f, %.2f]", alpha, global_min, global_max)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      legend.position = "top"
    )
  
  print(p2)
  

  ordering <- long %>%
    group_by(prot_name) %>%
    summarize(
      n_sig = sum(category %in% c("sig_down", "sig_up"), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # arrange((n_sig), desc(prot_name)) %>%   # sort by n_sig, then alphabetically
    arrange(desc(prot_name)) %>%   #sort alphabetically
    pull(prot_name)
  
  
  long_ord <- long %>%
    mutate(prot_name = factor(prot_name, levels = ordering))
  
  n_prots     <- length(unique(long_ord$prot_name))
  n_contrasts <- 6
  
  tile_size_cm <- 1.0
  margin_w_cm  <- 5
  margin_h_cm  <- 4
  
  plot_width  <- n_contrasts * tile_size_cm + margin_w_cm
  plot_height <- n_prots     * tile_size_cm + margin_h_cm
  
  p2_ord <- ggplot(long_ord, aes(x = contrast, y = prot_name)) +
    geom_tile(aes(fill = log2FC), color = "white", size = 0.2) +
    grad_colors +
    geom_point(
      data = subset(long_ord, significant),
      aes(x = contrast, y = prot_name),
      color = "black", size = dot_size
    ) +
    labs(
      x = NULL, y = NULL,
      title = paste0(i, "_", y),
      subtitle = sprintf("Dot = p < %.3f | Scale = [%.2f, %.2f]", alpha, global_min, global_max)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid     = element_blank(),
      legend.position = "top"
    ) +
    # scale_x_discrete(limits = c(
    #   paste0(y, "_LD_vs_DMSO"),
    #   paste0(y, "_HD_vs_DMSO"),
    #   paste0(y, "_UV_vs_NTUV")
    # )) +
    scale_x_discrete(limits = c(
      "mono_LD_vs_DMSO",
      "mono_HD_vs_DMSO",
      "mono_UV_vs_NTUV",
      "Di_LD_vs_DMSO",
      "Di_HD_vs_DMSO",
      "Di_UV_vs_NTUV"))+
    force_panelsizes(
      rows = unit(n_prots     * tile_size_cm, "cm"),
      cols = unit(n_contrasts * tile_size_cm, "cm")
    )
  
  print(p2_ord)
  
  z <- gsub(":", "_", i)
  ggsave(
    paste0(save_path, z, "_", y, "_", i, "_GO_heatmap.pdf"),
    p2_ord,
    width  = plot_width,
    height = plot_height,
    units  = "cm"
  )
}

#plot volcano contians proteins to color on volcano and their ~GO
list_for_figure <- read_excel("list_for_figure.txt")

#contains the GO for each volcano and the prot to label
list_to_color <- c("Chaperones & protein quality control", 
                   "Chromatin / transcription regulation", 
                   "cytoskeleton", 
                   "ER stress", 
                   "RQC", 
                   "Translation", 
                   "ubiquitin")

t_test_results$label <- "ns"

######### calculate mean expression, log2fc per replicate, ratio of log2fc


df <- long_df


df2 <- df %>%
  mutate(
    celltype = ifelse(str_detect(condition, "^Di"), "Di", "mono"),
    dose     = str_extract(condition, "LD|HD|UV"),
    treatment = case_when(
      str_detect(condition, "MSO")  ~ "control",
      str_detect(condition, "NTUV") ~ "control",
      TRUE                          ~ "treated"
    )
  )

df_wide <- df2 %>%
  dplyr::select(prot_name, celltype, dose, replicate, treatment, norm_counts) %>%
  pivot_wider(
    names_from  = treatment,
    values_from = norm_counts
  )

log2fc_rep <- df_wide %>%
  mutate(log2FC = treated - control)

log2fc_matrix <- log2fc_rep %>%
  unite(cond_rep, celltype, dose, replicate, sep = "_") %>%
  dplyr::select(prot_name, cond_rep, log2FC) %>%
  pivot_wider(
    names_from  = cond_rep,
    values_from = log2FC
  )

write.csv(log2fc_matrix, file=paste0(save_path, "log2fc_per_replicate.csv"))

log2fc_rep_wide <- log2fc_rep %>%
  pivot_wider(
    id_cols    = c(prot_name, dose, replicate),
    names_from = celltype,
    values_from = log2FC
  )

log2fc_di_vs_mono_matrix <- log2fc_rep_wide %>%
  mutate(log2FC_Di_vs_mono = Di - mono) %>%
  dplyr::select(prot_name, dose, replicate, log2FC_Di_vs_mono)%>%
  unite(cond_rep, dose, replicate, sep = "_") %>%
  pivot_wider(
    names_from  = cond_rep,
    values_from = log2FC_Di_vs_mono
  )


write.csv(log2fc_di_vs_mono_matrix, file=paste0(save_path, "log2fc_rep_diff.csv"))

df <- long_df

log2_mean_expression <- df %>% group_by(prot_name, condition)%>%
  summarise(mean_expression=mean(norm_counts)) %>%
  pivot_wider(names_from = condition, values_from = mean_expression)

write.csv(log2_mean_expression, file=paste0(save_path, "log2_mean_expression.csv"))


############
i="RQC"
for (i in list_to_color){
  t_test_results$label[t_test_results$prot_name %in% list_for_figure[list_for_figure$category==i,]$prot] <- i
  
}

t_test_results$label[t_test_results$prot_name=="Edf1"] <- "Edf1"
t_test_results$label[grepl("eif", t_test_results$prot_name, ignore.case = TRUE)] <- "initiation_factors"


y="mono_HD"
i="RQC"

rm(tmp_file)
i="Di_LD"
list_to_loop <- c("Di_LD","mono_LD","Di_HD","mono_HD","Di_UV","mono_UV" )
for(y in list_to_loop){

  
  tmp_file <- dplyr::select(t_test_results, prot_name, contains(y))%>%
    drop_na()
  
  tmp_file$label <- "ns"
  tmp_file[tmp_file[,2]<=0.05,]$label <- "sig"
 
  log2_col  <- grep("^log2FC_", names(tmp_file), value = TRUE)[1]
  ttest_col <- grep("^ttest_",  names(tmp_file), value = TRUE)[1]
  
  total_prot <- nrow(tmp_file)
  
  sig_up <- sum(tmp_file[[log2_col]]  > 0 &
                  tmp_file[[ttest_col]] <= 0.05,
                na.rm = TRUE)
  
  sig_down <- sum(tmp_file[[log2_col]]  < 0 &
                    tmp_file[[ttest_col]] <= 0.05,
                  na.rm = TRUE)
  
  subtitle_text <- paste0(
    "Total proteins: ", total_prot,
    " | significant up: ", sig_up,
    " | significant down: ", sig_down
  )
  
  
  targets_RQC <- c("Edf1", "Znf598", "Map3k20", "Serbp1", "Ascc1", "Trip4",
                   "Hbs1l", "Pelo", "Ltn1", "Nemf", "Ankzf1", "Ufl1", "Ddrgk1", "Cdk5rap3", "Ufm1", "Vcp", "Ube3c", 
                   "Eif5a", "Abcf1")
  targets_elongation <- c("Eef1b", "Eef1d", "Eef1g", "Eef2", "Hgh1", "Pa2g4")
  targets_termination <- c("Etf1", "Abce1")
  targets_eif1 <- c("Eif1", "Eif1a", "Eif1ax")
  targets_initation <- c("Eif2b2", "Eif2b3", "Eif2b4", "Eif2d", "Eif2s1", 
                         "Eif2s2", "Eif3a", "Eif3b", "Eif3c", "Eif3d", "Eif3e", "Eif3f", "Eif3g", "Eif3h", "Eif3i", 
                         "Eif3j1", "Eif3k", "Eif3l", "Eif3m", "Eif4a1", "Eif4b", "Eif4g1", "Eif4g2", "Eif5", "Eif5b", "Eif6")
  
  
  
  col_x <- colnames(tmp_file)[3]
  col_y <- colnames(tmp_file)[2]
  
  Volcano_MS <- ggplot(tmp_file, aes(x = .data[[col_x]], y = -log10(.data[[col_y]]))) +
    geom_point(aes(color = label), size = 1.5) +
        geom_point(data = subset(tmp_file, prot_name %in% targets_RQC), 
               color = "cornflowerblue", size = 1.5) +
    geom_text_repel(data = subset(tmp_file, prot_name %in% targets_RQC),
                    aes(label = prot_name),
                    color = "cornflowerblue", fontface = "bold",
                    box.padding = 0.5,
                    max.overlaps = 20) +
    geom_point(data = subset(tmp_file, prot_name %in% targets_eif1), 
               color = "purple2", size = 1.5) +
        geom_text_repel(data = subset(tmp_file, prot_name %in% targets_eif1),
                    aes(label = prot_name),
                    color = "purple2", fontface = "bold",
                    box.padding = 0.5,
                    max.overlaps = 20) +
    geom_point(data = subset(tmp_file, prot_name %in% targets_initation), 
               color = "purple4", size = 1.5) +
    geom_point(data = subset(tmp_file, prot_name %in% targets_elongation), 
               color = "blue", size = 1.5) +
    geom_point(data = subset(tmp_file, prot_name %in% targets_termination), 
               color = "red", size = 1.5) +
    theme_classic() +
    scale_color_manual(values = c("grey", "darkblue")) +
    scale_x_continuous(breaks = seq(-5.5, 5, 0.5), limits = c(-5.5, 5.5)) +
    scale_y_continuous(breaks = seq(0, 6, 1), limits = c(0, 5)) +
    theme(legend.position = "none")+
    xlab("Log2FC(treatment vs control)")+
    ylab("-log10(pval)")
  
  
  
  x=1

  ggsave(Volcano_MS, file = paste0(save_path, y,"_volcano_MS.pdf"), width = 6, height = 5, units = "in")
  
  assign(paste0(y,"_volcano_MS"), Volcano_MS)
  
  rm(tmp_file)
  
  tmp_file <- dplyr::select(t_test_results, prot_name, contains(y))%>%
    drop_na()
  
  tmp_file$label <- "ns"
  tmp_file[tmp_file[,2]<=0.05,]$label <- "sig"
  tmp_file$label[grepl("eif", tmp_file$prot_name, ignore.case = TRUE)] <- "initiation_factors"
  
  tmp_file_small <- subset(tmp_file, label == "initiation_factors")

  
}



ggplot()+
  geom_point(data=t_test_results,  aes(x = log2FC_Di_LD_vs_DMSO, y = -log10(ttest_Di_LD_vs_DMSO),  color = label))

head(t_test_results)

select_col <- str_split_fixed(grep(("ttest"), names(t_test_results), value=TRUE), "_", 2)[,2]

#GO analysis
for(i in select_col){
  
  col_select <- grep(i,names(t_test_results), value=TRUE)
  
  
  up_prot <- t_test_results[
    t_test_results[[col_select[1]]] < 0.05 & t_test_results[[col_select[2]]] > 0,
    "prot_name"
  ]
  
  
  down_prot <- t_test_results[
    t_test_results[[col_select[1]]] < 0.05 & t_test_results[[col_select[2]]] < 0,
    "prot_name"
  ]
  
  
  assign(paste0(i, "_up_prot"), up_prot)
  assign(paste0(i, "_down_prot"), down_prot)
  
  
  GO_up<-enrichGO(gene=unique(up_prot$prot_name),
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  # ont="BP",
                  ont="BP", # ont="MF" is the default
                  minGSSize = 30)
  
  GO_up <- simplify(GO_up, cutoff = 0.5, by = "p.adjust", select_fun = min)
  
  GO_down<-enrichGO(gene=unique(down_prot$prot_name),
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    # ont="BP",
                    ont="BP", # ont="MF" is the default
                    minGSSize = 30)
  
  GO_down <- simplify(GO_down, cutoff = 0.5, by = "p.adjust", select_fun = min)
  
  assign(paste0(i, "_GO_enrich_up"), GO_up@result)
  assign(paste0(i, "_GO_enrich_down"), GO_down@result)
  
  GO_plot_dot_up <- dotplot(GO_up, showCategory=10)+ggtitle(paste0(i, "_up"))
  GO_plot_dot_down <- dotplot(GO_down, showCategory=10)+ggtitle(paste0(i, "_down"))
  
  arranged_plots <- grid.arrange(GO_plot_dot_up, GO_plot_dot_down, nrow=1)
  
  ggsave(
    arranged_plots,
    file = paste0("~/BTF3 study/20250624_lib_ANSHD_UV30/figure/", i, "GO_plot_dot.pdf"),
    width = 12,
    height = 8,
    units = "in"
  )
  
}

#Plot PCA
for(i in ls(pattern = "ribo_df_log2_")){
  df <- get(i)
  
  intensity_cols <- grep("Quantity.log2_norm$", colnames(df), value = TRUE)
  
  cond_cols <- grep(cond, intensity_cols, value = TRUE)
  
  mat <- df[, intensity_cols, drop = FALSE]
  
  n_proteins <- nrow(mat)
  
  mat_t <- t(mat)
  
  pca_res <- prcomp(mat_t, scale. = TRUE)
  
  pca_df <- data.frame(
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2],
    sample = rownames(mat_t),
    condition = sub("[0-9]$","",str_split_fixed(intensity_cols,  "[.]",2)[,1])
  )
  
  p <- ggplot(pca_df, aes(PC1, PC2, label = sample, color=condition)) +
    geom_point(size = 4) +
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    geom_text_repel() +
    theme_minimal() +
    labs(
      title = paste("PCA for condition:", i),
      subtitle = paste0("n = ", n_proteins, " protéines"), # Ajout du nombre de protéines ici
      x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)")
    )
  
  print(p)
  
  ggsave(
    p,
    file = paste0("~/Path/to/figure/", i, "_PCA_MS.pdf"),
    width = 8,
    height = 6,
    units = "in"
  )
}

conditions <- c("mono",
                "Di")

#plot PCA
for (cond in conditions) {
  
  message("Running PCA for condition: ", cond)
  
  cond_cols <- grep(cond, intensity_cols, value = TRUE)
  
  mat <- df[, cond_cols, drop = FALSE]
  
  mat <- mat %>% mutate_all(~ ifelse(is.na(.), mean(., na.rm = TRUE), .))
  
  mat_t <- t(mat)
  
  pca_res <- prcomp(mat_t, scale. = TRUE)
  
  pca_df <- data.frame(
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2],
    sample = rownames(mat_t),
    condition = sub("[0-9]$","",str_split_fixed(cond_cols$sample,  "[.]",2)[,1])
  )
  
  
  p <- ggplot(pca_df, aes(PC1, PC2, label = sample, color=condition)) +
    geom_point(size = 4) +
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    geom_text_repel() +
    theme_minimal() +
    labs(
      title = paste("PCA for condition:", cond),
      x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)")
    )
  
  print(p)
  
  ggsave(
    p,
    file = paste0("~/path/to/figure/", cond, "_PCA_MS.pdf"),
    width = 14,
    height = 12,
    units = "in"
  )
}





