library(DESeq2)
library(ggplot2)
library(data.table)
library(tidyverse)
library(DESeq2)
library(dplyr)
library(tibble)
library(matrixStats)
library(ggplot2)
library(ggrepel)


##sample labeling
#A1HD = rep x of HD treatment
#HDMSO = HD matching DMSO mock cond
#AxLD = rep x LD treated
#LDMSO = LD matching mock cond
#uxshort = UV cond treated rep x
#n

rna_counts <- as.data.frame(fread("~/path/to/Counts_cDNA.txt", sep='\t', header=T))


rna_counts <- rna_counts %>% select(-contains("long"))

#rename UV columns
names(rna_counts) <- gsub("UV", "u", names(rna_counts))

#keep most expressed in RNA one per gene 
rna_counts_collapsed <- rna_counts %>%
  mutate(total_expression = rowSums(across(A1HD:u4short))) %>%
  group_by(gene_id) %>%
  slice_max(order_by = total_expression, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-total_expression)

rownames(rna_counts_collapsed) <- NULL

rna_counts <- rna_counts_collapsed %>% column_to_rownames(., "gene_id") %>%
  select(-transcript_id)

sum_per_condition <- colSums(rna_counts)

names(sum_per_condition) <- gsub("^u", "au", names(sum_per_condition)) #to have test always on top (alphabetical order)

df <- data.frame(
  sample = names(sum_per_condition),
  expression = as.numeric(sum_per_condition)
)

#create the experiement conditions column
df$experiment <- ifelse(grepl("HD", df$sample), "High Dose",
                        ifelse(grepl("LD", df$sample), "Low Dose", "UV (short)"))

#def condition
df$condition <- ifelse(grepl("^a", df$sample, ignore.case = TRUE), "Treated", "Control")

df$condition <- factor(df$condition, levels = c("Treated", "Control"))

df <- df[order(df$experiment, df$condition, df$sample), ]

RNA_qc <- ggplot(df, aes(x = experiment, y = expression, fill = sample )) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sample), position = position_stack(vjust = 0.5), color = "white", size = 3, fontface = "bold") +
  theme_classic(base_size = 14) +
  labs(
    title = "Total counts RNA-Seq",
    x = "Experiment",
    y = "Total Expression Count",
    fill = "Condition"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30")
  )

ggsave(RNA_qc, file=paste0("~/BTF3 study/20240318_library ANS and UV/figures/RNA_qc_plot.pdf"),width = 7, height = 5, units = "in")




#######
#Volcanos
########


hd_cols <- c("D1HD", "D2HD", "D3HD", "D4HD", "A1HD", "A2HD", "A3HD", "A4HD")

counts_mat <- rna_counts
hd_counts <- counts_mat[, hd_cols]

hd_metadata <- data.frame(
  row.names = hd_cols,
  condition = factor(c(rep("control", 4), rep("treated", 4)), levels = c("control", "treated"))
)

#DESeq2
dds_hd <- DESeqDataSetFromMatrix(countData = hd_counts, colData = hd_metadata, design = ~ condition)
dds_hd <- DESeq(dds_hd)
res_hd <- results(dds_hd, contrast = c("condition", "treated", "control"))

hd_vst <- assay(vst(dds_hd, blind = TRUE))

top_500_hd_indices <- head(order(rowVars(hd_vst), decreasing = TRUE), 500)
hd_vst_top500 <- hd_vst[top_500_hd_indices, ]

#PCA on top 500 most variable
pca_result_HD_500 <- prcomp(t(hd_vst_top500), scale. = TRUE)

pca_coords_HD_500 <- as.data.frame(pca_result_HD_500$x) %>%
  rownames_to_column("Sample")

plot_HD <- ggplot(pca_coords_HD_500, aes(x = PC1, y = PC2, label = Sample, color = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, show.legend = FALSE) + 
  theme_minimal() +
  labs(
    title = "PCA (VST + Top 500 Gènes Variables) - Groupe HD",
    x = paste0("PC1 (", round(100 * summary(pca_result_HD_500)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_result_HD_500)$importance[2,2], 1), "%)"),
    color = "Échantillon"
  )

print(plot_HD)



ld_cols <- c("D1LD", "D2LD", "D3LD",  "A1LD", "A2LD", "A3LD")
ld_counts <- counts_mat[, ld_cols]

ld_metadata <- data.frame(
  row.names = ld_cols,
  condition = factor(c(rep("control", 3), rep("treated", 3)), levels = c("control", "treated"))
)

dds_ld <- DESeqDataSetFromMatrix(countData = ld_counts, colData = ld_metadata, design = ~ condition)
dds_ld <- DESeq(dds_ld)
res_ld <- results(dds_ld, contrast = c("condition", "treated", "control"))


#same for UV
uv_cols <- c("N1short", "N2short", "N3short", "N4short" ,"u1short", "u2short", "u3short", "u4short") 
uv_counts <- counts_mat[, uv_cols]

uv_metadata <- data.frame(
  row.names = uv_cols,
  condition = factor(c(rep("control", 4), rep("treated", 4)), levels = c("control", "treated"))
)

dds_uv <- DESeqDataSetFromMatrix(countData = uv_counts, colData = uv_metadata, design = ~ condition)
dds_uv <- DESeq(dds_uv) 
res_uv <- results(dds_uv, contrast = c("condition", "treated", "control"))


for(i in ls(pattern = "^res_")){
  
  result_data <- get(i)
  result_data <- as.data.frame(result_data)
  result_data$gene_id <- rownames(result_data)
  
  result_data <- result_data %>% filter(!is.na(padj) & !is.na(log2FoldChange))
  
  result_data <- result_data %>%
    mutate(status = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up-regulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
      TRUE ~ "Not Significant"
    ))
  
  result_data$status <- factor(result_data$status, levels = c("Down-regulated", "Not Significant", "Up-regulated"))
  
  n_up    <- sum(result_data$status == "Up-regulated")
  n_down  <- sum(result_data$status == "Down-regulated")
  n_tot <- nrow(result_data)
  
  subtitle_text <- paste0("n_tot=", n_tot, ", n_up=", n_up, ", n_down=", n_down)
  
  top_genes <- result_data %>%
    filter(status != "Not Significant") %>%
    arrange(padj) %>%
    head(10)
  
  colcanoplot <- ggplot(result_data, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
    geom_point() +
    
    scale_color_manual(values = c("Down-regulated" = "#377eb8", 
                                  "Not Significant" = "#999999", 
                                  "Up-regulated" = "#e41a1c")) +
    
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.4) +
    
    # geom_text_repel(data = top_ld_genes, aes(label = gene_id), 
    #                 size = 3.5, color = "black", fontface = "bold",
    #                 box.padding = 0.4, point.padding = 0.3, segment.color = "grey50") +
    
    # Inject the calculated counts into the subtitle
    labs(
      title = paste0("RNAseq: ", i),
      subtitle = subtitle_text,
      x = expression("log"[2]*" Fold Change"),
      y = expression("-log"[10]*" Adjusted P-value"),
      color = "Gene Status"
    ) +
    xlim(-5.5,5.5)+
    theme_classic()
  
  ggsave(colcanoplot, file=paste0("~/BTF3 study/20240318_library ANS and UV/figures/", i, "_volcano.pdf"),width = 7, height = 5, units = "in")
}




counts_hd <- counts(dds_hd)
counts_ld <- counts(dds_ld)
counts_uv <- counts(dds_uv)

common_genes <- reduce(list(rownames(counts_hd), rownames(counts_ld), rownames(counts_uv)), intersect)

counts_combined <- cbind(
  counts_hd[common_genes, ],
  counts_ld[common_genes, ],
  counts_uv[common_genes, ]
)

coldata_combined <- data.frame(
  row.names = colnames(counts_combined),
  Group = factor(c(
    rep("HD", ncol(counts_hd)),
    rep("LD", ncol(counts_ld)),
    rep("UV", ncol(counts_uv))
  ))
)

dds_combined <- DESeqDataSetFromMatrix(
  countData = counts_combined,
  colData = coldata_combined,
  design = ~ Group
)

combined_vst <- assay(vst(dds_combined, blind = TRUE))

top_500_indices <- head(order(rowVars(combined_vst), decreasing = TRUE), 500)
vst_top500 <- combined_vst[top_500_indices, ]

pca_result <- prcomp(t(vst_top500), scale. = TRUE)

pca_coords <- as.data.frame(pca_result$x) %>%
  rownames_to_column("Sample") %>%
  mutate(Group = coldata_combined[Sample, "Group"])

plot_square <- ggplot(pca_coords, aes(x = PC1, y = PC2, label = Sample, color = Group)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 15) + 
  scale_color_manual(values = c("HD" = "royalblue", "LD" = "forestgreen", "UV" = "firebrick")) +
  coord_fixed() + # Force un ratio 1:1 (carré)
  theme_minimal() +
  theme(
    legend.position = "none" # Supprime complètement la légende
  ) +
  labs(
    title = "PCA (VST + Top 500 Gènes Variables)",
    x = paste0("PC1 (", round(100 * summary(pca_result)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_result)$importance[2,2], 1), "%)")
  )

ggsave(file=paste0("~/path/to/figures/plot_combined.pdf"), plot = plot_square, width = 7, height = 7, dpi = 300)

