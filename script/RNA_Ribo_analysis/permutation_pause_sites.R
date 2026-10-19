library(stringr)
library(ggthemes)
library(ggplot2)
library(ggtext)
library(org.Mm.eg.db)
library(tidyverse)
library(gridExtra)

codon_map <- data.frame(
  codon = c(
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", 
    "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", 
    "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", 
    "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", 
    "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", 
    "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", 
    "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", 
    "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
  ),
  aa = c(
    "Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", 
    "Arg", "Ser", "Arg", "Ser", "Ile", "Ile", "Met", "Ile", 
    "Gln", "His", "Gln", "His", "Pro", "Pro", "Pro", "Pro", 
    "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu", 
    "Glu", "Asp", "Glu", "Asp", "Ala", "Ala", "Ala", "Ala", 
    "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val", "Val", 
    "Stp", "Tyr", "Stp", "Tyr", "Ser", "Ser", "Ser", "Ser", 
    "Stp", "Cys", "Trp", "Cys", "Leu", "Phe", "Leu", "Phe"
  ),
  stringsAsFactors = FALSE
)




import_file <- list.files("~/path/file/fisher_tests/", pattern = "vs")
list.files("~/path/fisher_tests/")

link <- "~/path/file/fisher/"
FisherTest_tx_Codoncontent <- read.table(paste0(link, "FisherTest_tx_Codoncontent.txt"), sep="\t", header = TRUE)


list_histone_to_remove <- read.csv("~/path/file/histone/exclusion.csv", header = TRUE)

Mmusculus.GRCm39.111.startStop.APPRIS <- read.delim("~/path/Mmusculus.GRCm39.111.startStop.IDS.bed", header=FALSE)
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[1]="merge"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[2]="start"
colnames(Mmusculus.GRCm39.111.startStop.APPRIS)[3]="stop"
Mmusculus.GRCm39.111.startStop.APPRIS<-Mmusculus.GRCm39.111.startStop.APPRIS[-4]
split_merge <-as.data.frame(str_split_fixed(Mmusculus.GRCm39.111.startStop.APPRIS$merge, "[|]", 2))

Mmusculus.GRCm39.111.startStop.APPRIS<-cbind(Mmusculus.GRCm39.111.startStop.APPRIS, split_merge)
names(Mmusculus.GRCm39.111.startStop.APPRIS) <- c("merge", "start", "stop", "gene_id", "transcript_id")

for (i in import_file){
  
  tmp_file <- read.csv(paste0(link, i), sep="\t", header=TRUE)
  name_file <- str_split_fixed(i, "[.]", 2)[1]
  
  # Remove histones
  tmp_file <- subset(tmp_file,
                     !(transcript_id %in% list_histone_to_remove$ensembl_transcript_id))
  
  # Join annotation (NO filtering on relative_to_start here)
  tmp_file <- list(tmp_file, Mmusculus.GRCm39.111.startStop.APPRIS) %>%
    purrr::reduce(left_join, by="transcript_id") %>%
    dplyr::select(-gene_id.y) %>%
    dplyr::rename(gene_id = gene_id.x) %>%
    mutate(relative_to_start = pos - start)
  
  # Remove problematic transcripts
  tmp_file <- subset(tmp_file, !(transcript_id %in% c(
    "ENSMUST00000038863", # LARS2
    "ENSMUST00000004389", # GRCC10
    "ENSMUST00000094434"  # Ftl1
  )))
  
  # Stats before filtering
  n_gene <- length(unique(tmp_file$mgi_symbol))
  n_site <- length(tmp_file$pos)
  vec_all <- c(n_gene, n_site)
  
  # Clean odds ratios
  tmp_file$oddRatios[is.na(tmp_file$oddRatios)] <- -10
  tmp_file <- tmp_file[tmp_file$oddRatios > 0, ]  # keep only up sites
  
  # Initialize direction
  tmp_file$direction <- "ns"
  
  # Define TOP 1000 "up" sites
  top_idx <- order(tmp_file$pValue)[1:1000] ###CHECK PVAL OR PADJ
  tmp_file$direction[top_idx] <- "up"
  
  ###############################################
  #  SAME GENE SET LOGIC
  ###############################################
  
  # 1. Initial sig sites
  sites_oddratios_sig <- tmp_file[top_idx, ]
  
  # 2. Restrict to genes with ≥1 sig site
  gene_sig <- unique(sites_oddratios_sig$gene_id)
  tmp_gene_restricted <- tmp_file[tmp_file$gene_id %in% gene_sig, ]
  
  setdiff(unique(tmp_gene_restricted$transcript_id),unique(sites_oddratios_sig$transcript_id) )
  
  # 3. Recompute sig inside restricted universe
  top_idx2 <- order(tmp_gene_restricted$pValue)[1:1000]
  sites_oddratios_sig <- tmp_gene_restricted[top_idx2, ]
  
  # 4. Define ns from same genes
  sites_oddratios_ns_sig <- tmp_gene_restricted[
    tmp_gene_restricted$padj > 0.05, ] ###CHECK IF PADJ IS OK
  # 5. Enforce identical gene set (iterate until stable)
  
  
  common_genes <- intersect(unique(sites_oddratios_sig$gene_id),
                            unique(sites_oddratios_ns_sig$gene_id))
  
  sites_oddratios_sig    <- sites_oddratios_sig[sites_oddratios_sig$gene_id %in% common_genes, ]
  sites_oddratios_ns_sig <- sites_oddratios_ns_sig[sites_oddratios_ns_sig$gene_id %in% common_genes, ]
  
  
  
  ###############################################
  # Stats
  ###############################################
  
  n_gene_sig <- length(unique(sites_oddratios_sig$mgi_symbol))
  n_site_sig <- length(sites_oddratios_sig$pos)
  
  n_gene_ns_sig <- length(unique(sites_oddratios_ns_sig$mgi_symbol))
  n_site_ns_sig <- length(sites_oddratios_ns_sig$pos)
  
  vec_sig <- c(n_gene_sig, n_site_sig)
  vec_ns_sig <- c(n_gene_ns_sig, n_site_ns_sig)
  
  ###############################################
  # Annotation
  ###############################################
  
  sites_oddratios_sig$group <- ifelse(
    sites_oddratios_sig$oddRatios > 0, "up", "down")
  
  sites_oddratios_ns_sig$group <- ifelse(
    sites_oddratios_ns_sig$oddRatios > 0, "up", "down")
  
  ###############################################
  # Save objects
  ###############################################
  
  assign(paste0(name_file, "_permut"), tmp_gene_restricted)
  
  assign(paste0(i, "_sig_data"), sites_oddratios_sig)
  assign(paste0(i, "_ns_sig_data"), sites_oddratios_ns_sig)  
  
  tmp_df <- data.frame(all = vec_all,
                       sig = vec_sig,
                       ns_sig = vec_ns_sig)
  row.names(tmp_df) <- c("ngene", "nsite")
  
  assign(paste0(name_file, "_stats"), tmp_df)
  
  assign(name_file, tmp_gene_restricted)
  assign(paste0(name_file, "_sig"), sites_oddratios_sig)
  assign(paste0(name_file, "_ns_sig"), sites_oddratios_ns_sig)
  
}



###############################################
### 0. Setup & assumptions
###############################################
# Assumes you already have a data.frame 'df' with columns:
#   - codon      : codon at the site (e.g. "AAA")
#   - direction  : factor/character with at least "up" and "ns" (non-significant)
#   - padj       : adjusted p-value (optional for this script; used only if you redefine ns)

# If you need to read df from file, uncomment and adapt:
# df <- read.table("your_table.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


all_log2fc_list <- list()

##############################################
### 1. Define codon universe
###############################################
for(y in ls(pattern = "^FisherTest.*(mono|di)$")){
  
  df<- get(y)
  
  # All codons present in the dataset
  codons <- sort(unique(df$codon))
  
  
  ###############################################
  ### 2. Permutation: test enrichment of "up" sites per codon
  ###############################################
  
  # Sanity check: direction should have at least some "up" and some non-up
  print(table(df$direction))
  
  # perm_one: shuffle 'direction' across rows and count "up" codons
  perm_one <- function(df, codons) {
    # Shuffle labels
    df$perm_dir <- sample(df$direction)
    
    # Count codons among permuted "up" sites
    perm_tab <- table(df$codon[df$perm_dir == "up"])
    
    # Ensure we always return counts for ALL codons (fill missing with 0)
    out <- setNames(rep(0, length(codons)), codons)
    out[names(perm_tab)] <- as.numeric(perm_tab)
    
    return(out)
  }
  
  # Number of permutations
  n_perm <- 10000   # you can reduce to e.g. 1000 for testing   #CHECK IF WE NEED TO REDUCE
  
  # Allocate permutation matrix: rows = codons, columns = permutations
  perm_mat <- matrix(
    0,
    nrow = length(codons),
    ncol = n_perm,
    dimnames = list(codons, NULL)
  )
  
  set.seed(123)  # for reproducibility
  
  for (i in seq_len(n_perm)) {
    perm_mat[, i] <- perm_one(df, codons)
  }
  
  
  ###############################################
  ### 3. Observed codon counts among "up" sites
  ###############################################
  
  # Observed counts for codons in UP sites
  obs_counts <- table(df$codon[df$direction == "up"])
  
  # Align with full codon list, fill missing with 0
  obs_full <- setNames(rep(0, length(codons)), codons)
  obs_full[names(obs_counts)] <- as.numeric(obs_counts)
  
  # Observed proportions among UP sites
  obs_prop <- obs_full / sum(obs_full)
  
  
  ###############################################
  ### 4. Permutation summary and p-values
  ###############################################
  
  # Expected mean count under permutation
  perm_mean <- rowMeans(perm_mat)
  
  # Standard deviation across permutations
  perm_sd <- apply(perm_mat, 1, sd)
  
  # Z-score enrichment (may be NA/Inf for rare codons)
  z_scores <- (obs_full - perm_mean) / perm_sd
  z_scores[is.infinite(z_scores)] <- NA
  z_scores[is.nan(z_scores)] <- NA
  
  ###############################################
  ### 4. Permutation summary and *two-sided* enrichment p-values
  ###############################################
  
  # Up-enrichment p-values (observed >= perm)
  p_perm_up <- sapply(codons, function(cd) {   
    mean(perm_mat[cd, ] >= obs_full[cd])    ### DO WE HAVE TO PUT AN = (<=/=>)
  })
  
  # Down-enrichment p-values (observed <= perm)
  p_perm_down <- sapply(codons, function(cd) {
    mean(perm_mat[cd, ] < obs_full[cd])
  })
  
  
  # NEW results table
  results <- data.frame(
    codon      = codons,
    obs_count  = obs_full,
    obs_prop   = obs_prop,
    exp_count  = perm_mean,
    exp_sd     = perm_sd,
    z_score    = z_scores,
    p_up       = p_perm_up,
    p_down     = p_perm_down,
    stringsAsFactors = FALSE
  )
  
  # Order results by most significant enrichment (minimum p-value)
  results$best_p <- pmin(results$p_up, results$p_down)
  results <- results[order(results$best_p), ]
  
  
  ###############################################
  ### 5. Codon frequencies in up vs non-significant (ns) sites
  ###############################################
  
  # IMPORTANT:
  # Here we define n.s. sites as those with direction == "ns".
  # If instead you want padj >= 0.05 as n.s., replace the ns_counts line accordingly.
  
  # Counts in UP sites (re-use obs_counts)
  up_counts <- obs_counts
  
  # Counts in NS sites:
  ns_counts <- table(df$codon[df$direction == "ns"])
  # Alternative if you prefer padj-based NS:
  # ns_counts <- table(df$codon[df$padj >= 0.05])
  
  # Full set of codons (use same codons vector)
  all_codons <- codons
  
  # Expand UP counts to all codons
  up_full_for_plot <- setNames(rep(0, length(all_codons)), all_codons)
  up_full_for_plot[names(up_counts)] <- as.numeric(up_counts)
  
  # Expand NS counts to all codons
  ns_full_for_plot <- setNames(rep(0, length(all_codons)), all_codons)
  ns_full_for_plot[names(ns_counts)] <- as.numeric(ns_counts)
  
  # Convert to proportions (0 if no sites)
  up_prop_plot <- up_full_for_plot / sum(up_full_for_plot)
  ns_prop_plot <- ns_full_for_plot / sum(ns_full_for_plot)
  
  
  ###############################################
  ### 6. Build plotting data frame and annotate significance
  ###############################################
  
  plot_df <- data.frame(
    codon   = all_codons,
    prop_ns = as.numeric(ns_prop_plot),
    prop_up = as.numeric(up_prop_plot),
    stringsAsFactors = FALSE
  )
  

  ###############################################
  ### Annotate significance category
  ###############################################
  
  plot_df <- merge(
    plot_df,
    results[, c("codon", "p_up", "p_down")],
    by = "codon",
    all.x = TRUE
    # suffixes = c("", "_res")
  )
  
  # Ensure numeric
  plot_df$p_up   <- as.numeric(plot_df$p_up)
  plot_df$p_down <- as.numeric(plot_df$p_down)
  
  # Set significance category
  alpha <- 0.05
  plot_df$category <- "not significant"
  plot_df$category[plot_df$p_up   <= alpha] <- "up_enriched"
  plot_df$category[plot_df$p_down <= alpha] <- "down_enriched"
  
  ###############################################
  ### 7. Scatter plot: %codon in ns vs %codon in up
  ###############################################
  
  library(ggplot2)
  n_site_up  <- sum(df$direction == "up")
  n_site_ns  <- sum(df$direction == "ns")
  
  n_gene_up  <- length(unique(df$transcript_id[df$direction == "up"]))
  n_gene_ns  <- length(unique(df$transcript_id[df$direction == "ns"]))
  
  
  p <- ggplot(plot_df, aes(x = prop_ns * 100,
                           y = prop_up * 100,
                           color = category)) +
    geom_point(size = 3, alpha = 0.8) +
    # Diagonal: where up % == ns %
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey40") +
    # Label only significant codons
    geom_text_repel(
      data = subset(plot_df, category != "not significant"),
      aes(label = codon),
      size = 3,
      max.overlaps = 100,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey60",
      show.legend = FALSE
    )+
    scale_color_manual(values = c(
      "up_enriched"     = "blue",
      "down_enriched"   = "red",
      "not significant" = "lightgrey"
    ))+
    labs(
      x = paste0(
        "Codon frequency in non-significant sites (%)\n",
        "n = ", n_site_ns, " sites | ", n_gene_ns, " genes"
      ),
      y = paste0(
        "Codon frequency in up sites (%)\n",
        "n = ", n_site_up, " sites | ", n_gene_up, " genes"
      ),
      color = "Permutation\nsignificance",
      title = "Codon usage: up sites vs non-significant sites"
    )+
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )+
    xlim(0,15)+
    ylim(0,15)
  
  
  print(p)
  
  ggsave(p, file=paste0("~/save_file.pdf"),width = 8, height = 8, units = "in")
  
  
  ###############################################
  ### 8. Calcul log2FC codon usage (UP vs NS)
  ###############################################
  
  ###############################################
  ### 8. Calcul log2FC codon usage (UP vs NS)
  ###############################################
  
plot_df$log2FC <- with(
    plot_df,
    ifelse(prop_ns == 0, NA, log2(prop_up / prop_ns))
  )
  
  plot_df$condition <- y
  
  log2fc_df <- plot_df[, c("condition", "codon", "log2FC", "category")]

    all_log2fc_list[[y]] <- log2fc_df
  
  
}

###############################################
### 9. Heatmap global : conditions × codons
###############################################

all_log2fc <- do.call(rbind, all_log2fc_list)

all_log2fc$condition_clean <- all_log2fc$condition %>%
  stringr::str_replace("_permut$", "")  # par ex. enlever "_permut" si tu veux

all_log2fc$codon <- factor(all_log2fc$codon, levels = sort(unique(all_log2fc$codon)))

all_log2fc$condition_clean <- factor(all_log2fc$condition_clean,
                                     levels = unique(all_log2fc$condition_clean))

###############################################
### order conditions : *_mono then *_di
###############################################

conds <- unique(all_log2fc$condition_clean)

conds_mono <- conds[grepl("mono", conds)]
conds_di   <- conds[grepl("di", conds)]

condition_levels <- c(conds_mono, conds_di)

all_log2fc$condition_clean <- factor(
  all_log2fc$condition_clean,
  levels = condition_levels
)

#order by aa
all_log2fc <- merge(all_log2fc, codon_map, by = "codon")

all_log2fc$codon_ordered <- with(
  all_log2fc, 
  factor(codon, levels = unique(codon[order(aa, codon)]))
)

p_heat_global <- ggplot(all_log2fc, aes(x = codon_ordered, y = condition_clean, fill = log2FC)) +
  geom_tile(color = "white") +
  
  # Significant dots
  geom_point(data = subset(all_log2fc, category != "not significant"), 
             aes(x = codon_ordered, y = condition_clean), 
             color = "black", size = 1, shape = 16) +
  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  
  # Grouping by Amino Acid
  facet_grid(. ~ aa, scales = "free_x", space = "free_x") +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    panel.grid = element_blank(),
    # Styling for the facet headers (AA names)
    strip.text.x = element_text(angle = 90, size = 8, face = "bold"),
    panel.spacing = unit(0.1, "lines") 
  ) +
  labs(x = "Codon (grouped by AA)", y = "Condition")

print(p_heat_global)


ggsave(
  p_heat_global,
  file = "~/path/save.pdf",
  width = 11, height = 6
)

