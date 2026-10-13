library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggforce)
library(GenomicRanges)
library(stringr)
library(annotatr)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(VennDiagram)

setwd("/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2025.10_RRBS_CpG_methylation/Results/Analysis")
methylation_matrix <- read.csv("rrbs_CpG_matrix.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("rrbs_meta.csv", header = TRUE, stringsAsFactors = FALSE)
stopifnot(all(colnames(methylation_matrix) %in% metadata$sample))
annots <- build_annotations(
  genome = "mm10", 
  annotations = c("mm10_basicgenes", "mm10_enhancers_fantom")
)
promoter_ranges <- annots[annots$type == "mm10_genes_promoters"]
gene_body_ranges <- annots[annots$type %in% c("mm10_genes_exons","mm10_genes_introns",
                                              "mm10_genes_5UTRs","mm10_genes_3UTRs")]
enhancer_ranges <- annots[annots$type == "mm10_enhancers_fantom"]
coords <- rownames(methylation_matrix)
df <- str_split_fixed(coords, ":", 2) %>% as.data.frame(stringsAsFactors = FALSE)
colnames(df) <- c("chr","pos")
df$pos <- as.integer(df$pos)
gr_cpgs <- GRanges(
  seqnames = paste0("chr", df$chr),
  ranges   = IRanges(start = df$pos, end = df$pos),
  strand   = "*"
)
hits_promoter    <- findOverlaps(gr_cpgs, promoter_ranges)
hits_gene_body   <- findOverlaps(gr_cpgs, gene_body_ranges)
hits_enhancer    <- findOverlaps(gr_cpgs, enhancer_ranges)

promoter_meth_list <- split(queryHits(hits_promoter), subjectHits(hits_promoter))
meth_promoter_avg <- sapply(promoter_meth_list, function(cpg_rows) {
  colMeans(methylation_matrix[cpg_rows, , drop = FALSE], na.rm = TRUE)
})
meth_promoter_avg <- t(meth_promoter_avg)
cat("Promoter matrix (averaged):", nrow(meth_promoter_avg), "promoters x", ncol(meth_promoter_avg), "samples\n")

gene_body_meth_list <- split(queryHits(hits_gene_body), subjectHits(hits_gene_body))
meth_gene_body_avg <- sapply(gene_body_meth_list, function(cpg_rows) {
  colMeans(methylation_matrix[cpg_rows, , drop = FALSE], na.rm = TRUE)
})
meth_gene_body_avg <- t(meth_gene_body_avg)
cat("Gene-body matrix (averaged per feature):", nrow(meth_gene_body_avg), "features x", ncol(meth_gene_body_avg), "samples\n")
gene_symbols <- gene_body_ranges$symbol[subjectHits(hits_gene_body)]
gene_cpg_map <- data.frame(
  cpg_idx = queryHits(hits_gene_body),
  gene = gene_symbols,
  stringsAsFactors = FALSE
)
gene_cpg_list <- split(gene_cpg_map$cpg_idx, gene_cpg_map$gene)
meth_gene_avg <- sapply(gene_cpg_list, function(cpg_rows) {
  colMeans(methylation_matrix[cpg_rows, , drop = FALSE], na.rm = TRUE)
})
meth_gene_avg <- t(meth_gene_avg)
cat("Gene-body matrix (averaged per gene):", nrow(meth_gene_avg), "genes x", ncol(meth_gene_avg), "samples\n")

enhancer_meth_list <- split(queryHits(hits_enhancer), subjectHits(hits_enhancer))
meth_enhancer_avg <- sapply(enhancer_meth_list, function(cpg_rows) {
  colMeans(methylation_matrix[cpg_rows, , drop = FALSE], na.rm = TRUE)
})
meth_enhancer_avg <- t(meth_enhancer_avg)
cat("Enhancer matrix (averaged):", nrow(meth_enhancer_avg), "enhancers x", ncol(meth_enhancer_avg), "samples\n")

avg_global <- colMeans(methylation_matrix, na.rm = TRUE)
avg_promoter <- colMeans(meth_promoter_avg, na.rm = TRUE)
avg_gene_body <- colMeans(meth_gene_avg, na.rm = TRUE)
avg_enhancer <- colMeans(meth_enhancer_avg, na.rm = TRUE)
avg_meth_df <- data.frame(
  sample = colnames(methylation_matrix),
  global    = avg_global,
  promoter  = avg_promoter,
  gene_body = avg_gene_body,
  enhancer  = avg_enhancer
)
avg_meth_df <- left_join(avg_meth_df, metadata, by = "sample")
avg_meth_df <- avg_meth_df %>%
  mutate(Group = sub("(_\\d+)$", "", sample))
avg_meth_long <- avg_meth_df %>%
  pivot_longer(
    cols = c(global, promoter, gene_body, enhancer),
    names_to = "feature_type",
    values_to = "avg_methylation"
  )
avg_group <- avg_meth_long %>%
  group_by(Group, feature_type) %>%
  summarise(
    mean_methylation = mean(avg_methylation, na.rm = TRUE),
    sd_methylation   = sd(avg_methylation, na.rm = TRUE),
    n                = n(),
    .groups          = "drop"
  )
plot_avg_methylation_box <- function(avg_meth_df, conditions_to_plot = NULL, title = "Average methylation per group") {
    if(!"Group" %in% colnames(avg_meth_df)) {
    avg_meth_df <- avg_meth_df %>%
      mutate(Group = sub("(_\\d+)$", "", sample))
  }
    avg_meth_long <- avg_meth_df %>%
    pivot_longer(cols = c(global, promoter, gene_body, enhancer),
                 names_to = "feature_type",
                 values_to = "avg_methylation")
    if(!is.null(conditions_to_plot)) {
    # trim spaces just in case
    avg_meth_long <- avg_meth_long %>%
      mutate(Group = trimws(Group)) %>%
      dplyr::filter(Group %in% conditions_to_plot)
  }
    custom_colors <- c(
    "YY1_ctrl" = "dodgerblue4",
    "YY1_tev" = "dodgerblue1",
    "YY12_ctrl" = "darkolivegreen4",
    "YY12_tev" = "darkolivegreen2",
    "OO12_ctrl" = "tomato4",
    "OO12_tev" = "tomato1"
  )
  custom_colors <- custom_colors[names(custom_colors) %in% unique(avg_meth_long$Group)]
    ggplot(avg_meth_long, aes(x = Group, y = avg_methylation, fill = Group)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +        # boxplot without showing outliers twice
    geom_jitter(aes(color = Group), width = 0.15, size = 3) + # show replicates
    facet_wrap(~feature_type, scales = "free_y") +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +
    labs(
      title = title,
      x = "Group",
      y = "Average methylation",
      fill = "Group",
      color = "Group"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
conditions_to_plot <- c("OO12_ctrl", "OO12_tev", "YY1_ctrl", "YY1_tev", "YY12_ctrl", "YY12_tev")
plot_avg_methylation_box(avg_meth_df, conditions_to_plot,
                         title = "Average methylation per group (selected conditions)")
export_df <- avg_meth_df %>% 
  dplyr::filter(Group %in% conditions_to_plot)
write.csv(export_df, "methylation_per_replicate.csv", row.names = FALSE)

pval_thresh <- 0.05
fc_thresh   <- 0.0
metadata <- metadata %>%
  mutate(
    Group2 = sub("(_\\d+)$", "", sample),
    AnalysisGroup = ifelse(Group2 %in% c("YY1_ctrl", "YY12_ctrl"), "Young_ctrl", Group2)
  )
groups_to_compare <- c("YY1_tev", "Young_ctrl")
meta_promoter <- metadata %>% dplyr::filter(AnalysisGroup %in% groups_to_compare)
epsilon <- 1e-6
beta_global_adj <- pmin(pmax(methylation_matrix, epsilon), 1 - epsilon)
M_global_all <- log2(beta_global_adj / (1 - beta_global_adj))
keep <- complete.cases(M_global_all)
M_global_all <- M_global_all[keep, ]
beta_global_all <- methylation_matrix[keep, ]
cat("Total global CpGs after removing NAs:", nrow(M_global_all), "\n")
beta_global <- M_global_all[, meta_promoter$sample, drop = FALSE]
M_global <- beta_global
group_factor_global <- factor(meta_promoter$AnalysisGroup, levels = groups_to_compare)
design <- model.matrix(~0 + group_factor_global)
colnames(design) <- make.names(levels(group_factor_global))
contr <- makeContrasts(
  Young_tev_vs_Young = YY1_tev - Young_ctrl,
  levels = design
)
fit <- lmFit(M_global, design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
res_global <- topTable(fit2, coef = "Young_tev_vs_Young", number = Inf)
beta_young <- rowMeans(beta_global_all[, metadata$AnalysisGroup == "Young_ctrl", drop = FALSE])
beta_young_tev   <- rowMeans(beta_global_all[, metadata$AnalysisGroup == "YY1_tev", drop = FALSE])
res_global$delta_beta <- beta_young_tev[rownames(res_global)] - beta_young[rownames(res_global)]
res_global$cpg_id <- rownames(res_global)
cat("Total global CpGs tested:", nrow(M_global), "\n")
sig_global <- res_global %>%
  dplyr::filter(P.Value < pval_thresh & abs(delta_beta) >= fc_thresh)
M_global_sig <- M_global_all[rownames(M_global_all) %in% sig_global$cpg_id, ]
cat("Number of significant global CpGs (P < 0.05):", nrow(M_global_sig), "\n")
samples_to_use <- c(
  "YY1_tev_1","YY1_tev_2","YY1_tev_3","YY1_tev_4",
  "YY1_ctrl_1","YY1_ctrl_2","YY1_ctrl_3","YY1_ctrl_4",
  "YY12_ctrl_1","YY12_ctrl_2","YY12_ctrl_3","YY12_ctrl_4",
  "YY12_tev_1","YY12_tev_2","YY12_tev_3","YY12_tev_4"
)
samples_to_use <- samples_to_use[samples_to_use %in% colnames(M_global_sig)]
M_global_heat <- M_global_sig[, samples_to_use, drop = FALSE]
annotation_col <- data.frame(
  AnalysisGroup = metadata$AnalysisGroup,
  row.names = metadata$sample
)
annotation_col <- annotation_col[colnames(M_global_heat), , drop = FALSE] 
ann_colors <- list(AnalysisGroup = c(
  "YY1_ctrl"   = "dodgerblue1",
  "YY1_tev"    = "dodgerblue4",
  "YY12_ctrl"  = "darkolivegreen2",
  "YY12_tev"   = "darkolivegreen4",
  "OO12_ctrl"  = "tomato4",
  "OO12_tev"   = "tomato1",
  "Young_ctrl" = "lightblue"          
))
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(
  M_global_heat,            
  kmeans_k = 10,           
  scale = "row",
  cluster_rows = TRUE,     
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,     
  color = heat_colors,
  main = "Top global CpG methylation (K-means 10 clusters; P < 0.05)"
)
write.csv(M_global_heat, "titin_cleavage_associated_average_methylation.csv", row.names = TRUE)

M_global_pca <- M_global_heat
M_t <- t(M_global_pca)
pca_res <- prcomp(M_t, center = TRUE, scale. = TRUE)
pca_coords <- as.data.frame(pca_res$x)
pca_coords$sample <- rownames(pca_coords)
pca_coords <- left_join(
  pca_coords, 
  metadata %>% dplyr::select(sample, AnalysisGroup), 
  by = "sample"
)
pca_var <- pca_res$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)
cat("PC1:", pca_var_perc[1], "%, PC2:", pca_var_perc[2], "%\n")
ggplot(pca_coords, aes(x = PC1, y = PC2, color = AnalysisGroup)) + 
  geom_point(size = 4) +
  geom_mark_ellipse(aes(fill = AnalysisGroup),  
                    alpha = 0.2,          
                    show.legend = FALSE,  
                    expand = unit(0.5, "mm")) + 
  coord_cartesian(xlim = range(pca_coords$PC1) + c(-1,1),
                  ylim = range(pca_coords$PC2) + c(-5,5)) +
  theme_minimal() +
  labs(
    title = "PCA of differential global CpGs",
    x = paste0("PC1 (", pca_var_perc[1], "%)"),
    y = paste0("PC2 (", pca_var_perc[2], "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c(
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue"  
  )) +
  scale_fill_manual(values = c(            
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue"  
  ))


avg_sig_global <- colMeans(methylation_matrix[sig_global$cpg_id, , drop = FALSE], na.rm = TRUE)
avg_meth_sig_df <- data.frame(
  sample    = colnames(methylation_matrix),
  global    = avg_sig_global
)
avg_meth_sig_df <- left_join(avg_meth_sig_df, metadata, by = "sample")
if(!"Group" %in% colnames(avg_meth_sig_df)) {
  avg_meth_sig_df <- avg_meth_sig_df %>%
    mutate(Group = sub("(_\\d+)$", "", sample))
}
conditions_to_plot <- c("YY1_ctrl", "YY1_tev", "YY12_ctrl", "YY12_tev")
export_sig_df <- avg_meth_sig_df %>% 
  dplyr::filter(Group2 %in% conditions_to_plot)
write.csv(export_sig_df, "mKO_associated_methylation_per_replicate.csv", row.names = FALSE)

detect_restoration <- function(beta_matrix, metadata, feature_ids) {
  
  young_samples <- metadata$sample[metadata$AnalysisGroup == "Young_ctrl"]
  old_samples   <- metadata$sample[metadata$AnalysisGroup == "YY1_tev"]
  rej_samples   <- metadata$sample[metadata$AnalysisGroup == "YY12_tev"]
  
  beta_young <- rowMeans(beta_matrix[, young_samples, drop = FALSE])
  beta_old   <- rowMeans(beta_matrix[, old_samples, drop = FALSE])
  beta_rej   <- rowMeans(beta_matrix[, rej_samples, drop = FALSE])
  
  d_old     <- beta_old - beta_young         
  d_rej_old <- beta_rej - beta_old           
  d_rej_yng <- beta_rej - beta_young         
  
  restored <- (sign(d_old) != sign(d_rej_old)) &  
    (abs(d_rej_yng) < abs(d_old))                 
  out <- data.frame(
    feature_id = feature_ids,
    beta_young = beta_young,
    beta_old   = beta_old,
    beta_rej   = beta_rej,
    delta_old  = d_old,
    delta_rej  = d_rej_old,
    restored   = restored,
    restoration_type = dplyr::case_when(
      restored & d_old > 0 ~ "hyper_in_mKO_reversed",
      restored & d_old < 0 ~ "hypo_in_mKO_reversed",
      TRUE ~ "not_restored"
    )
  )
  return(out)
}

detect_restoration_categories <- function(df, aging_thresh = 0.05, rej_thresh = 0.05) {
  required_cols <- c("delta_old", "delta_rej")
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns: delta_old or delta_rej")
  }
  df <- df %>%
    mutate(
      aging_direction = case_when(
        delta_old >  aging_thresh ~ "hyper_in_mKO",
        delta_old < -aging_thresh ~ "hypo_in_mKO",
        TRUE ~ "no_change"
      ),
      restored_simple = case_when(
        aging_direction == "hyper_in_mKO" & delta_rej < -rej_thresh ~ "restored",
        aging_direction == "hypo_in_mKO"  & delta_rej >  rej_thresh ~ "restored",
        aging_direction %in% c("hyper_in_mKO","hypo_in_mKO") ~ "not_restored",
        TRUE ~ "no_change"
      ),
      restoration_category = paste0(aging_direction, "_", restored_simple)
    )
  return(df)
}
rej_thresh <- 0.05
global_sig_mat <- methylation_matrix[sig_global$cpg_id, , drop = FALSE]
rest_global_sig <- detect_restoration(global_sig_mat, metadata, rownames(global_sig_mat)) %>%
  detect_restoration_categories(aging_thresh = 0.00, rej_thresh = rej_thresh)
all_rest <- bind_rows(
  rest_global_sig    %>% mutate(feature_type = "GlobalCpG")
)
summary_table <- all_rest %>%
  group_by(feature_type, restoration_category) %>%
  summarise(n_features = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = restoration_category, values_from = n_features, values_fill = 0)
print(summary_table)

pval_thresh <- 0.05
fc_thresh   <- 0.0
metadata <- metadata %>%
  mutate(
    Group2 = sub("(_\\d+)$", "", sample),
    AnalysisGroup = ifelse(Group2 %in% c("YY1_ctrl", "YY12_ctrl"), "Young_ctrl", Group2)
  )

groups_to_compare <- c("OO12_ctrl", "Young_ctrl")
meta_promoter <- metadata %>% dplyr::filter(AnalysisGroup %in% groups_to_compare)
epsilon <- 1e-6
beta_global_adj <- pmin(pmax(methylation_matrix, epsilon), 1 - epsilon)
M_global_all <- log2(beta_global_adj / (1 - beta_global_adj))
keep <- complete.cases(M_global_all)
M_global_all <- M_global_all[keep, ]
beta_global_all <- methylation_matrix[keep, ]
cat("Total global CpGs after removing NAs:", nrow(M_global_all), "\n")
beta_global <- M_global_all[, meta_promoter$sample, drop = FALSE]
M_global <- beta_global
group_factor_global <- factor(meta_promoter$AnalysisGroup, levels = groups_to_compare)
design <- model.matrix(~0 + group_factor_global)
colnames(design) <- make.names(levels(group_factor_global))
contr <- makeContrasts(
  Old_vs_Young = OO12_ctrl - Young_ctrl,
  levels = design
)
fit <- lmFit(M_global, design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
res_global <- topTable(fit2, coef = "Old_vs_Young", number = Inf)
beta_young <- rowMeans(beta_global_all[, metadata$AnalysisGroup == "Young_ctrl", drop = FALSE])
beta_old   <- rowMeans(beta_global_all[, metadata$AnalysisGroup == "OO12_ctrl", drop = FALSE])
res_global$delta_beta <- beta_old[rownames(res_global)] - beta_young[rownames(res_global)]
res_global$cpg_id <- rownames(res_global)
cat("Total global CpGs tested:", nrow(M_global), "\n")
sig_global <- res_global %>%
  dplyr::filter(P.Value < pval_thresh & abs(delta_beta) >= fc_thresh)
M_global_sig <- M_global_all[rownames(M_global_all) %in% sig_global$cpg_id, ]
cat("Number of significant global CpGs (P < 0.05):", nrow(M_global_sig), "\n")
samples_to_use <- c(
  "OO12_ctrl_1","OO12_ctrl_2","OO12_ctrl_3","OO12_ctrl_4",
  "YY1_ctrl_1","YY1_ctrl_2","YY1_ctrl_3","YY1_ctrl_4",
  "YY12_ctrl_1","YY12_ctrl_2","YY12_ctrl_3","YY12_ctrl_4",
  "OO12_tev_1","OO12_tev_2","OO12_tev_3","OO12_tev_4"
)
samples_to_use <- samples_to_use[samples_to_use %in% colnames(M_global_sig)]
M_global_heat <- M_global_sig[, samples_to_use, drop = FALSE]
annotation_col <- data.frame(
  AnalysisGroup = metadata$AnalysisGroup,
  row.names = metadata$sample
)
annotation_col <- annotation_col[colnames(M_global_heat), , drop = FALSE] 
ann_colors <- list(AnalysisGroup = c(
  "YY1_ctrl"   = "dodgerblue1",
  "YY1_tev"    = "dodgerblue4",
  "YY12_ctrl"  = "darkolivegreen2",
  "YY12_tev"   = "darkolivegreen4",
  "OO12_ctrl"  = "tomato4",
  "OO12_tev"   = "tomato1",
  "Young_ctrl" = "lightblue"          
))
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(
  M_global_heat,            
  kmeans_k = 10,           
  scale = "row",
  cluster_rows = TRUE,     
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,     
  color = heat_colors,
  main = "Top global CpG methylation (K-means 10 clusters; P < 0.05)"
)
write.csv(M_global_heat, "ageing_associated_average_methylation.csv", row.names = TRUE)

M_global_pca <- M_global_heat
M_t <- t(M_global_pca)
pca_res <- prcomp(M_t, center = TRUE, scale. = TRUE)
pca_coords <- as.data.frame(pca_res$x)
pca_coords$sample <- rownames(pca_coords)
pca_coords <- left_join(
  pca_coords, 
  metadata %>% dplyr::select(sample, AnalysisGroup), 
  by = "sample"
)
pca_var <- pca_res$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)
cat("PC1:", pca_var_perc[1], "%, PC2:", pca_var_perc[2], "%\n")
ggplot(pca_coords, aes(x = PC1, y = PC2, color = AnalysisGroup)) + 
  geom_point(size = 4) +
  geom_mark_ellipse(aes(fill = AnalysisGroup),  
                    alpha = 0.2,          
                    show.legend = FALSE,  
                    expand = unit(0.5, "mm")) + 
  coord_cartesian(xlim = range(pca_coords$PC1) + c(-1,1),
                  ylim = range(pca_coords$PC2) + c(-5,5)) +
  theme_minimal() +
  labs(
    title = "PCA of differential global CpGs",
    x = paste0("PC1 (", pca_var_perc[1], "%)"),
    y = paste0("PC2 (", pca_var_perc[2], "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c(
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue"  
  )) +
  scale_fill_manual(values = c(            
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue"  
  ))

metadata <- metadata %>%
  mutate(
    Group2 = sub("(_\\d+)$", "", sample),
    # Combine young controls into 'Young_ctrl', keep 'OO12_ctrl' separate
    AnalysisGroup = ifelse(Group2 %in% c("YY1_ctrl", "YY12_ctrl"), "Young_ctrl", Group2)
  )
groups_to_compare <- c("OO12_ctrl", "Young_ctrl")
meta_promoter <- metadata %>% dplyr::filter(AnalysisGroup %in% groups_to_compare)
epsilon <- 1e-6
beta_promoter_adj <- pmin(pmax(meth_promoter_avg, epsilon), 1 - epsilon)
M_promoter_all <- log2(beta_promoter_adj / (1 - beta_promoter_adj))
keep <- complete.cases(M_promoter_all)
M_promoter_all <- M_promoter_all[keep, ]
beta_promoter_all <- meth_promoter_avg[keep, ]
cat("Total promoters after removing NAs:", nrow(M_promoter_all), "\n")
beta_promoter <- M_promoter_all[, meta_promoter$sample, drop = FALSE]
M_promoter <- beta_promoter
group_factor_promoter <- factor(meta_promoter$AnalysisGroup, levels = groups_to_compare)
design <- model.matrix(~0 + group_factor_promoter)
colnames(design) <- make.names(levels(group_factor_promoter))
contr <- makeContrasts(
  Old_vs_Young = OO12_ctrl - Young_ctrl,
  levels = design
)
fit <- lmFit(M_promoter, design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
res_promoter <- topTable(fit2, coef = "Old_vs_Young", number = Inf)
beta_young <- rowMeans(beta_promoter_all[, metadata$AnalysisGroup == "Young_ctrl", drop = FALSE])
beta_old   <- rowMeans(beta_promoter_all[, metadata$AnalysisGroup == "OO12_ctrl", drop = FALSE])
res_promoter$delta_beta <- beta_old[rownames(res_promoter)] - beta_young[rownames(res_promoter)]
res_promoter$promoter_id <- rownames(res_promoter)
res_promoter$gene_symbol <- promoter_ranges$symbol[
  match(res_promoter$promoter_id, names(promoter_ranges))
]
cat("Total promoters tested:", nrow(M_promoter), "\n")
sig_promoters <- res_promoter %>%
  dplyr::filter(P.Value < pval_thresh & abs(delta_beta) >= fc_thresh)
M_promoter_sig <- M_promoter_all[rownames(M_promoter_all) %in% sig_promoters$promoter_id, ]
cat("Number of significant promoters (P < 0.05):", nrow(M_promoter_sig), "\n")
samples_to_use <- c(
  "OO12_ctrl_1","OO12_ctrl_2","OO12_ctrl_3","OO12_ctrl_4",
  "YY1_ctrl_1","YY1_ctrl_2","YY1_ctrl_3","YY1_ctrl_4",
  "YY12_ctrl_1","YY12_ctrl_2","YY12_ctrl_3","YY12_ctrl_4",
  "OO12_tev_1","OO12_tev_2","OO12_tev_3","OO12_tev_4"
)
samples_to_use <- samples_to_use[samples_to_use %in% colnames(M_promoter_sig)]
M_promoter_heat <- M_promoter_sig[, samples_to_use, drop = FALSE]
annotation_col <- data.frame(
  AnalysisGroup = metadata$AnalysisGroup,
  row.names = metadata$sample
)
annotation_col <- annotation_col[colnames(M_promoter_heat), , drop = FALSE] 
ann_colors <- list(AnalysisGroup = c(
  "YY1_ctrl"   = "dodgerblue1",
  "YY1_tev"    = "dodgerblue4",
  "YY12_ctrl"  = "darkolivegreen2",
  "YY12_tev"   = "darkolivegreen4",
  "OO12_ctrl"  = "tomato4",
  "OO12_tev"   = "tomato1",
  "Young_ctrl" = "lightblue"          
))
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(
  M_promoter_heat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,     
  color = heat_colors,
  main = "Top promoter methylation (P < 0.05)"
)

M_promoter_pca <- M_promoter_heat
M_t <- t(M_promoter_pca)
pca_res <- prcomp(M_t, center = TRUE, scale. = TRUE)
pca_coords <- as.data.frame(pca_res$x)
pca_coords$sample <- rownames(pca_coords)
pca_coords <- left_join(
  pca_coords, 
  metadata %>% dplyr::select(sample, AnalysisGroup), # Changed from Group2
  by = "sample"
)
pca_var <- pca_res$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)
cat("PC1:", pca_var_perc[1], "%, PC2:", pca_var_perc[2], "%\n")
ggplot(pca_coords, aes(x = PC1, y = PC2, color = AnalysisGroup)) + 
  geom_point(size = 4) +
  geom_mark_ellipse(aes(fill = AnalysisGroup),
                    alpha = 0.2,          
                    show.legend = FALSE,  
                    expand = unit(0.5, "mm")) + 
  coord_cartesian(xlim = range(pca_coords$PC1) + c(-1,1),
                  ylim = range(pca_coords$PC2) + c(-5,5)) +
  theme_minimal() +
  labs(
    title = "PCA of differential promoters",
    x = paste0("PC1 (", pca_var_perc[1], "%)"),
    y = paste0("PC2 (", pca_var_perc[2], "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c(
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue" 
  )) +
  scale_fill_manual(values = c(            
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue"  
  ))

metadata <- metadata %>%
  mutate(
    Group2 = sub("(_\\d+)$", "", sample),
    AnalysisGroup = ifelse(Group2 %in% c("YY1_ctrl", "YY12_ctrl"), "Young_ctrl", Group2)
  )
groups_to_compare <- c("OO12_ctrl", "Young_ctrl")
meta_gene <- metadata %>% dplyr::filter(AnalysisGroup %in% groups_to_compare)
epsilon <- 1e-6
beta_gene_adj <- pmin(pmax(meth_gene_avg, epsilon), 1 - epsilon)   # gene-level beta matrix
M_gene_all <- log2(beta_gene_adj / (1 - beta_gene_adj))
keep <- complete.cases(M_gene_all)
M_gene_all <- M_gene_all[keep, ]
beta_gene_all <- meth_gene_avg[keep, ]
cat("Total genes after removing NAs:", nrow(M_gene_all), "\n")
M_gene <- M_gene_all[, meta_gene$sample, drop = FALSE]
group_factor_gene <- factor(meta_gene$AnalysisGroup, levels = groups_to_compare)
design <- model.matrix(~0 + group_factor_gene)
colnames(design) <- make.names(levels(group_factor_gene))
contr <- makeContrasts(
  Old_vs_Young = OO12_ctrl - Young_ctrl,
  levels = design
)
fit <- lmFit(M_gene, design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
res_gene <- topTable(fit2, coef = "Old_vs_Young", number = Inf)
beta_young <- rowMeans(beta_gene_all[, metadata$AnalysisGroup == "Young_ctrl", drop = FALSE])
beta_old   <- rowMeans(beta_gene_all[, metadata$AnalysisGroup == "OO12_ctrl", drop = FALSE])
res_gene$delta_beta <- beta_old[rownames(res_gene)] - beta_young[rownames(res_gene)]
res_gene$gene_id <- rownames(res_gene)
res_gene$gene_symbol <- rownames(res_gene)
cat("Total genes tested:", nrow(M_gene), "\n")
sig_genes <- res_gene %>%
  dplyr::filter(P.Value < pval_thresh & abs(delta_beta) >= fc_thresh)

M_gene_sig <- M_gene_all[rownames(M_gene_all) %in% sig_genes$gene_id, ]
cat("Number of significant genes (P < 0.05):", nrow(M_gene_sig), "\n")
samples_to_use <- c(
  "OO12_ctrl_1","OO12_ctrl_2","OO12_ctrl_3","OO12_ctrl_4",
  "YY1_ctrl_1","YY1_ctrl_2","YY1_ctrl_3","YY1_ctrl_4",
  "YY12_ctrl_1","YY12_ctrl_2","YY12_ctrl_3","YY12_ctrl_4",
  "OO12_tev_1","OO12_tev_2","OO12_tev_3","OO12_tev_4"
)
samples_to_use <- samples_to_use[samples_to_use %in% colnames(M_gene_sig)]
M_gene_heat <- M_gene_sig[, samples_to_use, drop = FALSE]
annotation_col <- data.frame(
  AnalysisGroup = metadata$AnalysisGroup,
  row.names = metadata$sample
)
annotation_col <- annotation_col[colnames(M_gene_heat), , drop = FALSE] 
ann_colors <- list(AnalysisGroup = c(
  "YY1_ctrl"   = "dodgerblue1",
  "YY1_tev"    = "dodgerblue4",
  "YY12_ctrl"  = "darkolivegreen2",
  "YY12_tev"   = "darkolivegreen4",
  "OO12_ctrl"  = "tomato4",
  "OO12_tev"   = "tomato1",
  "Young_ctrl" = "lightblue"
))
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(
  M_gene_heat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heat_colors,
  main = "Top gene methylation (P < 0.05)"
)

M_gene_pca <- M_gene_heat
M_t <- t(M_gene_pca)
pca_res <- prcomp(M_t, center = TRUE, scale. = TRUE)
pca_coords <- as.data.frame(pca_res$x)
pca_coords$sample <- rownames(pca_coords)
pca_coords <- left_join(
  pca_coords,
  metadata %>% dplyr::select(sample, AnalysisGroup),
  by = "sample"
)
pca_var <- pca_res$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)
cat("PC1:", pca_var_perc[1], "%, PC2:", pca_var_perc[2], "%\n")
ggplot(pca_coords, aes(x = PC1, y = PC2, color = AnalysisGroup)) + 
  geom_point(size = 4) +
  geom_mark_ellipse(aes(fill = AnalysisGroup),   
                    alpha = 0.2,          
                    show.legend = FALSE,  
                    expand = unit(0.5, "mm")) + 
  coord_cartesian(xlim = range(pca_coords$PC1) + c(-1,1),
                  ylim = range(pca_coords$PC2) + c(-5,5)) +
  theme_minimal() +
  labs(
    title = "PCA of differential gene body",
    x = paste0("PC1 (", pca_var_perc[1], "%)"),
    y = paste0("PC2 (", pca_var_perc[2], "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c(
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue"
  )) +
  scale_fill_manual(values = c(            
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue"
  ))

metadata <- metadata %>%
  mutate(
    Group2 = sub("(_\\d+)$", "", sample),
    AnalysisGroup = ifelse(Group2 %in% c("YY1_ctrl", "YY12_ctrl"), "Young_ctrl", Group2)
  )
groups_to_compare <- c("OO12_ctrl", "Young_ctrl")
meta_enhancer <- metadata %>% dplyr::filter(AnalysisGroup %in% groups_to_compare)
epsilon <- 1e-6
beta_enhancer_adj <- pmin(pmax(meth_enhancer_avg, epsilon), 1 - epsilon)
M_enhancer_all <- log2(beta_enhancer_adj / (1 - beta_enhancer_adj))
keep <- complete.cases(M_enhancer_all)
M_enhancer_all <- M_enhancer_all[keep, ]
beta_enhancer_all <- meth_enhancer_avg[keep, ]
cat("Total enhancers after removing NAs:", nrow(M_enhancer_all), "\n")
beta_enhancer <- M_enhancer_all[, meta_enhancer$sample, drop = FALSE]
M_enhancer <- beta_enhancer
group_factor_enhancer <- factor(meta_enhancer$AnalysisGroup, levels = groups_to_compare)
design <- model.matrix(~0 + group_factor_enhancer)
colnames(design) <- make.names(levels(group_factor_enhancer))
contr <- makeContrasts(
  Old_vs_Young = OO12_ctrl - Young_ctrl,
  levels = design
)
fit <- lmFit(M_enhancer, design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
res_enhancer <- topTable(fit2, coef = "Old_vs_Young", number = Inf)
beta_young <- rowMeans(beta_enhancer_all[, metadata$AnalysisGroup == "Young_ctrl", drop = FALSE])
beta_old   <- rowMeans(beta_enhancer_all[, metadata$AnalysisGroup == "OO12_ctrl", drop = FALSE])
res_enhancer$delta_beta <- beta_old[rownames(res_enhancer)] - beta_young[rownames(res_enhancer)]
res_enhancer$enhancer_id <- rownames(res_enhancer)
res_enhancer$gene_symbol <- enhancer_ranges$symbol[
  match(res_enhancer$enhancer_id, names(enhancer_ranges))
]
cat("Total enhancers tested:", nrow(M_enhancer), "\n")
sig_enhancers <- res_enhancer %>%
  dplyr::filter(P.Value < pval_thresh & abs(delta_beta) >= fc_thresh)
M_enhancer_sig <- M_enhancer_all[rownames(M_enhancer_all) %in% sig_enhancers$enhancer_id, ]
cat("Number of significant enhancers (P < 0.05):", nrow(M_enhancer_sig), "\n")
samples_to_use <- c(
  "OO12_ctrl_1","OO12_ctrl_2","OO12_ctrl_3","OO12_ctrl_4",
  "YY1_ctrl_1","YY1_ctrl_2","YY1_ctrl_3","YY1_ctrl_4",
  "YY12_ctrl_1","YY12_ctrl_2","YY12_ctrl_3","YY12_ctrl_4",
  "OO12_tev_1","OO12_tev_2","OO12_tev_3","OO12_tev_4"
)
samples_to_use <- samples_to_use[samples_to_use %in% colnames(M_enhancer_sig)]
M_enhancer_heat <- M_enhancer_sig[, samples_to_use, drop = FALSE]
annotation_col <- data.frame(
  AnalysisGroup = metadata$AnalysisGroup,
  row.names = metadata$sample
)
annotation_col <- annotation_col[colnames(M_enhancer_heat), , drop = FALSE] 
ann_colors <- list(AnalysisGroup = c(
  "YY1_ctrl"   = "dodgerblue1",
  "YY1_tev"    = "dodgerblue4",
  "YY12_ctrl"  = "darkolivegreen2",
  "YY12_tev"   = "darkolivegreen4",
  "OO12_ctrl"  = "tomato4",
  "OO12_tev"   = "tomato1",
  "Young_ctrl" = "lightblue"
))
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(
  M_enhancer_heat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heat_colors,
  main = "Top enhancer methylation (P < 0.05)"
)

M_enhancer_pca <- M_enhancer_heat
M_t <- t(M_enhancer_pca)
pca_res <- prcomp(M_t, center = TRUE, scale. = TRUE)
pca_coords <- as.data.frame(pca_res$x)
pca_coords$sample <- rownames(pca_coords)
pca_coords <- left_join(
  pca_coords,
  metadata %>% dplyr::select(sample, AnalysisGroup),
  by = "sample"
)
pca_var <- pca_res$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)
cat("PC1:", pca_var_perc[1], "%, PC2:", pca_var_perc[2], "%\n")
ggplot(pca_coords, aes(x = PC1, y = PC2, color = AnalysisGroup)) + 
  geom_point(size = 4) +
  geom_mark_ellipse(aes(fill = AnalysisGroup),   
                    alpha = 0.2,          
                    show.legend = FALSE,  
                    expand = unit(0.5, "mm")) + 
  coord_cartesian(xlim = range(pca_coords$PC1) + c(-5,3),
                  ylim = range(pca_coords$PC2) + c(-3,3)) +
  theme_minimal() +
  labs(
    title = "PCA of differential enhancers",
    x = paste0("PC1 (", pca_var_perc[1], "%)"),
    y = paste0("PC2 (", pca_var_perc[2], "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c(
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue"
  )) +
  scale_fill_manual(values = c(            
    "YY1_ctrl"   = "dodgerblue1",
    "YY1_tev"    = "dodgerblue4",
    "YY12_ctrl"  = "darkolivegreen2",
    "YY12_tev"   = "darkolivegreen4",
    "OO12_ctrl"  = "tomato4",
    "OO12_tev"   = "tomato1",
    "Young_ctrl" = "lightblue"
  ))


avg_sig_global <- colMeans(methylation_matrix[sig_global$cpg_id, , drop = FALSE], na.rm = TRUE)
avg_sig_promoter <- colMeans(meth_promoter_avg[sig_promoters$promoter_id, , drop = FALSE], na.rm = TRUE)
avg_sig_gene_body <- colMeans(meth_gene_avg[sig_genes$gene_symbol, , drop = FALSE], na.rm = TRUE)
avg_sig_enhancer <- colMeans(meth_enhancer_avg[sig_enhancers$enhancer_id, , drop = FALSE], na.rm = TRUE)
avg_meth_sig_df <- data.frame(
  sample    = colnames(methylation_matrix),
  global    = avg_sig_global,
  promoter  = avg_sig_promoter,
  gene_body = avg_sig_gene_body,
  enhancer  = avg_sig_enhancer
)
avg_meth_sig_df <- left_join(avg_meth_sig_df, metadata, by = "sample")
if(!"Group" %in% colnames(avg_meth_sig_df)) {
  avg_meth_sig_df <- avg_meth_sig_df %>%
    mutate(Group = sub("(_\\d+)$", "", sample))
}
conditions_to_plot <- c("OO12_ctrl", "OO12_tev", "YY1_ctrl", "YY1_tev", "YY12_ctrl", "YY12_tev")
export_sig_df <- avg_meth_sig_df %>% 
  dplyr::filter(Group2 %in% conditions_to_plot)
write.csv(export_sig_df, "aging_associated_methylation_per_replicate.csv", row.names = FALSE)

detect_restoration <- function(beta_matrix, metadata, feature_ids) {
  
  young_samples <- metadata$sample[metadata$AnalysisGroup == "Young_ctrl"]
  old_samples   <- metadata$sample[metadata$AnalysisGroup == "OO12_ctrl"]
  rej_samples   <- metadata$sample[metadata$AnalysisGroup == "OO12_tev"]
  
  beta_young <- rowMeans(beta_matrix[, young_samples, drop = FALSE])
  beta_old   <- rowMeans(beta_matrix[, old_samples, drop = FALSE])
  beta_rej   <- rowMeans(beta_matrix[, rej_samples, drop = FALSE])
  
  d_old     <- beta_old - beta_young         
  d_rej_old <- beta_rej - beta_old           
  d_rej_yng <- beta_rej - beta_young         
  
  restored <- (sign(d_old) != sign(d_rej_old)) & 
    (abs(d_rej_yng) < abs(d_old))                  
  
  out <- data.frame(
    feature_id = feature_ids,
    beta_young = beta_young,
    beta_old   = beta_old,
    beta_rej   = beta_rej,
    delta_old  = d_old,
    delta_rej  = d_rej_old,
    restored   = restored,
    restoration_type = dplyr::case_when(
      restored & d_old > 0 ~ "hyper_in_aging_reversed",
      restored & d_old < 0 ~ "hypo_in_aging_reversed",
      TRUE ~ "not_restored"
    )
  )
  return(out)
}

detect_restoration_categories <- function(df, aging_thresh = 0.05, rej_thresh = 0.05) {
  required_cols <- c("delta_old", "delta_rej")
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns: delta_old or delta_rej")
  }
  df <- df %>%
    mutate(
      aging_direction = case_when(
        delta_old >  aging_thresh ~ "hyper_in_aging",
        delta_old < -aging_thresh ~ "hypo_in_aging",
        TRUE ~ "no_change"
      ),
      restored_simple = case_when(
        aging_direction == "hyper_in_aging" & delta_rej < -rej_thresh ~ "restored",
        aging_direction == "hypo_in_aging"  & delta_rej >  rej_thresh ~ "restored",
        aging_direction %in% c("hyper_in_aging","hypo_in_aging") ~ "not_restored",
        TRUE ~ "no_change"
      ),
      restoration_category = paste0(aging_direction, "_", restored_simple)
    )
  return(df)
}

rej_thresh <- 0.05
global_sig_mat <- methylation_matrix[sig_global$cpg_id, , drop = FALSE]
rest_global_sig <- detect_restoration(global_sig_mat, metadata, rownames(global_sig_mat)) %>%
  detect_restoration_categories(aging_thresh = 0.00, rej_thresh = rej_thresh)

promoter_sig_mat <- meth_promoter_avg[sig_promoters$promoter_id, , drop = FALSE]
rest_promoters_sig <- detect_restoration(promoter_sig_mat, metadata, rownames(promoter_sig_mat)) %>%
  detect_restoration_categories(aging_thresh = 0.00, rej_thresh = rej_thresh)

gene_sig_mat <- meth_gene_avg[sig_genes$gene_symbol, , drop = FALSE]
rest_genes_sig <- detect_restoration(gene_sig_mat, metadata, rownames(gene_sig_mat)) %>%
  detect_restoration_categories(aging_thresh = 0.00, rej_thresh = rej_thresh)

enh_sig_mat <- meth_enhancer_avg[sig_enhancers$enhancer_id, , drop = FALSE]
rest_enhancers_sig <- detect_restoration(enh_sig_mat, metadata, rownames(enh_sig_mat)) %>%
  detect_restoration_categories(aging_thresh = 0.00, rej_thresh = rej_thresh)

all_rest <- bind_rows(
  rest_global_sig    %>% mutate(feature_type = "GlobalCpG"),
  rest_promoters_sig %>% mutate(feature_type = "Promoter"),
  rest_genes_sig     %>% mutate(feature_type = "GeneBody"),
  rest_enhancers_sig %>% mutate(feature_type = "Enhancer")
)

summary_table <- all_rest %>%
  group_by(feature_type, restoration_category) %>%
  summarise(n_features = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = restoration_category, values_from = n_features, values_fill = 0)
print(summary_table)
