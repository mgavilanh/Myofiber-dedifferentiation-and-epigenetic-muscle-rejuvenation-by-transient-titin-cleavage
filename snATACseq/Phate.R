library(Seurat)
library(phateR)
library(reticulate)
library(ggplot2)
library(ggridges)
library(dplyr)

setwd("/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/")

MuSCs_Myonuclei <- readRDS("/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_MuSCs_Myonuclei.rds")

DefaultAssay(MuSCs_Myonuclei) <- "RNA"
rna_mat <- t(as.matrix(GetAssayData(MuSCs_Myonuclei, slot = "data")))
set.seed(42)
phate_obj <- phate(rna_mat, ndim = 2, knn = 10, t = 40)

phate_embed <- as.matrix(phate_obj)
colnames(phate_embed) <- c("PHATE1","PHATE2")

MuSCs_Myonuclei[["phate"]] <- CreateDimReducObject(
  embeddings = phate_embed,
  key = "PHATE_",
  assay = "RNA"
)

DimPlot(
  MuSCs_Myonuclei,
  reduction = "phate",
  label = TRUE,
  repel = TRUE,
  pt.size = 1
) + ggplot2::ggtitle("PHATE (colored by Idents)")

phate_df <- as.data.frame(Embeddings(MuSCs_Myonuclei, "phate"))
print(colnames(phate_df))
phate_df$Ident <- Idents(MuSCs_Myonuclei)
phate_export <- phate_df[, c(colnames(phate_df)[1], "Ident")]
write.csv(phate_export, "/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/excels/PHATE1_MuSCs_Myonuclei.csv", row.names = TRUE)

phate1_col <- colnames(phate_df)[1]
ggplot(phate_df, aes_string(x = phate1_col, fill = "Ident")) +
  geom_density(alpha = 0.4) +
  theme_classic() +
  ggtitle("Distribution of PHATE1 per Ident") +
  xlab("PHATE1") +
  ylab("Density")

ggplot(phate_df, aes(x = Ident, y = .data[[phate1_col]], fill = Ident)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  ggtitle("PHATE1 Distribution per Ident") +
  ylab("PHATE1")

ggplot(phate_df, aes(x = .data[[phate1_col]], y = Ident, fill = Ident)) +
  geom_density_ridges(alpha = 0.6) +
  theme_classic() +
  ggtitle("PHATE1 Ridgeline Distribution per Ident") +
  xlab("PHATE1")

analyze_phate_threshold <- function(phate_df, phate1_col, threshold, all_levels) {
  
  make_subset <- function(data, phate1_col, condition_fn, all_levels) {
    df <- data %>% filter(condition_fn(.data[[phate1_col]]))

    missing <- setdiff(all_levels, unique(df$Ident))
    if(length(missing) > 0){
      dummy <- data.frame(matrix(NA, nrow = length(missing), 
                                 ncol = ncol(df)-1))
      colnames(dummy) <- colnames(df)[1:(ncol(df)-1)]
      dummy$Ident <- factor(missing, levels = all_levels)
      df <- bind_rows(df, dummy)
    }
    
    df$Ident <- factor(df$Ident, levels = all_levels)
    df
  }
  
  subset_low  <- make_subset(phate_df, phate1_col, function(x) x < threshold, all_levels)
  subset_high <- make_subset(phate_df, phate1_col, function(x) x > threshold, all_levels)
  
  total_counts <- phate_df %>%
    group_by(Ident) %>% summarise(total_cells = n(), .groups = "drop")
  
  pct_tbl <- function(subset_df) {
    subset_df %>%
      group_by(Ident) %>% summarise(subset_cells = n(), .groups = "drop") %>%
      right_join(total_counts, by = "Ident") %>%
      mutate(subset_cells = ifelse(is.na(subset_cells), 0, subset_cells),
             percent = 100 * subset_cells / total_cells)
  }
  
  percent_low  <- pct_tbl(subset_low)
  percent_high <- pct_tbl(subset_high)
  
  combined <- percent_low %>%
    rename(subset_cells_low = subset_cells, percent_low = percent) %>%
    left_join(
      percent_high %>% 
        rename(subset_cells_high = subset_cells, percent_high = percent),
      by = "Ident"
    )
  
  plot_density <- function(df, title) {
    ggplot(df, aes_string(x = phate1_col, fill = "Ident")) +
      geom_density(alpha = 0.4, na.rm = TRUE) +
      theme_classic() + ggtitle(title)
  }
  
  plot_violin <- function(df, title) {
    ggplot(df, aes(x = Ident, y = .data[[phate1_col]], fill = Ident)) +
      geom_violin(trim = TRUE, na.rm = TRUE) +
      theme_classic() + ggtitle(title)
  }
  
  plot_ridge <- function(df, title) {
    ggplot(df, aes(x = .data[[phate1_col]], y = Ident, fill = Ident)) +
      geom_density_ridges(alpha = 0.6, na.rm = TRUE) +
      theme_classic() + ggtitle(title)
  }
  
  list(
    subset_low   = subset_low,
    subset_high  = subset_high,
    percent_table = combined,
    density_low  = plot_density(subset_low,  "PHATE1 < threshold"),
    density_high = plot_density(subset_high, "PHATE1 > threshold"),
    violin_low   = plot_violin(subset_low,   "PHATE1 < threshold"),
    violin_high  = plot_violin(subset_high,  "PHATE1 > threshold"),
    ridge_low    = plot_ridge(subset_low,    "PHATE1 < threshold"),
    ridge_high   = plot_ridge(subset_high,   "PHATE1 > threshold")
  )
}

threshold = -0.025

result <- analyze_phate_threshold(
  phate_df,
  phate1_col = phate1_col,
  threshold  = threshold,
  all_levels = levels(Idents(MuSCs_Myonuclei))
)

result$percent_table
result$density_low
result$density_high
result$violin_low
result$violin_high
result$ridge_low
result$ridge_high

combined_violin_df <- bind_rows(
  result$subset_low  %>% filter(!is.na(.data[[phate1_col]])) %>% mutate(Group = paste0("PHATE1 < ", threshold)),
  result$subset_high %>% filter(!is.na(.data[[phate1_col]])) %>% mutate(Group = paste0("PHATE1 > ", threshold))
)

ggplot(combined_violin_df,
       aes(x = Ident, y = .data[[phate1_col]], fill = Ident)) +
  geom_violin(trim = TRUE, na.rm = TRUE, scale = "width") +
  facet_wrap(~ Group, ncol = 2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("PHATE1 Distribution Below and Above Threshold") +
  ylab("PHATE1")
