library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(ggplot2)
library(magick)
library(viridis)
library(metap)
library(writexl)
library(viridis)
library(pheatmap)
library(ComplexHeatmap)
library(reticulate)
library(phateR)
library(readr)
library(Rmagic)
library(devtools)
library(monocle3)
library(cicero)
library(SeuratWrappers)
library(gplots)
library(tidyr)
set.seed(123)

ATAC_Myonuclei_chromVar_ <- "/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_Myonuclei_chromVar.rds"
ATAC_Myonuclei_ <- "/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_Myonuclei.rds"
ATAC_Myonuclei_chromVar <- readRDS(ATAC_Myonuclei_chromVar_)
ATAC_Myonuclei <- readRDS(ATAC_Myonuclei_)

input_cds <- as.cell_data_set(ATAC_Myonuclei_chromVar)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, num_dim =50)
input_cds <- cluster_cells(input_cds,resolution = 0.0001)
input_cds <- learn_graph(input_cds)
input_cds <- order_cells(input_cds) 
pseudotime_values <- pseudotime(input_cds)

pseudo <- plot_cells(input_cds,
          color_cells_by = "pseudotime",
          cell_size = 2.0,
          cell_stroke = 0.,
          alpha = 0.8,
          label_branch_points = FALSE,
          label_roots = FALSE,
          label_leaves = FALSE,
          genes = NULL,
          show_trajectory_graph = TRUE,
          trajectory_graph_color = "white",
          trajectory_graph_segment_size = 1.
)

custom_palette <- viridis(100, option = "magma", direction = -1)  # Use the magma color palette from viridis
pseudo <- pseudo + scale_color_gradientn(colors = custom_palette)
pseudo_sample_ID <- plot_cells(input_cds,
                     color_cells_by = "sample_ID",
                     cell_size = 2.0,
                     cell_stroke = 0.,
                     alpha = 0.8,
                     label_branch_points = FALSE,
                     label_roots = FALSE,
                     label_leaves = FALSE,
                     genes = NULL,
                     show_trajectory_graph = TRUE,
                     trajectory_graph_color = "white",
                     trajectory_graph_segment_size = 1.,
)
pseudo_sample_ID <- pseudo_sample_ID + scale_color_manual(values = c("grey", "deepskyblue2")) + 
  ggplot2::facet_wrap(~sample_ID) 
pseudo | pseudo_sample_ID
pseudo_sample_ID

DefaultAssay(ATAC_Myonuclei) <- "RNA"
Sarcomeric_gene_set <- c(
  'Myl1', 'Myh1', 'Acta1',
  'Myom1', 'Tpm2', 'Myl3', 'Tnni2', 'Actn2', 'Tnnt3', 'Myh2',
  'Ttn', 'Neb'
)
Dediff_gene_set <- c(
  'Myh3', 'Myh8', 'Tnnt2'
)
gene <- as.data.frame(t(ATAC_Myonuclei@assays$RNA$data))
gene <- gene[Sarcomeric_gene_set]
gene_df <- data.frame(Pseudotime = pseudotime_values, gene)
gene_df <- gene_df[order(gene_df$Pseudotime), ]
gene_df <- as.matrix(t(gene_df[, -1]))
gene_df <- t(apply(gene_df, 1, function(x) smooth.spline(x, df = 3)$y))
gene_df <- t(apply(gene_df, 1, function(x) (x - mean(x)) / sd(x)))
gene_peak_time <- apply(gene_df, 1, which.max)
gene_order <- order(gene_peak_time)
gene_df <- gene_df[gene_order, ]

myColor <- colorRampPalette(c('steelblue3', "azure", 'darkorange'))(100)
pheatmap(
  gene_df,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = myColor,
  main = "Sarcomeric Gene Accessibility Along Pseudotime",
  use_raster = FALSE
)
write.csv(gene_df, "Sarcomeric_pseudotime.csv", row.names = FALSE)

diff_test_res <- graph_test(input_cds, neighbor_graph="principal_graph", cores=4)
sig_peaks <- subset(diff_test_res, q_value < 0.05)
ACRs_ <- as.data.frame(t(ATAC_Myonuclei_chromVar@assays$peaks$counts))

ACRs <- ACRs_[rownames(sig_peaks)]
ACRs <- t(scale(t(ACRs)))
ACRs_df <- data.frame(
  Pseudotime = pseudotime_values,
  ACRs
)
n_rows <- min(10000, nrow(ACRs))
random_rows <- sample(1:nrow(ACRs), n_rows)
n_cols <- min(2500, ncol(ACRs))
random_columns <- sample(1:ncol(ACRs), n_cols)
ACRs_df <- ACRs_df[random_rows, c(1, random_columns)]
ACRs_df <- ACRs_df[order(ACRs_df$Pseudotime), ] 
ACRs_df <- as.matrix(t(ACRs_df[, -1]))
ACRs_df <- t(apply(ACRs_df,1,function(x){smooth.spline(x,df=3)$y}))
ACRs_df <- t(apply(ACRs_df,1,function(x){(x-mean(x))/sd(x)}))
t0_values <- ACRs_df[, 2]
t0_values_ord <- order(t0_values, decreasing = TRUE)
ACRs_df <- ACRs_df[t0_values_ord, ]  
ACRs_df__ <- t(apply(ACRs_df,2,function(x){smooth.spline(x,df=3)$y}))
ACRs_df__ <- t(apply(ACRs_df__,2,function(x){(x-mean(x))/sd(x)}))


myColor <- colorRampPalette(c('gray60', "azure", 'steelblue3'))(100)
pheatmap(
  ACRs_df__, 
  cluster_cols = FALSE,     
  cluster_rows = FALSE,   
  scale = "row",         
  show_rownames = FALSE,   
  show_colnames = FALSE, 
  color = myColor, 
  main = "Motif Activity Along Pseudotime",
  use_raster = FALSE
)
ACRs_mean_pseudotime <- colMeans(ACRs_df)
plot(ACRs_mean_pseudotime, type = "l", col = "gray50", lwd = 2,
     xlab = "Pseudotime", ylab = "Average Value", 
     main = "Average ACRs")
grid()

mean_by_sign <- function(mat, sign = c("positive", "negative")) {
  sign <- match.arg(sign)
  
  apply(mat, 2, function(x) {
    if(sign == "positive") {
      vals <- x[x > 0]
    } else {
      vals <- x[x < 0]
    }
    if(length(vals) == 0) return(NA)
    mean(vals)
  })
}
ACRs_mean_pos <- mean_by_sign(ACRs_df, "positive")
ACRs_mean_neg <- mean_by_sign(ACRs_df, "negative")
y_min <- min(c(ACRs_mean_pos, ACRs_mean_neg), na.rm = TRUE)
y_max <- max(c(ACRs_mean_pos, ACRs_mean_neg), na.rm = TRUE)
plot(ACRs_mean_pos, type = "l", col = "steelblue3", lwd = 2,
     xlab = "Pseudotime", ylab = "Mean Value",
     main = "Positive vs Negative ACRs Along Pseudotime",
     ylim = c(y_min, y_max))

lines(ACRs_mean_neg, col = "gray50", lwd = 2)
grid()

motif_list <- c("MA0499.2", "MA1641.1", "MA0667.1","MA0500.2")
motif <- as.data.frame(t(ATAC_Myonuclei_chromVar@assays$chromvar@data)) 
motif <- motif[motif_list] 
motif_df <- data.frame(
  Pseudotime = pseudotime_values,
  motif
)
motif_df <- motif_df[order(motif_df$Pseudotime), ]
motif_df <- as.matrix(t(motif_df[, -1])) 
motif_df <- t(apply(motif_df,1,function(x){smooth.spline(x,df=3)$y}))
motif_df <- t(apply(motif_df,1,function(x){(x-mean(x))/sd(x)}))
write.csv(motif_df, "Motifs_MRFs_pseudotime.csv", row.names = FALSE)

pheatmap(
  motif_df, 
  cluster_cols = FALSE,      
  cluster_rows = FALSE,     
  scale = "row",           
  show_rownames = TRUE,   
  show_colnames = FALSE,    
  color = myColor, 
  main = "Motif Activity Along Pseudotime",
)
