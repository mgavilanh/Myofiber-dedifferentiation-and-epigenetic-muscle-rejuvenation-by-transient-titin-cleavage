## Motif enrichment analysis.
library(Seurat) 
library(Signac) 
library(pheatmap)
library(openxlsx)
library(here)
library(patchwork)
set.seed(1234)

ATAC_Myonuclei_chromVar_ <- "/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_Myonuclei_chromVar.rds"
motif_ <- "/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_Myonuclei_motif.rds"

ATAC_Myonuclei_chromVar <- readRDS(ATAC_Myonuclei_chromVar_)
motif <- readRDS(motif_)


GetMotifs <- function(cluster, seurat_aggregate) {
  print(paste0("Finding motifs for: ",cluster))
  atac <- seurat_aggregate
  DefaultAssay(atac) <- 'peaks'
  dac <- FindMarkers(atac,
                     ident.1 = cluster,
                     test.use = 'LR',
                     latent.vars = "nCount_peaks") # find all cluster-specific degs
  enriched.motifs <- FindMotifs(object = atac, features = rownames(dac[dac$p_val < 0.05, ]))
  return(enriched.motifs)
}

GetChromvarActivities <- function(cluster, seurat_aggregate, motif) {
  print(paste0("Finding chromVAR activities for: ",cluster))
  atac <- seurat_aggregate
  DefaultAssay(atac) <- 'chromvar'
  dam <- FindMarkers(atac, 
                     ident.1 = cluster,
                     test.use = 'LR',
                     latent.vars = "nCount_peaks",
                     logfc.threshold = 0,) # find all cluster-specific degs
  motifLookup <- rownames(dam)
  motifNames <- sapply(motifLookup, function(x) motif@motif.names[[x]])
  return(cbind(dam, gene = motifNames))
}

idents <- levels(ATAC_Myonuclei_chromVar)
list.cluster.dac <- lapply(idents, function(x) GetMotifs(x, seurat_aggregate = ATAC_Myonuclei_chromVar))
write.xlsx(list.cluster.dac, file = "atac_myonuclei_cluster_motifs.xlsx", sheetName = idents, rowNames = T)

list.cluster.dam <- lapply(idents, function(x) GetChromvarActivities(x, seurat_aggregate = ATAC_Myonuclei_chromVar, motif = motif))
write.xlsx(list.cluster.dam, file = "atac_myonuclei_cluster_motifs_chromVAR.xlsx", sheetName = idents, rowNames = T)

Selected_motifs<-c("MA0139.1",
                   "MA1102.2",
                   "MA0117.2",
                   "MA0842.2",
                   "MA1568.1",
                   "MA0495.3",
                   "MA0751.1",
                   "MA0696.1",
                   "MA1596.1",
                   "MA1653.1",
                   "MA1640.1",
                   "MA1584.1",
                   "MA1625.1",
                   "MA0507.1",
                   "MA0757.1",
                   "MA0592.3",
                   "MA0643.1",
                   "MA0141.3",
                   "MA1141.1",
                   "MA1134.1",
                   "MA1137.1",
                   "MA1130.1",
                   "MA0491.2",
                   "MA0478.1",
                   "MA1128.1",
                   "MA1144.1",
                   "MA0477.2",
                   "MA0490.2",
                   "MA1142.1",
                   "MA1135.1",
                   "MA0521.1",
                   "MA0816.1",
                   "MA0500.2",
                   "MA1641.1",
                   "MA1467.1",
                   "MA1619.1",
                   "MA0832.1",
                   "MA0665.1",
                   "MA1100.2",
                   "MA0667.1",
                   "MA0048.2",
                   "MA0691.1",
                   "MA0499.2",
                   "MA1635.1",
                   "MA0091.1")

aver_chromvar <- AverageExpression(ATAC_Myonuclei_chromVar, assays = "chromvar", features = Selected_motifs)

str(aver_chromvar)
aver_chromvar_matrix <- aver_chromvar$chromvar
ordered_aver_chromvar <- aver_chromvar_matrix[order(aver_chromvar_matrix[,1], decreasing = TRUE), ]
head(ordered_aver_chromvar)
aver_chromvar <- as.data.frame(ordered_aver_chromvar) 
aver_chromvar$max <- max.col(aver_chromvar)
aver_chromvar <- aver_chromvar[order(aver_chromvar$max), ]
aver_chromvar <- dplyr::select(aver_chromvar, -max)
colnames(aver_chromvar) <- levels(Idents(ATAC_Myonuclei_chromVar))
head(aver_chromvar)

aver_chromvar$max <- max.col(aver_chromvar)
aver_chromvar <- aver_chromvar[order(aver_chromvar$max), ]
aver_chromvar <- dplyr::select(aver_chromvar, -max)
colnames(aver_chromvar) <- levels(Idents(ATAC_Myonuclei_chromVar))

aver_chromvar_matrix <- as.matrix(aver_chromvar)


paletteLength <- 50
myColor <- colorRampPalette(c('gray60', "azure", 'steelblue3'))(100)
myColor <- viridis(100, option = "mako", direction = -1)  # Use the magma color palette from viridis

TF_by_cluster<-pheatmap(aver_chromvar_matrix,scale = "row", cluster_cols=F,cluster_rows = F,show_rownames=T,color = myColor)

TF_by_cluster

DefaultAssay(ATAC_Myonuclei_chromVar) <- 'chromvar'


TF_Myod_feature <- FeaturePlot(
  object = ATAC_Myonuclei_chromVar,
  features = c("MA0499.2"),
  cols = c('lightgrey', 'blue'),
  split.by = 'sample_ID',
  pt.size = 1.,
  alpha = 0.7,
  min.cutoff = 'q10',
  max.cutoff = 'q90'
)

TF_Myog_feature <- FeaturePlot(
  object = ATAC_Myonuclei_chromVar,
  features = c("MA0500.2"),
  cols = c('lightgrey', 'blue'),
  split.by = 'sample_ID',
  pt.size = 1.,
  alpha = 0.7,
  min.cutoff = 'q10',
  max.cutoff = 'q90'
)

TF_Myf5_feature <- FeaturePlot(
  object = ATAC_Myonuclei_chromVar,
  features = c("MA1641.1"),
  cols = c('lightgrey', 'blue'),
  split.by = 'sample_ID',
  pt.size = 1.,
  alpha = 0.7,
  min.cutoff = 'q10',
  max.cutoff = 'q90'
)

TF_Myf6_feature <- FeaturePlot(
  object = ATAC_Myonuclei_chromVar,
  features = c("MA0667.1"),
  cols = c('lightgrey', 'blue'),
  split.by = 'sample_ID',
  pt.size = 1.,
  alpha = 0.7,
  min.cutoff = 'q10',
  max.cutoff = 'q90'
)

(TF_Myod_feature | TF_Myog_feature) / (TF_Myf5_feature | TF_Myf6_feature)