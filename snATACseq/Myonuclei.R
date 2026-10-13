library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2016)
library(dplyr)
library(tibble)
library(ggplot2)
library(viridis)
library(metap)
library(writexl)
library(viridis)
library(pheatmap)
library(ComplexHeatmap)
library(reticulate)
library(phateR)
library(readr)
library(openxlsx)
library(Rmagic)
library(devtools)
library(monocle3)

ATAC_Myonuclei_ <- "/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_Myonuclei.rds"
Myonuclei <- readRDS(ATAC_Myonuclei_)

Cluster_colors <-c('coral2', 'burlywood', 'darkseagreen')
plot_newclusters_default_Myonuclei <- DimPlot(object = Myonuclei, label = FALSE, pt.size = 1.2, cols = Cluster_colors) + NoLegend()
plot_newclusters_condition_Myonuclei <- DimPlot(object = Myonuclei, label = FALSE, pt.size = 1.2, , group.by = 'sample_ID', cols = c('GFP' = 'grey', 'TEVp' = 'deepskyblue2')) 
plot_newclusters_default_Myonuclei | plot_newclusters_condition_Myonuclei

gfp_myonuclei <- WhichCells(Myonuclei, expression = sample_ID == "GFP")
tevp_myonuclei <- WhichCells(Myonuclei, expression = sample_ID == "TEVp")

gfp_myonuclei_plot <- DimPlot(Myonuclei, reduction = "umap", cells = gfp_myonuclei, pt.size = 2.0, cols = Cluster_colors) + NoLegend()
tevp_myonuclei_plot <- DimPlot(Myonuclei, reduction = "umap", cells = tevp_myonuclei, pt.size = 2.0, cols = Cluster_colors)
gfp_myonuclei_plot | tevp_myonuclei_plot