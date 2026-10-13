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
library(metap)
library(writexl)
library(viridis)
library(pheatmap)
library(reticulate)
library(phateR)
library(readr)
library(Rmagic)
library(openxlsx)

ATAC_MuSCs_ <- "/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_MuSCs.rds"
MuSCs <- readRDS(ATAC_MuSCs_)

Cluster_colors <-c('coral2', 'skyblue1')
plot_clusters_default_MuSCs<- DimPlot(object = MuSCs, label = TRUE, pt.size = 3.0, cols = Cluster_colors) + NoLegend()
plot_clusters_condition_MuSCs <- DimPlot(object = MuSCs, label = FALSE, pt.size = 3.0, group.by = 'sample_ID', cols = c('GFP' = 'grey', 'TEVp' = 'deepskyblue2')) + ggtitle('MuSCs') + theme(legend.position = c(0.6, 0.1))
plot_clusters_default_MuSCs
plot_clusters_default_MuSCs | plot_clusters_condition_MuSCs
gfp_MuSCs <- WhichCells(MuSCs, expression = sample_ID == "GFP")
tevp_MuSCs <- WhichCells(MuSCs, expression = sample_ID == "TEVp")
gfp_MuSCs_plot <- DimPlot(MuSCs, cells = gfp_MuSCs, pt.size = 3.0, cols = Cluster_colors)
tevp_MuSCs_plot <- DimPlot(MuSCs, cells = tevp_MuSCs, pt.size = 3.0, cols = Cluster_colors)
gfp_MuSCs_plot
tevp_MuSCs_plot
gfp_MuSCs_plot | tevp_MuSCs_plot


cluster_counts <- table(MuSCs@meta.data$seurat_clusters, MuSCs@meta.data$sample_ID)
cluster_df <- as.data.frame(cluster_counts)
colnames(cluster_df) <- c("seurat_cluster", "sample_ID", "count")
cluster_df <- cluster_df %>%
  group_by(sample_ID) %>%
  mutate(percentage = count / sum(count) * 100)
ggplot(cluster_df, aes(x = "", y = percentage, fill = as.factor(seurat_cluster))) +
  geom_bar(stat = "identity", width = 10) +
  coord_polar("y") +
  facet_wrap(~ sample_ID) +  
  labs(fill = "Cluster", y = "Percentage", x = "") +
  theme_void() +  
  scale_fill_manual(values = Cluster_colors)
ggtitle("-")
