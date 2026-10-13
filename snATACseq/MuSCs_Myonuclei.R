library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)

MuSCs <- readRDS("/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_MuSCs.rds")
Myonuclei <- readRDS("/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_Myonuclei.rds")

MuSCs$dataset <- "MuSCs"
Myonuclei$dataset <- "Myonuclei"

MuSCs$combined_ident <- paste0("MuSCs_", Idents(MuSCs))
Myonuclei$combined_ident <- paste0("Myonuclei_", Idents(Myonuclei))

combined <- merge(MuSCs, y = Myonuclei, add.cell.ids = c("MuSCs", "Myonuclei"))
Idents(combined) <- combined$combined_ident
DefaultAssay(combined) <- "RNA"
markers <- c("Pax7", "Calcr", "Itga7", "Myod1", "Myog", "Myf6")

DotPlot(
  object = combined,
  features = markers,
  dot.scale = 12,
  cols = "RdBu"
) +
  scale_color_viridis(option = "mako", direction = -1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    x = "Marker genes",
    y = "Clusters",
    title = "Expression of muscle-related genes across MuSC and Myonuclei clusters"
  )
