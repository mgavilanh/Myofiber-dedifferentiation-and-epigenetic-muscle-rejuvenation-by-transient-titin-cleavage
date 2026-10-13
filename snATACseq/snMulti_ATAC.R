library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2016)
library(dplyr)
library(tidyr)
library(ggplot2)
library(metap)
library(writexl)
library(viridis)
library(pheatmap)
library(phateR)
library(readr)
library(Rmagic)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

setwd("/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics")

ATAC_counts <- Read10X_h5(filename = "./ATAC_results/filtered_peak_bc_matrix.h5")
ATAC_chrom_assay <- CreateChromatinAssay(
  counts = ATAC_counts,
  sep = c(":", "-"),
  fragments = './ATAC_results/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
ATAC_metadata <- read.csv(
  file = "./ATAC_results/singlecell.csv",
  header = TRUE,
  row.names = 1
)
ATAC <- CreateSeuratObject(
  counts = ATAC_chrom_assay,
  assay = "peaks",
  meta.data = ATAC_metadata
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"
Annotation(ATAC) <- annotations
ATAC <- NucleosomeSignal(ATAC)
ATAC <- TSSEnrichment(object = ATAC, fast = FALSE)
ATAC$blacklist_ratio <- ATAC$blacklist_region_fragments / ATAC$peak_region_fragments * 100
ATAC$pct_reads_in_peaks <- ATAC$peak_region_fragments / ATAC$passed_filters * 100
VlnPlot(object = ATAC,
        features = c('nCount_peaks', 'nFeature_peaks'),
        pt.size = 0.1,
        ncol = 2)
ATAC <- subset(
  x = ATAC,
    nCount_peaks > 500 &
    nCount_peaks < 70000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
VlnPlot(object = ATAC,
        features = c('nCount_peaks', 'nFeature_peaks'),
        pt.size = 0.1,
        ncol = 2)
ATAC@meta.data <- ATAC@meta.data %>%
  mutate(sample_ID = case_when(
    grepl("-1$", row.names(ATAC@meta.data)) ~ "GFP",
    grepl("-2$", row.names(ATAC@meta.data)) ~ "TEVp",
    TRUE ~ NA_character_ # Optionally handle other cases
  ))

ATAC <- RunTFIDF(ATAC)
ATAC <- FindTopFeatures(ATAC, min.cutoff = 'q0')
ATAC <- RunSVD(ATAC)
ATAC <- RunUMAP(object = ATAC, reduction = 'lsi', dims = 2:30)
ATAC <- FindNeighbors(object = ATAC, reduction = 'lsi', dims = 2:30)
ATAC <- FindClusters(object = ATAC, verbose = FALSE, algorithm = 3, resolution = 0.5)
plot_clusters_default <- DimPlot(object = ATAC, label = TRUE, pt.size = 0.7) + NoLegend()
plot_clusters_condition <- DimPlot(object = ATAC, label = FALSE, pt.size = 0.7, group.by = 'sample_ID', cols = c('GFP' = 'grey', 'TEVp' = 'deepskyblue2')) + ggtitle('') + theme(legend.position = c(0.6, 0.1))
plot_clusters_default | plot_clusters_condition
plot_clusters_default
plot_clusters_condition

new.cluster.default <- c("Myonuclei", 
                     "Myonuclei",
                     "FAPs",
                     "Myonuclei",
                     "Myeloid cell",
                     "FAPs",
                     "FAPs", 
                     "Smooth muscle",
                     "MuSCs",
                     "Endothelial",
                     "Lymphocyte",
                     "Myonuclei"
)
names(new.cluster.default) <- levels(ATAC)
ATAC <- RenameIdents(ATAC, new.cluster.default)
plot_newclusters_default <- DimPlot(object = ATAC, label = TRUE, pt.size = 0.7) + NoLegend()
plot_newclusters_default | plot_clusters_condition
saveRDS(ATAC, file='/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC.rds')
ATAC_gene.activities <- GeneActivity(ATAC)
ATAC[['RNA']] <- CreateAssayObject(counts = ATAC_gene.activities)
ATAC <- NormalizeData(
  object = ATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC$nCount_RNA)
)

cells_per_cluster_sample <- table(cluster = Idents(ATAC), sample = ATAC$sample_ID)
frags <- ATAC$peak_region_fragments
frag_summary <- aggregate(frags, by = list(cluster = as.character(Idents(ATAC)), sample = ATAC$sample_ID), FUN = median)
colnames(frag_summary)[3] <- "median_fragments"
peaks <- GetAssayData(ATAC, assay = "peaks", slot = "counts")
peaks_per_cell <- colSums(peaks > 0)
ATAC$peaks_per_cell <- peaks_per_cell
meta <- data.frame(
  cell = colnames(ATAC),
  cluster = Idents(ATAC),
  sample = ATAC$sample_ID,
  peaks_per_cell = peaks_per_cell
)
peaks_per_cell_summary <- meta %>%
  group_by(cluster, sample) %>%
  summarize(
    n_cells = n(),
    median_peaks_per_cell = median(peaks_per_cell),
    mean_peaks_per_cell = mean(peaks_per_cell)
  ) %>%
  arrange(cluster, sample)

ggplot(peaks_per_cell_summary,
       aes(x = cluster,
           y = mean_peaks_per_cell,
           fill = sample)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.75)) +
  labs(
    x = "Cluster",
    y = "Mean peaks per nucleus",
    fill = "Condition"
  ) +
  scale_fill_manual(values = c("GFP" = "grey", "TEVp" = "deepskyblue2")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak_mat <- GetAssayData(ATAC, assay="peaks", slot="counts")
thresh <- 0.05
cells_gfp  <- WhichCells(ATAC, expression = sample_ID == "GFP")
cells_tevp <- WhichCells(ATAC, expression = sample_ID == "TEVp")
anno_list <- list()
if (length(cells_gfp) > 20) {
  peaks_gfp <- granges(ATAC)[
    rowSums(peak_mat[, cells_gfp] > 0) / length(cells_gfp) > thresh
  ]
  anno_list[["GFP"]] <- annotatePeak(
    peaks_gfp,
    TxDb = txdb,
    tssRegion = c(-2000, 2000),
    verbose = FALSE
  )
}
if (length(cells_tevp) > 20) {
  peaks_tevp <- granges(ATAC)[
    rowSums(peak_mat[, cells_tevp] > 0) / length(cells_tevp) > thresh
  ]
  anno_list[["TEVp"]] <- annotatePeak(
    peaks_tevp,
    TxDb = txdb,
    tssRegion = c(-2000, 2000),
    verbose = FALSE
  )
}
plotAnnoBar(anno_list)
anno_summary <- lapply(names(anno_list), function(cond) {
  df <- anno_list[[cond]]@annoStat
  df$Condition <- cond
    if ("Feature" %in% colnames(df)) {
    colnames(df)[colnames(df) == "Feature"] <- "Annotation"
  } else if ("Category" %in% colnames(df)) {
    colnames(df)[colnames(df) == "Category"] <- "Annotation"
  }
  df
}) %>%
  bind_rows()
anno_table <- anno_summary %>%
  select(Condition, Annotation, Frequency) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = Frequency, values_fill = 0) %>%
  arrange(desc(rowMeans(select(., -Annotation)))) %>%
  mutate(across(where(is.numeric), round, 2))
anno_table

anno_list <- list()
for(cl in levels(ATAC)){
  cells_cl <- WhichCells(ATAC, idents = cl)
  cells_gfp  <- cells_cl[ ATAC$sample_ID[cells_cl]  == "GFP" ]
  cells_tevp <- cells_cl[ ATAC$sample_ID[cells_cl] == "TEVp" ]
  if(length(cells_gfp) > 20){
    peaks_gfp <- granges(ATAC)[ rowSums( peak_mat[,cells_gfp] > 0 )/length(cells_gfp) > thresh ]
    anno_list[[ paste0(cl,"_GFP") ]] <- annotatePeak(peaks_gfp,
                                                     TxDb = txdb,
                                                     tssRegion = c(-2000,2000))
  }
  if(length(cells_tevp) > 20){
    peaks_tevp <- granges(ATAC)[ rowSums( peak_mat[,cells_tevp] > 0 )/length(cells_tevp) > thresh ]
    anno_list[[ paste0(cl,"_TEVp") ]] <- annotatePeak(peaks_tevp,
                                                      TxDb = txdb,
                                                      tssRegion = c(-2000,2000))
  }
}
plotAnnoBar(anno_list)

ATAC$cluster_sample <- paste0(Idents(ATAC), "_", ATAC$sample_ID)
plot1_MuSCs <- CoveragePlot(
  object = ATAC,
  region = "Pax7",
  features = "Pax7",
  group.by = "cluster_sample",
  idents = c("Myonuclei_GFP", "Myonuclei_TEVp", "FAPs_GFP", "FAPs_TEVp", "Myeloid cell_GFP", "Myeloid cell_TEVp", "Smooth muscle_GFP", "Smooth muscle_TEVp", "MuSCs_GFP", "MuSCs_TEVp", "Endothelial_GFP", "Endothelial_TEVp", "Lymphocyte_GFP", "Lymphocyte_TEVp"),
  annotation = TRUE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
plot1_MuSCs <- plot1_MuSCs + theme_void() + NoLegend() + ggtitle('Pax7') + theme(plot.title = element_text(hjust = 0.5, vjust = 0.4))
plot1_MuSCs

plot2_Myonuclei <- CoveragePlot(
  object = ATAC,
  region = "Ttn",
  features = "Ttn",
  group.by = "cluster_sample",
  idents = c("Myonuclei_GFP", "Myonuclei_TEVp", "FAPs_GFP", "FAPs_TEVp", "Myeloid cell_GFP", "Myeloid cell_TEVp", "Smooth muscle_GFP", "Smooth muscle_TEVp", "MuSCs_GFP", "MuSCs_TEVp", "Endothelial_GFP", "Endothelial_TEVp", "Lymphocyte_GFP", "Lymphocyte_TEVp"),
  annotation = TRUE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
plot2_Myonuclei <- plot2_Myonuclei + theme_void() + NoLegend() + ggtitle('Ttn') + theme(plot.title = element_text(hjust = 0.5, vjust = 0.4))
plot2_Myonuclei

plot3_Endothelial <- CoveragePlot(
  object = ATAC,
  region = "Flt1",
  features = "Flt1",
  group.by = "cluster_sample",
  idents = c("Myonuclei_GFP", "Myonuclei_TEVp", "FAPs_GFP", "FAPs_TEVp", "Myeloid cell_GFP", "Myeloid cell_TEVp", "Smooth muscle_GFP", "Smooth muscle_TEVp", "MuSCs_GFP", "MuSCs_TEVp", "Endothelial_GFP", "Endothelial_TEVp", "Lymphocyte_GFP", "Lymphocyte_TEVp"),
  annotation = TRUE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
plot3_Endothelial <- plot3_Endothelial + theme_void() + NoLegend() + ggtitle('Flt1') + theme(plot.title = element_text(hjust = 0.5, vjust = 0.4))
plot3_Endothelial

plot4_SMC <- CoveragePlot(
  object = ATAC,
  region = "Myh11",
  features = "Myh11",
  group.by = "cluster_sample",
  idents = c("Myonuclei_GFP", "Myonuclei_TEVp", "FAPs_GFP", "FAPs_TEVp", "Myeloid cell_GFP", "Myeloid cell_TEVp", "Smooth muscle_GFP", "Smooth muscle_TEVp", "MuSCs_GFP", "MuSCs_TEVp", "Endothelial_GFP", "Endothelial_TEVp", "Lymphocyte_GFP", "Lymphocyte_TEVp"),
  annotation = TRUE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
plot4_SMC <- plot4_SMC + theme_void() + NoLegend() + ggtitle('Myh11') + theme(plot.title = element_text(hjust = 0.5, vjust = 0.4))
plot4_SMC

plot5_FAP <- CoveragePlot(
  object = ATAC,
  region = "Pdgfra",
  features = "Pdgfra",
  group.by = "cluster_sample",
  idents = c("Myonuclei_GFP", "Myonuclei_TEVp", "FAPs_GFP", "FAPs_TEVp", "Myeloid cell_GFP", "Myeloid cell_TEVp", "Smooth muscle_GFP", "Smooth muscle_TEVp", "MuSCs_GFP", "MuSCs_TEVp", "Endothelial_GFP", "Endothelial_TEVp", "Lymphocyte_GFP", "Lymphocyte_TEVp"),
  annotation = TRUE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
plot5_FAP <- plot5_FAP + theme_void() + NoLegend() + ggtitle('Pdgfra') + theme(plot.title = element_text(hjust = 0.5, vjust = 0.4))
plot5_FAP

plot6_Lymphocyte <- CoveragePlot(
  object = ATAC,
  region = "Cd3d",
  features = "Cd3d",
  group.by = "cluster_sample",
  idents = c("Myonuclei_GFP", "Myonuclei_TEVp", "FAPs_GFP", "FAPs_TEVp", "Myeloid cell_GFP", "Myeloid cell_TEVp", "Smooth muscle_GFP", "Smooth muscle_TEVp", "MuSCs_GFP", "MuSCs_TEVp", "Endothelial_GFP", "Endothelial_TEVp", "Lymphocyte_GFP", "Lymphocyte_TEVp"),
  annotation = TRUE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
plot6_Lymphocyte <- plot6_Lymphocyte + theme_void() + NoLegend() + ggtitle('Cd3d') + theme(plot.title = element_text(hjust = 0.5, vjust = 0.4))
plot6_Lymphocyte

plot7_Myeloid <- CoveragePlot(
  object = ATAC,
  region = "Mrc1",
  features = "Mrc1",
  group.by = "cluster_sample",
  idents = c("Myonuclei_GFP", "Myonuclei_TEVp", "FAPs_GFP", "FAPs_TEVp", "Myeloid cell_GFP", "Myeloid cell_TEVp", "Smooth muscle_GFP", "Smooth muscle_TEVp", "MuSCs_GFP", "MuSCs_TEVp", "Endothelial_GFP", "Endothelial_TEVp", "Lymphocyte_GFP", "Lymphocyte_TEVp"),
  annotation = TRUE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
plot7_Myeloid <- plot7_Myeloid + theme_void() + NoLegend() + ggtitle('Mrc1') + theme(plot.title = element_text(hjust = 0.5, vjust = 0.4))
plot7_Myeloid

MuSCs <- subset(ATAC, idents = c("MuSCs"))
Cluster_colors <-c('deepskyblue4', 'cadetblue', 'darkslategray3')

MuSCs <- RunTFIDF(MuSCs)
MuSCs <- FindTopFeatures(MuSCs, min.cutoff = 'q0')
MuSCs <- RunSVD(MuSCs)

MuSCs <- RunUMAP(MuSCs, reduction = "lsi", dims = 2:30)
MuSCs <- FindNeighbors(object = MuSCs, reduction = 'lsi', dims = 2:30)
MuSCs <- FindClusters(object = MuSCs, verbose = FALSE, algorithm = 3)
plot_clusters_default_MuSCs<- DimPlot(object = MuSCs, label = TRUE, pt.size = 2.5, cols = Cluster_colors) + NoLegend()
plot_clusters_condition_MuSCs <- DimPlot(object = MuSCs, label = FALSE, pt.size = 2.5, group.by = 'sample_ID', cols = c('GFP' = 'grey', 'TEVp' = 'deepskyblue2')) + ggtitle('MuSCs') + theme(legend.position = c(0.6, 0.1))
plot_clusters_default_MuSCs | plot_clusters_condition_MuSCs
plot_clusters_condition_MuSCs <- DimPlot(object = MuSCs, label = FALSE, pt.size = 2.5, group.by = 'sample_ID', cols = c('GFP' = 'grey', 'TEVp' = 'deepskyblue2')) + ggtitle('MuSCs') + theme(legend.position = c(0.6, 0.1))
new.cluster.default_MuSCs <- c("MuSCs_2", 
                               "MuSCs_1",
                               "MuSCs_1"
                               
)
names(new.cluster.default_MuSCs) <- levels(MuSCs)
MuSCs <- RenameIdents(MuSCs, new.cluster.default_MuSCs)
MuSCs@meta.data$seurat_clusters <- Idents(MuSCs)
plot_newclusters_default_MuSCs <- DimPlot(object = MuSCs, label = TRUE, pt.size = 2.5) + NoLegend()
plot_newclusters_default_MuSCs | plot_clusters_condition_MuSCs
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
saveRDS(MuSCs, file='/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_MuSCs.rds')

Myonuclei <- subset(ATAC, idents = c("Myonuclei"))
Myonuclei <- RunTFIDF(Myonuclei)
Myonuclei <- FindTopFeatures(Myonuclei, min.cutoff = 'q0')
Myonuclei <- RunSVD(Myonuclei)
Myonuclei <- RunUMAP(Myonuclei, reduction = "lsi", dims = 2:30)
Myonuclei <- FindNeighbors(object = Myonuclei, reduction = 'lsi', dims = 2:30)
Myonuclei <- FindClusters(object = Myonuclei, verbose = FALSE, algorithm = 3)
plot_clusters_default_Myonuclei<- DimPlot(object = Myonuclei, label = TRUE, pt.size = 1.2) + NoLegend()
plot_clusters_condition_Myonuclei <- DimPlot(object = Myonuclei, label = FALSE, pt.size = 1.2, group.by = 'sample_ID', cols = c('GFP' = 'grey', 'TEVp' = 'deepskyblue2')) + ggtitle('Myonuclei') + theme(legend.position = c(0.6, 0.1))
plot_clusters_default_Myonuclei | plot_clusters_condition_Myonuclei
new.cluster.default_Myonuclei <- c("Myonuclei_1", 
                         "Myonuclei_2",
                         "Myonuclei_3",
                         "Myonuclei_2",
                         "na",
                         "na",
                         "na",
                         "na"
                         
)
names(new.cluster.default_Myonuclei) <- levels(Myonuclei)
Myonuclei <- RenameIdents(Myonuclei, new.cluster.default_Myonuclei)
plot_newclusters_default_Myonuclei <- DimPlot(object = Myonuclei, label = TRUE, pt.size = 1.2) + NoLegend()
plot_newclusters_default_Myonuclei | plot_clusters_condition_Myonuclei
Myonuclei <- subset(Myonuclei, idents = c('Myonuclei_1', 'Myonuclei_2', 'Myonuclei_3'))
Cluster_colors <-c('coral2', 'burlywood', 'darkseagreen')
plot_newclusters_default_Myonuclei <- DimPlot(object = Myonuclei, label = FALSE, pt.size = 1.2, cols = Cluster_colors) + NoLegend()
plot_newclusters_condition_Myonuclei <- DimPlot(object = Myonuclei, label = FALSE, pt.size = 1.2, , group.by = 'sample_ID', cols = c('GFP' = 'grey', 'TEVp' = 'deepskyblue2')) 
plot_newclusters_default_Myonuclei | plot_newclusters_condition_Myonuclei

saveRDS(Myonuclei, file='/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_Myonuclei.rds')
DefaultAssay(Myonuclei) <- 'RNA'
Myonuclei_markers_Cluster1 <-FindMarkers(Myonuclei, ident.1 = 'Myonuclei_1', ident.2 = c('Myonuclei_2', 'Myonuclei_3'),
                                              only.pos = FALSE)
Myonuclei_markers_Cluster2 <-FindMarkers(Myonuclei, ident.1 = 'Myonuclei_2', ident.2 = c('Myonuclei_1', 'Myonuclei_3'),
                                              only.pos = FALSE)
Myonuclei_markers_Cluster3 <-FindMarkers(Myonuclei, ident.1 = 'Myonuclei_3', ident.2 = c('Myonuclei_1', 'Myonuclei_2'),
                                              only.pos = FALSE)
write.csv(Myonuclei_markers_Cluster1, "./ATAC_results/excels/Myonuclei_markers_cluster1.csv", row.names = TRUE)
write.csv(Myonuclei_markers_Cluster2, "./ATAC_results/excels/Myonuclei_markers_cluster2.csv", row.names = TRUE)
write.csv(Myonuclei_markers_Cluster3, "./ATAC_results/excels/Myonuclei_markers_cluster3.csv", row.names = TRUE)
Myonuclei_markers <- c(
  'Myh4',
  'Mylk4',
  'Pfkfb3',
  'Mybpc2',
  'Tpm1',
  'Pvalb',
  'Gpd2',
  'Ampd1',
  'Trdn',
  'Phkb',
  'Actn2',
  'Myh1',
  'Ptpn3',
  'Agbl1',
  'Fgf1',
  'Sorbs1',
  'Csrp3',
  'Fhl1',
  'Osbpl3',
  'Lmod2',
  'Crlf2',
  'Sln',
  'Runx1',
  'Ttc9',
  'Lipc',
  'Ptprj',
  'Ltbp1',
  'Cd82',
  'Musk',
  'Ncam1')

Dot_myonuclei_clusters <- DotPlot(
  object = Myonuclei,
  features = rev(Myonuclei_markers),
  dot.scale=10, 
  cols="RdBu"
) + coord_flip()
Dot_myonuclei_clusters


MuSCs_Myonuclei <- merge(Myonuclei, y = MuSCs)
MuSCs_Myonuclei <- RunTFIDF(MuSCs_Myonuclei)
MuSCs_Myonuclei <- FindTopFeatures(MuSCs_Myonuclei, min.cutoff = 'q0')
MuSCs_Myonuclei <- RunSVD(MuSCs_Myonuclei)
MuSCs_Myonuclei <- RunUMAP(MuSCs_Myonuclei, reduction = "lsi", dims = 2:30)
MuSCs_Myonuclei <- FindNeighbors(object = MuSCs_Myonuclei, reduction = 'lsi', dims = 2:30)
MuSCs_Myonuclei <- FindClusters(object = MuSCs_Myonuclei, verbose = FALSE, algorithm = 3, resolution = 0.95)
plot_clusters_default_MuSCs_Myonuclei<- DimPlot(object = MuSCs_Myonuclei, label = TRUE, pt.size = 1.2) + NoLegend()
plot_clusters_condition_MuSCs_Myonuclei <- DimPlot(object = MuSCs_Myonuclei, label = FALSE, pt.size = 1.2, group.by = 'sample_ID', cols = c('GFP' = 'grey', 'TEVp' = 'deepskyblue2')) + ggtitle('MuSCs_Myonuclei') + theme(legend.position = c(0.6, 0.1))
plot_clusters_default_MuSCs_Myonuclei | plot_clusters_condition_MuSCs_Myonuclei
DefaultAssay(MuSCs_Myonuclei) <- 'RNA'
new.cluster.default_MuSCs_Myonuclei <- c("Myonuclei_II", 
                                   "Myonuclei_II",
                                   "Myonuclei_II",
                                   "Myonuclei_II",
                                   "Myonuclei_TEVp",
                                   "Myonuclei_TEVp",
                                   "Myonuclei_II",
                                   "MuSCs"
)
names(new.cluster.default_MuSCs_Myonuclei) <- levels(MuSCs_Myonuclei)
MuSCs_Myonuclei <- RenameIdents(MuSCs_Myonuclei, new.cluster.default_MuSCs_Myonuclei)
plot_newclusters_default_MuSCs_Myonuclei <- DimPlot(object = MuSCs_Myonuclei, label = TRUE, pt.size = 1.2) + NoLegend()
plot_newclusters_default_MuSCs_Myonuclei | plot_clusters_condition_MuSCs_Myonuclei
MuSCs_Myonuclei <- subset(MuSCs_Myonuclei, idents = c('Myonuclei_IIb','Myonuclei_IIx','Myonuclei_TEVp', 'MuSCs'))
Cluster_colors <-c('coral2', 'burlywood', 'darkseagreen' ,'darkturquoise')
plot_newclusters_default_MuSCs_Myonuclei <- DimPlot(object = MuSCs_Myonuclei, label = FALSE, pt.size = 1.2, cols = Cluster_colors) + NoLegend()
plot_newclusters_condition_MuSCs_Myonuclei<- DimPlot(object = MuSCs_Myonuclei, label = FALSE, pt.size = 1.2, , group.by = 'sample_ID', cols = c('TEVp' = 'deepskyblue2', 'GFP' = 'grey')) 
plot_newclusters_default_MuSCs_Myonuclei | plot_newclusters_condition_MuSCs_Myonuclei
saveRDS(MuSCs_Myonuclei, file='/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_MuSCs_Myonuclei.rds')
