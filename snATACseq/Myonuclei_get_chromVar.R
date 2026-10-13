library(Seurat) 
library(Signac) 
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(chromVAR)
library(here)
library(patchwork)
set.seed(1234)

ATAC_Myonuclei_ <- "/Volumes/SC/LAB_JAC/LAB/ManuG/12_Omics/2024.03_snMultiome/ManuelGavilan_SC_Multiome_GeX_10XGenomics/ATAC_results/ATAC_Myonuclei.rds"
ATAC_Myonuclei <- readRDS(ATAC_Myonuclei_)

features <- rownames(ATAC_Myonuclei)
granges_features <- StringToGRanges(features)
feature_seqlevels <- seqlevels(granges_features)
genome_seqlevels <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
missing_seqlevels <- setdiff(feature_seqlevels, genome_seqlevels)
granges_filtered <- dropSeqlevels(granges_features, missing_seqlevels, pruning.mode="coarse")
filtered_features <- paste0(seqnames(granges_filtered), "-", start(granges_filtered), "-", end(granges_filtered))

ATAC_Myonuclei <- subset(ATAC_Myonuclei, features = filtered_features)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(ATAC_Myonuclei)),
  pwm = pfm,
  genome = 'BSgenome.Mmusculus.UCSC.mm10',
)

motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

ATAC_Myonuclei <- AddMotifs(
  object = ATAC_Myonuclei,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

ATAC_Myonuclei <- RunChromVAR(
  object = ATAC_Myonuclei,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  motif.matrix = motif.matrix
)

saveRDS(ATAC_Myonuclei, file = "ATAC_Myonuclei_chromVar.rds")
saveRDS(motif, file = "ATAC_Myonuclei_motif.rds")
