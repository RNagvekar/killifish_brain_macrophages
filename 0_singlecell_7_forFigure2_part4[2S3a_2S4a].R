rm(list = ls())
library(Seurat)
library(dplyr)
library(readxl)
library(ggplot2)
library(biomaRt)
library(enrichR)
library(presto)
library(tidyr)
library(scales)
library(qs)
library(DoubletFinder)

update.packages(ask = FALSE)

# Set working directory
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202")

# Read Helena's immune object
HBimmune <- qread("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/Immune.qs")

# Cell type labels
DimPlot(HBimmune, reduction = "umap", group.by = "Broad_Cell_Annotation", label = TRUE) +
  theme(legend.position = "none")

# Circadian time labels
DimPlot(HBimmune, reduction = "umap", group.by = "time", label = FALSE)

# Highlight specific cells
DimPlot(HBimmune,
        cells.highlight = WhichCells(HBimmune, expression = Broad_Cell_Annotation == "Microglia"),
        cols.highlight = "red",
        cols = "grey",
        pt.size = 0.1,
        sizes.highlight = 0.1
)

# Feature plot
FeaturePlot(HBimmune, "Apoe", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(HBimmune, features = c("Apoe")) +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none")

# Subset to ZT12 (beginning of active phase - matched to killifish circadian timepoint)
HBimmune_ZT12 <- subset(HBimmune, time == "ZT12")

# Subset to microglia, PBM [BAM], and MDM, ignore IFN-M
HB_microglia_BAM_MDM_ZT12 <- subset(HBimmune_ZT12, Broad_Cell_Annotation %in% c("Microglia", "PBM", "MDM"))

# Markers distinguishing microglia, PBM, MDM
celltypemarkers_HB_microglia_BAM_MDM_ZT12 <- FindAllMarkers(
  object = HB_microglia_BAM_MDM_ZT12,
  group.by = "Broad_Cell_Annotation",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  min.pct = 0.5,
  logfc.threshold = 3,
  only.pos = TRUE
)
celltypemarkers_HB_microglia_BAM_MDM_ZT12_df <- as.data.frame(celltypemarkers_HB_microglia_BAM_MDM_ZT12)
celltypemarkers_HB_microglia_BAM_MDM_ZT12_df$gene <- rownames(celltypemarkers_HB_microglia_BAM_MDM_ZT12_df)

# Top markers enriched in each
Microglia_strongmarkers_df <- celltypemarkers_HB_microglia_BAM_MDM_ZT12_df %>% dplyr::filter(cluster == "Microglia")
Microglia_strongmarkers <- Microglia_strongmarkers_df$gene
Microglia_strongmarkers_killifish <- c()
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Mouse[i] %in% Microglia_strongmarkers) {
    Microglia_strongmarkers_killifish <- append(Microglia_strongmarkers_killifish, gene_names_df$`N. furzeri (NCBI)`[i])}
}

PBM_strongmarkers_df <- celltypemarkers_HB_microglia_BAM_MDM_ZT12_df %>% dplyr::filter(cluster == "PBM")
PBM_strongmarkers <- PBM_strongmarkers_df$gene
PBM_strongmarkers_killifish <- c()
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Mouse[i] %in% PBM_strongmarkers) {
    PBM_strongmarkers_killifish <- append(PBM_strongmarkers_killifish, gene_names_df$`N. furzeri (NCBI)`[i])}
}

MDM_strongmarkers_df <- celltypemarkers_HB_microglia_BAM_MDM_ZT12_df %>% dplyr::filter(cluster == "MDM")
MDM_strongmarkers <- MDM_strongmarkers_df$gene
MDM_strongmarkers_killifish <- c()
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Mouse[i] %in% MDM_strongmarkers) {
    MDM_strongmarkers_killifish <- append(MDM_strongmarkers_killifish, gene_names_df$`N. furzeri (NCBI)`[i])}
}

strongmarkers <- Reduce(union, list(Microglia_strongmarkers, PBM_strongmarkers, MDM_strongmarkers))
strongmarkers_killifish <- Reduce(union, list(Microglia_strongmarkers_killifish, PBM_strongmarkers_killifish, MDM_strongmarkers_killifish))

# Load YSR object
YSR_20251202 <- readRDS("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/YSR_20251202.rds")

# Order cells
cells_order_YSR <- colnames(YSR_20251202)

# Extract UMAP coordinates from the Seurat object
YSR_umap_coords <- Embeddings(YSR_20251202, "umap")
YSR_xlims <- range(YSR_umap_coords[,1])
YSR_ylims <- range(YSR_umap_coords[,2])

# Feature plot
FeaturePlot(YSR_20251202, "LOC107379398", cols = c("yellow", "purple"))

# Fig 2S3a
# Heatmap of microglia markers with legend
# Save at 300W x 600H
DoHeatmap(YSR_20251202,
          features = Microglia_strongmarkers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_YSR,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0, 6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_text(size = 5, face = "italic"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Fig 2S3a
# Heatmap of PBM markers with legend
# Save at 300W x 500H
DoHeatmap(YSR_20251202,
          features = PBM_strongmarkers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_YSR,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0, 6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_text(size = 5, face = "italic"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Fig 2S3a
# Heatmap of MDM markers with legend
# Save at 300W x 200H
DoHeatmap(YSR_20251202,
          features = MDM_strongmarkers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_YSR,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0, 6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_text(size = 5, face = "italic"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Load Ayana young myeloid object (SRR24058857_myeloid)
SRR24058857_myeloid_20260315 <- readRDS("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/SRR24058857_myeloid_20260315.rds")

# Order cells
cells_order_SRR24058857_myeloid <- colnames(SRR24058857_myeloid_20260315)

# Extract UMAP coordinates from the Seurat object
SRR24058857_myeloid_umap_coords <- Embeddings(SRR24058857_myeloid_20260315, "umap")
SRR24058857_myeloid_xlims <- range(SRR24058857_myeloid_umap_coords[,1])
SRR24058857_myeloid_ylims <- range(SRR24058857_myeloid_umap_coords[,2])

# Feature plot
FeaturePlot(SRR24058857_myeloid_20260315, "f13a1", cols = c("yellow", "purple"))

# Fig 2S4a
# Heatmap of microglia markers with legend
# Save at 300W x 600H
DoHeatmap(SRR24058857_myeloid_20260315,
          features = Microglia_strongmarkers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_SRR24058857_myeloid,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0,6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_text(size = 5, face = "italic"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Fig 2S4a
# Heatmap of BAM markers with legend
# Save at 300W x 500H
DoHeatmap(SRR24058857_myeloid_20260315,
          features = PBM_strongmarkers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_SRR24058857_myeloid,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0,6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_text(size = 5, face = "italic"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Fig 2S4a
# Heatmap of MDM markers with legend
# Save at 300W x 200H
DoHeatmap(SRR24058857_myeloid_20260315,
          features = MDM_strongmarkers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_SRR24058857_myeloid,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0,6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_text(size = 5, face = "italic"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))