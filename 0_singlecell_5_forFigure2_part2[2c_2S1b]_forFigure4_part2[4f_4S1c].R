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
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.5
)
celltypemarkers_HB_microglia_BAM_MDM_ZT12_df <- as.data.frame(celltypemarkers_HB_microglia_BAM_MDM_ZT12)
celltypemarkers_HB_microglia_BAM_MDM_ZT12_df$gene <- rownames(celltypemarkers_HB_microglia_BAM_MDM_ZT12_df)

# Top markers enriched in each
Microglia_markers_df <- celltypemarkers_HB_microglia_BAM_MDM_ZT12_df %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::filter(cluster == "Microglia")
Microglia_markers <- Microglia_markers_df$gene
Microglia_markers_killifish <- c()
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Mouse[i] %in% Microglia_markers) {
    Microglia_markers_killifish <- append(Microglia_markers_killifish, gene_names_df$`N. furzeri (NCBI)`[i])}
}

PBM_markers_df <- celltypemarkers_HB_microglia_BAM_MDM_ZT12_df %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::filter(cluster == "PBM")
PBM_markers <- PBM_markers_df$gene
PBM_markers_killifish <- c()
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Mouse[i] %in% PBM_markers) {
    PBM_markers_killifish <- append(PBM_markers_killifish, gene_names_df$`N. furzeri (NCBI)`[i])}
}

MDM_markers_df <- celltypemarkers_HB_microglia_BAM_MDM_ZT12_df %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::filter(cluster == "MDM")
MDM_markers <- MDM_markers_df$gene
MDM_markers_killifish <- c()
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Mouse[i] %in% MDM_markers) {
    MDM_markers_killifish <- append(MDM_markers_killifish, gene_names_df$`N. furzeri (NCBI)`[i])}
}

markers <- Reduce(union, list(Microglia_markers, PBM_markers, MDM_markers))
markers_killifish <- Reduce(union, list(Microglia_markers_killifish, PBM_markers_killifish, MDM_markers_killifish))

# PCA based on markers
HB_microglia_BAM_MDM_ZT12_markers <- subset(HB_microglia_BAM_MDM_ZT12, features = markers)
HB_microglia_BAM_MDM_ZT12_markers <- NormalizeData(HB_microglia_BAM_MDM_ZT12_markers, normalization.method = "LogNormalize", scale.factor = 10000)
HB_microglia_BAM_MDM_ZT12_markers <- FindVariableFeatures(HB_microglia_BAM_MDM_ZT12_markers, selection.method = "vst", nfeatures = 2000)
markers <- rownames(HB_microglia_BAM_MDM_ZT12_markers)
HB_microglia_BAM_MDM_ZT12_markers <- ScaleData(HB_microglia_BAM_MDM_ZT12_markers, features = markers)
HB_microglia_BAM_MDM_ZT12_markers <- RunPCA(HB_microglia_BAM_MDM_ZT12_markers, features = VariableFeatures(object = HB_microglia_BAM_MDM_ZT12_markers))
custom_colors <- c("Microglia" = "#5bb450", "PBM" = "skyblue","MDM" = "blue")

# Fig 2c
# PCA (mouse markers) for Helena ZT12 microglia, BAM, MDM (PC1 vs. PC2)
# Save at 500W x 500H
DimPlot(HB_microglia_BAM_MDM_ZT12_markers, reduction = "pca", group.by = "Broad_Cell_Annotation", cols = custom_colors) +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1b
# PCA (mouse markers) for Helena ZT12 microglia, BAM, MDM (PC1 vs. PC3)
# Save at 500W x 500H
DimPlot(HB_microglia_BAM_MDM_ZT12_markers, reduction = "pca", dims = c(1, 3), group.by = "Broad_Cell_Annotation", cols = custom_colors) +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Get the standard deviation of each PC
stdev_HB_microglia_BAM_MDM_ZT12_markers <- Stdev(object = HB_microglia_BAM_MDM_ZT12_markers, reduction = "pca")

# Calculate variance explained (eigenvalues are squared sdev)
var_explained_HB_microglia_BAM_MDM_ZT12_markers <- (stdev_HB_microglia_BAM_MDM_ZT12_markers^2)/sum(stdev_HB_microglia_BAM_MDM_ZT12_markers^2)

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

# PCA based on markers_killifish
YSR_20251202_markers <- subset(YSR_20251202, features = markers_killifish)
YSR_20251202_markers <- NormalizeData(YSR_20251202_markers, normalization.method = "LogNormalize", scale.factor = 10000)
YSR_20251202_markers <- FindVariableFeatures(YSR_20251202_markers, selection.method = "vst", nfeatures = 2000)
markers <- rownames(YSR_20251202_markers)
YSR_20251202_markers <- ScaleData(YSR_20251202_markers, features = markers)
YSR_20251202_markers <- RunPCA(YSR_20251202_markers, features = VariableFeatures(object = YSR_20251202_markers))

# Fig 2c
# PCA (mouse markers) for YSR (PC1 vs. PC2)
# Save at 500W x 500H
DimPlot(YSR_20251202_markers, reduction = "pca", group.by = "orig.ident", cols = "red") +
  ggtitle(NULL) +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1b
# PCA (mouse markers) for YSR (PC1 vs. PC3)
# Save at 500W x 500H
DimPlot(YSR_20251202_markers, reduction = "pca", dims = c(1, 3), group.by = "orig.ident", cols = "red") +
  ggtitle(NULL) +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Get the standard deviation of each PC
stdev_YSR_20251202_markers <- Stdev(object = YSR_20251202_markers, reduction = "pca")

# Calculate variance explained (eigenvalues are squared sdev)
var_explained_YSR_20251202_markers <- (stdev_YSR_20251202_markers^2)/sum(stdev_YSR_20251202_markers^2)

# Heatmap of microglia markers with legend
DoHeatmap(YSR_20251202,
          features = Microglia_markers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_YSR,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0, 6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Heatmap of PBM markers with legend
DoHeatmap(YSR_20251202,
          features = PBM_markers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_YSR,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0, 6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Heatmap of MDM markers with legend
DoHeatmap(YSR_20251202,
          features = MDM_markers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_YSR,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0, 6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Load Ayana young object (SRR24058857)
SRR24058857_20260315 <- readRDS("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/SRR24058857_20260315.rds")

# UMAP
SRR24058857_20260315 <- FindClusters(SRR24058857_20260315, resolution = 0.1)
DimPlot(SRR24058857_20260315, reduction = "umap", label = TRUE)

# Feature plot
FeaturePlot(SRR24058857_20260315, "LOC107379395") # apoeb
FeaturePlot(SRR24058857_20260315, "LOC107378674") # iba1
FeaturePlot(SRR24058857_20260315, "lcp1")
FeaturePlot(SRR24058857_20260315, "csf1r")

# Violin plot
VlnPlot(SRR24058857_20260315, features = c("LOC107379395")) # apoeb
VlnPlot(SRR24058857_20260315, features = c("LOC107378674")) # iba1
VlnPlot(SRR24058857_20260315, features = c("lcp1"))
VlnPlot(SRR24058857_20260315, features = c("csf1r"))

# Subset to myeloid (cluster 4 at resolution = 0.1)
SRR24058857_myeloid_20260315 <- subset(SRR24058857_20260315, seurat_clusters == "4")

# Normalization
SRR24058857_myeloid_20260315 <- NormalizeData(SRR24058857_myeloid_20260315, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
SRR24058857_myeloid_20260315 <- FindVariableFeatures(SRR24058857_myeloid_20260315, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(SRR24058857_myeloid_20260315)
SRR24058857_myeloid_20260315 <- ScaleData(SRR24058857_myeloid_20260315, features = all.genes)

# PCA
SRR24058857_myeloid_20260315 <- RunPCA(SRR24058857_myeloid_20260315, features = VariableFeatures(object = SRR24058857_myeloid_20260315))

# UMAP
SRR24058857_myeloid_20260315 <- RunUMAP(SRR24058857_myeloid_20260315, dims = 1:30)
SRR24058857_myeloid_20260315 <- FindNeighbors(SRR24058857_myeloid_20260315, dims = 1:30)
SRR24058857_myeloid_20260315 <- FindClusters(SRR24058857_myeloid_20260315, resolution = 0.5)
DimPlot(SRR24058857_myeloid_20260315, reduction = "umap", label = TRUE)

# Order cells
cells_order_SRR24058857_myeloid <- colnames(SRR24058857_myeloid_20260315)

# Extract UMAP coordinates from the Seurat object
SRR24058857_myeloid_umap_coords <- Embeddings(SRR24058857_myeloid_20260315, "umap")
SRR24058857_myeloid_xlims <- range(SRR24058857_myeloid_umap_coords[,1])
SRR24058857_myeloid_ylims <- range(SRR24058857_myeloid_umap_coords[,2])

# Feature plot
FeaturePlot(SRR24058857_myeloid_20260315, "f13a1", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(SRR24058857_myeloid_20260315, features = c("sall3"))

# Cluster markers
SRR24058857_myeloid_20260315_clustermarkers_df <- as.data.frame(FindAllMarkers(SRR24058857_myeloid_20260315)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(SRR24058857_myeloid_20260315_clustermarkers_df)) {
  killifish_NCBI_gene_name = SRR24058857_myeloid_20260315_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {SRR24058857_myeloid_20260315_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {SRR24058857_myeloid_20260315_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {SRR24058857_myeloid_20260315_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {SRR24058857_myeloid_20260315_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {SRR24058857_myeloid_20260315_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {SRR24058857_myeloid_20260315_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# PCA based on markers_killifish
SRR24058857_myeloid_20260315_markers <- subset(SRR24058857_myeloid_20260315, features = markers_killifish)
SRR24058857_myeloid_20260315_markers <- NormalizeData(SRR24058857_myeloid_20260315_markers, normalization.method = "LogNormalize", scale.factor = 10000)
SRR24058857_myeloid_20260315_markers <- FindVariableFeatures(SRR24058857_myeloid_20260315_markers, selection.method = "vst", nfeatures = 2000)
markers <- rownames(SRR24058857_myeloid_20260315_markers)
SRR24058857_myeloid_20260315_markers <- ScaleData(SRR24058857_myeloid_20260315_markers, features = markers)
SRR24058857_myeloid_20260315_markers <- RunPCA(SRR24058857_myeloid_20260315_markers, features = VariableFeatures(object = SRR24058857_myeloid_20260315_markers))

# Fig 2c
# PCA (mouse markers) for SRR24058857_myeloid (PC1 vs. PC2)
# Save at 500W x 500H
DimPlot(SRR24058857_myeloid_20260315_markers, reduction = "pca", group.by = "orig.ident", cols = "black", pt.size = 2) +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1b
# PCA (mouse markers) for SRR24058857_myeloid (PC1 vs. PC3)
# Save at 500W x 500H
DimPlot(SRR24058857_myeloid_20260315_markers, reduction = "pca", dims = c(1, 3), group.by = "orig.ident", cols = "black", pt.size = 2) +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Get the standard deviation of each PC
stdev_SRR24058857_myeloid_20260315_markers <- Stdev(object = SRR24058857_myeloid_20260315_markers, reduction = "pca")

# Calculate variance explained (eigenvalues are squared sdev)
var_explained_SRR24058857_myeloid_20260315_markers <- (stdev_SRR24058857_myeloid_20260315_markers^2)/sum(stdev_SRR24058857_myeloid_20260315_markers^2)

# Heatmap of microglia markers with legend
DoHeatmap(SRR24058857_myeloid_20260315,
          features = Microglia_markers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_SRR24058857_myeloid,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0,6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Heatmap of BAM markers with legend
DoHeatmap(SRR24058857_myeloid_20260315,
          features = PBM_markers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_SRR24058857_myeloid,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0,6)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# Heatmap of MDM markers with legend
DoHeatmap(SRR24058857_myeloid_20260315,
          features = MDM_markers_killifish,
          slot = "data",
          group.by = "orig.ident",
          cells = cells_order_SRR24058857_myeloid,
          label = FALSE,
          group.bar = FALSE) +
  scale_fill_gradientn(colors = c("yellow", "purple"), limits = c(0,5)) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

# N.B. For the creation of SR_WL_myeloid object please refer to script 0_singlecell_9_forFigure4_part3[4g-i_4S1d-e].R

# Load SR_WL_myeloid object
# N.B. These are all old cells
SR_WL_myeloid_20251202 <- readRDS("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/SR_WL_myeloid_20251202.rds")

# PCA based on markers_killifish
SR_WL_myeloid_20251202_markers <- subset(SR_WL_myeloid_20251202, features = markers_killifish)
SR_WL_myeloid_20251202_markers <- NormalizeData(SR_WL_myeloid_20251202_markers, normalization.method = "LogNormalize", scale.factor = 10000)
SR_WL_myeloid_20251202_markers <- FindVariableFeatures(SR_WL_myeloid_20251202_markers, selection.method = "vst", nfeatures = 2000)
markers <- rownames(SR_WL_myeloid_20251202_markers)
SR_WL_myeloid_20251202_markers <- ScaleData(SR_WL_myeloid_20251202_markers, features = markers)
SR_WL_myeloid_20251202_markers <- RunPCA(SR_WL_myeloid_20251202_markers, features = VariableFeatures(object = SR_WL_myeloid_20251202_markers))

# Fig 4f
# PCA (mouse markers) for SR_WL_myeloid (PC1 vs. PC2)
# Save at 500W x 500H
DimPlot(SR_WL_myeloid_20251202_markers, reduction = "pca", group.by = "orig.ident", cols = c("purple", "darkgrey")) +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 4S1c
# PCA (mouse markers) for SR_WL_myeloid (PC1 vs. PC3)
# Save at 500W x 500H
DimPlot(SR_WL_myeloid_20251202_markers, reduction = "pca", dims = c(1, 3), group.by = "orig.ident", cols = c("purple", "darkgrey")) +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Get the standard deviation of each PC
stdev_SR_WL_myeloid_20251202_markers <- Stdev(object = SR_WL_myeloid_20251202_markers, reduction = "pca")

# Calculate variance explained (eigenvalues are squared sdev)
var_explained_SR_WL_myeloid_20251202_markers <- (stdev_SR_WL_myeloid_20251202_markers^2)/sum(stdev_SR_WL_myeloid_20251202_markers^2)

# Subset to SR_WL_myeloid to wildtype (WL_myeloid)
# N.B. The cells selected would be slightly different if subset only from WL, subsetting from a combined object with SR cells helps exclude potential myeloid/non-myeloid doublets from WL_myeloid
WL_myeloid_20251202 <- subset(SR_WL_myeloid_20251202, orig.ident %in% c("WL"))

# PCA based on markers_killifish
WL_myeloid_20251202_markers <- subset(WL_myeloid_20251202, features = markers_killifish)
WL_myeloid_20251202_markers <- NormalizeData(WL_myeloid_20251202_markers, normalization.method = "LogNormalize", scale.factor = 10000)
WL_myeloid_20251202_markers <- FindVariableFeatures(WL_myeloid_20251202_markers, selection.method = "vst", nfeatures = 2000)
markers <- rownames(WL_myeloid_20251202_markers)
WL_myeloid_20251202_markers <- ScaleData(WL_myeloid_20251202_markers, features = markers)
WL_myeloid_20251202_markers <- RunPCA(WL_myeloid_20251202_markers, features = VariableFeatures(object = WL_myeloid_20251202_markers))

# PCA (mouse markers) for WL_myeloid (PC1 vs. PC2)
# Save at 500W x 500H
DimPlot(WL_myeloid_20251202_markers, reduction = "pca", group.by = "orig.ident", cols = "grey") +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# PCA (mouse markers) for WL_myeloid (PC1 vs. PC3)
# Save at 500W x 500H
DimPlot(WL_myeloid_20251202_markers, reduction = "pca", dims = c(1, 3), group.by = "orig.ident", cols = "grey") +
  xlim(c(-70, 70)) +
  ylim(c(-70, 70)) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Get the standard deviation of each PC
stdev_WL_myeloid_20251202_markers <- Stdev(object = WL_myeloid_20251202_markers, reduction = "pca")

# Calculate variance explained (eigenvalues are squared sdev)
var_explained_WL_myeloid_20251202_markers <- (stdev_WL_myeloid_20251202_markers^2)/sum(stdev_WL_myeloid_20251202_markers^2)

# Save Seurat objects
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(HBimmune_ZT12, file = "HBimmune_ZT12.rds")
saveRDS(HB_microglia_BAM_MDM_ZT12, file = "HB_microglia_BAM_MDM_ZT12.rds")
saveRDS(HB_microglia_BAM_MDM_ZT12_markers, file = "HB_microglia_BAM_MDM_ZT12_markers.rds")
saveRDS(YSR_20251202_markers, file = "YSR_20251202_markers.rds")
saveRDS(SRR24058857_myeloid_20260315, file = "SRR24058857_myeloid_20260315.rds")
saveRDS(SRR24058857_myeloid_20260315_markers, file = "SRR24058857_myeloid_20260315_markers.rds")
saveRDS(WL_myeloid_20251202, file = "WL_myeloid_20251202.rds")
saveRDS(WL_myeloid_20251202_markers, file = "WL_myeloid_20251202_markers.rds")

# Reload Seurat objects
HB_microglia_BAM_MDM_ZT12 <- readRDS(file = "HB_microglia_BAM_MDM_ZT12.rds")
SRR24058857_myeloid_20260315 <- readRDS(file = "SRR24058857_myeloid_20260315.rds")
WL_myeloid_20251202 <- readRDS(file = "WL_myeloid_20251202.rds")