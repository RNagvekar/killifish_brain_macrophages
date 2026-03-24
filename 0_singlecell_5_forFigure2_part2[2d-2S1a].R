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

# Read Helena's immune object (for microglia vs. BAM comparison)
HBimmune <- qread("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/Immune.qs")

# Cell type labels
DimPlot(HBimmune, reduction = "umap", group.by = "Broad_Cell_Annotation", label = FALSE) +
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
FeaturePlot(HBimmune, "Dab2", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(HBimmune, features = c("Apoe")) +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none")

# Subset to ZT12 (beginning of active phase - matched to killifish circadian timepoint)
HBimmune_ZT12 <- subset(HBimmune, time == "ZT12")

# Define macrophage types (BAM = PBM + MDM + IFN-M)
HBimmune_ZT12$macrophagetype <- case_when(
  HBimmune_ZT12$Broad_Cell_Annotation %in% c("Microglia") ~ "Microglia",
  HBimmune_ZT12$Broad_Cell_Annotation %in% c("PBM", "MDM", "IFN-M") ~ "BAM",
  TRUE ~ "Non-macrophage")

# Markers distinguishing microglia from BAMs
DEG_Microglia_BAM_HBimmune_ZT12 <- FindMarkers(
  object = HBimmune_ZT12,
  ident.1 = "Microglia",
  ident.2 = "BAM",
  group.by = "macrophagetype",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.5
)
DEG_Microglia_BAM_HBimmune_ZT12_df <- as.data.frame(DEG_Microglia_BAM_HBimmune_ZT12)
DEG_Microglia_BAM_HBimmune_ZT12_df$gene <- rownames(DEG_Microglia_BAM_HBimmune_ZT12_df)

# Top 100 enriched in each
Microglia_top100markers_vsBAM_df <- DEG_Microglia_BAM_HBimmune_ZT12_df %>% dplyr::filter(avg_log2FC > 0.5) %>% dplyr::filter(pct.1 > 0.1) %>% dplyr::filter(p_val_adj < 0.05)
Microglia_top100markers <- head(Microglia_top100markers_vsBAM_df$gene, 100)

BAM_top100markers_vsMicroglia_df <- DEG_Microglia_BAM_HBimmune_ZT12_df %>% dplyr::filter(avg_log2FC < -0.5) %>% dplyr::filter(pct.2 > 0.1) %>% dplyr::filter(p_val_adj < 0.05)
BAM_top100markers <- head(BAM_top100markers_vsMicroglia_df$gene, 100)

top200markers <- union(Microglia_top100markers, BAM_top100markers)
top200markers_killifish <- c()
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Mouse[i] %in% top200markers) {
    top200markers_killifish <- append(top200markers_killifish, gene_names_df$`N. furzeri (NCBI)`[i])}
}

# Subset HBimmune_ZT12 to only Microglia and BAMs
HBmacrophage_ZT12 <- subset(HBimmune_ZT12, macrophagetype != "Non-macrophage")

# PCA based on top200markers
HBmacrophage_ZT12_topmarkers <- subset(HBmacrophage_ZT12, features = top200markers)
HBmacrophage_ZT12_topmarkers <- NormalizeData(HBmacrophage_ZT12_topmarkers, normalization.method = "LogNormalize", scale.factor = 10000)
HBmacrophage_ZT12_topmarkers <- FindVariableFeatures(HBmacrophage_ZT12_topmarkers, selection.method = "vst", nfeatures = 2000)
topmarkers <- rownames(HBmacrophage_ZT12_topmarkers)
HBmacrophage_ZT12_topmarkers <- ScaleData(HBmacrophage_ZT12_topmarkers, features = topmarkers)
HBmacrophage_ZT12_topmarkers <- RunPCA(HBmacrophage_ZT12_topmarkers, features = VariableFeatures(object = HBmacrophage_ZT12_topmarkers))
custom_colors <- c("Microglia" = "#5bb450", "BAM" = "#0076b6")

# Fig 2d
# PCA (mouse markers) for Helena ZT12 macrophages
# Save at 500W x 500H
DimPlot(HBmacrophage_ZT12_topmarkers, reduction = "pca", group.by = "macrophagetype", cols = custom_colors) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Get the standard deviation of each PC
stdev_HBmacrophage_ZT12_topmarkers <- Stdev(object = HBmacrophage_ZT12_topmarkers, reduction = "pca")

# Calculate variance explained (eigenvalues are squared sdev)
var_explained_HBmacrophage_ZT12_topmarkers <- (stdev_HBmacrophage_ZT12_topmarkers^2)/sum(stdev_HBmacrophage_ZT12_topmarkers^2)

# Re-process HBmacrophage_ZT12 for feature plots (standard processing)
HBmacrophage_ZT12 <- NormalizeData(HBmacrophage_ZT12, normalization.method = "LogNormalize", scale.factor = 10000)
HBmacrophage_ZT12 <- FindVariableFeatures(HBmacrophage_ZT12, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HBmacrophage_ZT12)
HBmacrophage_ZT12 <- ScaleData(HBmacrophage_ZT12, features = all.genes)
HBmacrophage_ZT12 <- RunPCA(HBmacrophage_ZT12, features = VariableFeatures(object = HBmacrophage_ZT12))
HBmacrophage_ZT12 <- RunUMAP(HBmacrophage_ZT12, dims = 1:30)

# Extract UMAP coordinates from the Seurat object
HBmacrophage_ZT12_umap_coords <- Embeddings(HBmacrophage_ZT12, "umap")
HBmacrophage_ZT12_xlims <- range(HBmacrophage_ZT12_umap_coords[,1])
HBmacrophage_ZT12_ylims <- range(HBmacrophage_ZT12_umap_coords[,2])

# UMAP
DimPlot(HBmacrophage_ZT12, group.by = "Broad_Cell_Annotation")

# Fig 2e
# Tmem119 feature plot
# Save at 500W x 500H
FeaturePlot(HBmacrophage_ZT12, "Tmem119", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2e
# Tmem119 feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Tmem119", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2e
# Mrc1 feature plot
# Save at 500W x 500H
FeaturePlot(HBmacrophage_ZT12, "Mrc1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2e
# Mrc1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Mrc1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# Lyve1 feature plot
# Save at 500W x 500H
FeaturePlot(HBmacrophage_ZT12, "Lyve1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# Lyve1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Lyve1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# Cd74 feature plot
# Save at 500W x 500H
FeaturePlot(HBmacrophage_ZT12, "Cd74", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# Cd74 feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Cd74", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# Csf1r feature plot
# Save at 500W x 500H
FeaturePlot(HBmacrophage_ZT12, "Csf1r", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# Csf1r feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Csf1r", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# Apoe feature plot
# Save at 500W x 500H
FeaturePlot(HBmacrophage_ZT12, "Apoe", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# Apoe feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Apoe", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Load YSR object
YSR_20251202 <- readRDS("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/YSR_20251202.rds")

# Extract UMAP coordinates from the Seurat object
YSR_umap_coords <- Embeddings(YSR_20251202, "umap")
YSR_xlims <- range(YSR_umap_coords[,1])
YSR_ylims <- range(YSR_umap_coords[,2])

# PCA based on top200markers_killifish
YSR_20251202_topmarkers <- subset(YSR_20251202, features = top200markers_killifish)
YSR_20251202_topmarkers <- NormalizeData(YSR_20251202_topmarkers, normalization.method = "LogNormalize", scale.factor = 10000)
YSR_20251202_topmarkers <- FindVariableFeatures(YSR_20251202_topmarkers, selection.method = "vst", nfeatures = 2000)
topmarkers <- rownames(YSR_20251202_topmarkers)
YSR_20251202_topmarkers <- ScaleData(YSR_20251202_topmarkers, features = topmarkers)
YSR_20251202_topmarkers <- RunPCA(YSR_20251202_topmarkers, features = VariableFeatures(object = YSR_20251202_topmarkers))

# Fig 2d
# PCA (mouse markers) for YSR
# Save at 500W x 500H
DimPlot(YSR_20251202_topmarkers, reduction = "pca", group.by = "orig.ident", cols = "red") +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Get the standard deviation of each PC
stdev_YSR_20251202_topmarkers <- Stdev(object = YSR_20251202_topmarkers, reduction = "pca")

# Calculate variance explained (eigenvalues are squared sdev)
var_explained_YSR_20251202_topmarkers <- (stdev_YSR_20251202_topmarkers^2)/sum(stdev_YSR_20251202_topmarkers^2)

# Fig 2e
# tmem119 feature plot
# Save at 500W x 500H
FeaturePlot(YSR_20251202, "tmem119", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2e
# tmem119 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "tmem119", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2e
# mrc1 feature plot
# Save at 500W x 500H
FeaturePlot(YSR_20251202, "LOC107383768", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2e
# mrc1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "LOC107383768", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# lyve1 feature plot
# Save at 500W x 500H
FeaturePlot(YSR_20251202, "lyve1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# lyve1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "lyve1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# cd74 feature plot
# Save at 500W x 500H
FeaturePlot(YSR_20251202, "cd74", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# cd74 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "cd74", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# csf1r feature plot
# Save at 500W x 500H
FeaturePlot(YSR_20251202, "csf1r", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# csf1r feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "csf1r", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# apoeb feature plot
# Save at 500W x 500H
FeaturePlot(YSR_20251202, "LOC107379395", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# apoeb feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "LOC107379395", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

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

# Extract UMAP coordinates from the Seurat object
SRR24058857_myeloid_umap_coords <- Embeddings(SRR24058857_myeloid_20260315, "umap")
SRR24058857_myeloid_xlims <- range(SRR24058857_myeloid_umap_coords[,1])
SRR24058857_myeloid_ylims <- range(SRR24058857_myeloid_umap_coords[,2])

# Feature plot
FeaturePlot(SRR24058857_myeloid_20260315, "cldn5", cols = c("yellow", "purple"))

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

# PCA based on top200markers_killifish
SRR24058857_myeloid_20260315_topmarkers <- subset(SRR24058857_myeloid_20260315, features = top200markers_killifish)
SRR24058857_myeloid_20260315_topmarkers <- NormalizeData(SRR24058857_myeloid_20260315_topmarkers, normalization.method = "LogNormalize", scale.factor = 10000)
SRR24058857_myeloid_20260315_topmarkers <- FindVariableFeatures(SRR24058857_myeloid_20260315_topmarkers, selection.method = "vst", nfeatures = 2000)
topmarkers <- rownames(SRR24058857_myeloid_20260315_topmarkers)
SRR24058857_myeloid_20260315_topmarkers <- ScaleData(SRR24058857_myeloid_20260315_topmarkers, features = topmarkers)
SRR24058857_myeloid_20260315_topmarkers <- RunPCA(SRR24058857_myeloid_20260315_topmarkers, features = VariableFeatures(object = SRR24058857_myeloid_20260315_topmarkers))

# Fig 2d
# PCA (mouse markers) for SRR24058857_myeloid
# Save at 500W x 500H
DimPlot(SRR24058857_myeloid_20260315_topmarkers, reduction = "pca", group.by = "orig.ident", cols = "black", pt.size = 2) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Get the standard deviation of each PC
stdev_SRR24058857_myeloid_20260315_topmarkers <- Stdev(object = SRR24058857_myeloid_20260315_topmarkers, reduction = "pca")

# Calculate variance explained (eigenvalues are squared sdev)
var_explained_SRR24058857_myeloid_20260315_topmarkers <- (stdev_SRR24058857_myeloid_20260315_topmarkers^2)/sum(stdev_SRR24058857_myeloid_20260315_topmarkers^2)

# Fig 2e
# tmem119 feature plot
# Save at 500W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "tmem119", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2e
# tmem119 feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "tmem119", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2e
# mrc1 feature plot
# Save at 500W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "LOC107383768", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2e
# mrc1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "LOC107383768", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# lyve1 feature plot
# Save at 500W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "lyve1", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# lyve1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "lyve1", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# cd74 feature plot
# Save at 500W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "cd74", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# cd74 feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "cd74", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# csf1r feature plot
# Save at 500W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "csf1r", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# csf1r feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "csf1r", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S1a
# apoeb feature plot
# Save at 500W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "LOC107379395", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S1a
# apoeb feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "LOC107379395", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Save Seurat objects
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(HBimmune_ZT12, file = "HBimmune_ZT12.rds")
saveRDS(HBmacrophage_ZT12, file = "HBmacrophage_ZT12.rds")
saveRDS(HBmacrophage_ZT12_topmarkers, file = "HBmacrophage_ZT12_topmarkers.rds")
saveRDS(YSR_20251202_topmarkers, file = "HBmacrophage_ZT12_topmarkers.rds")
saveRDS(SRR24058857_myeloid_20260315, file = "SRR24058857_myeloid_20260315.rds")
saveRDS(SRR24058857_myeloid_20260315_topmarkers, file = "SRR24058857_myeloid_20260315_topmarkers.rds")

# Reload Seurat objects
HBmacrophage_ZT12 <- readRDS(file = "HBmacrophage_ZT12.rds")
SRR24058857_myeloid_20260315 <- readRDS(file = "SRR24058857_myeloid_20260315.rds")