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
FeaturePlot(HBimmune, "Ybx1", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(HBimmune, features = c("Apoe")) +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none")

# Subset to ZT12 (beginning of active phase - matched to killifish circadian timepoint)
HBimmune_ZT12 <- subset(HBimmune, time == "ZT12")

# Subset HBimmune_ZT12 to only macrophages
HBmacrophage_ZT12 <- subset(HBimmune_ZT12, Broad_Cell_Annotation %in% c("Microglia", "PBM", "MDM", "IFN-M"))

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

# Fig 2d
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

# Fig 2d
# Tmem119 feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Tmem119", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2d
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

# Fig 2d
# Cd74 feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Cd74", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
# Cldn5 feature plot
# Save at 500W x 500H
FeaturePlot(HBmacrophage_ZT12, "Cldn5", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S3a
# Cldn5 feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Cldn5", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
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

# Fig 2S3a
# Mrc1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Mrc1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
# F13a1 feature plot
# Save at 500W x 500H
FeaturePlot(HBmacrophage_ZT12, "F13a1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S3a
# F13a1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "F13a1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
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

# Fig 2S3a
# Csf1r feature plot with legend
# Save at 600W x 500H
FeaturePlot(HBmacrophage_ZT12, "Csf1r", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = HBmacrophage_ZT12_xlims, ylim = HBmacrophage_ZT12_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
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

# Fig 2S3a
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

# Order cells
cells_order_YSR <- colnames(YSR_20251202)

# Extract UMAP coordinates from the Seurat object
YSR_umap_coords <- Embeddings(YSR_20251202, "umap")
YSR_xlims <- range(YSR_umap_coords[,1])
YSR_ylims <- range(YSR_umap_coords[,2])

# Fig 2d
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

# Fig 2d
# tmem119 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "tmem119", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2d
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

# Fig 2d
# cd74 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "cd74", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
# cldn5 feature plot
# Save at 500W x 500H
FeaturePlot(YSR_20251202, "cldn5", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S3a
# cldn5 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "cldn5", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
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

# Fig 2S3a
# mrc1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "LOC107383768", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
# f13a1 feature plot
# Save at 500W x 500H
FeaturePlot(YSR_20251202, "f13a1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S3a
# f13a1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "f13a1", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
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

# Fig 2S3a
# csf1r feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "csf1r", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
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

# Fig 2S3a
# apoeb feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_20251202, "LOC107379395", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_xlims, ylim = YSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Load Ayana young myeloid object (SRR24058857_myeloid)
SRR24058857_myeloid_20260315 <- readRDS("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/SRR24058857_myeloid_20260315.rds")

# Order cells
cells_order_SRR24058857_myeloid <- colnames(SRR24058857_myeloid_20260315)

# Extract UMAP coordinates from the Seurat object
SRR24058857_myeloid_umap_coords <- Embeddings(SRR24058857_myeloid_20260315, "umap")
SRR24058857_myeloid_xlims <- range(SRR24058857_myeloid_umap_coords[,1])
SRR24058857_myeloid_ylims <- range(SRR24058857_myeloid_umap_coords[,2])

# Fig 2d
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

# Fig 2d
# tmem119 feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "tmem119", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2d
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

# Fig 2d
# cd74 feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "cd74", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
# cldn5 feature plot
# Save at 500W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "cldn5", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S3a
# cldn5 feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "cldn5", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
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

# Fig 2S3a
# mrc1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "LOC107383768", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
# f13a1 feature plot
# Save at 500W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "f13a1", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 2S3a
# f13a1 feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "f13a1", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
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

# Fig 2S3a
# csf1r feature plot with legend
# Save at 600W x 500H
FeaturePlot(SRR24058857_myeloid_20260315, "csf1r", cols = c("yellow", "purple"), pt.size = 2) +
  coord_cartesian(xlim = SRR24058857_myeloid_xlims, ylim = SRR24058857_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 2S3a
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

# Fig 2S3a
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

# Reload Seurat objects
HBmacrophage_ZT12 <- readRDS(file = "HBmacrophage_ZT12.rds")