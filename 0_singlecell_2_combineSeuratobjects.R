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

# ----------------------------------------------------------------------------------------------------
# YSR_YSNR
# ----------------------------------------------------------------------------------------------------

# Read in Seurat objects
YSR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/YSR_20251202.rds")
YSNR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/YSNR_20251202.rds")

# Merge
YSR_YSNR_20251202 <- merge(YSR_20251202, y = YSNR_20251202)
table(YSR_YSNR_20251202@meta.data$orig.ident)

# Join layers
YSR_YSNR_20251202 <- JoinLayers(YSR_YSNR_20251202)

# Normalization
YSR_YSNR_20251202 <- NormalizeData(YSR_YSNR_20251202, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
YSR_YSNR_20251202 <- FindVariableFeatures(YSR_YSNR_20251202, selection.method = "vst", nfeatures = 2000)

# Scaling 
all.genes <- rownames(YSR_YSNR_20251202)
YSR_YSNR_20251202 <- ScaleData(YSR_YSNR_20251202, features = all.genes)

# PCA
YSR_YSNR_20251202 <- RunPCA(YSR_YSNR_20251202, features = VariableFeatures(object = YSR_YSNR_20251202))

# UMAP
YSR_YSNR_20251202 <- RunUMAP(YSR_YSNR_20251202, dims = 1:30)
YSR_YSNR_20251202 <- FindNeighbors(YSR_YSNR_20251202, dims = 1:30)
YSR_YSNR_20251202 <- FindClusters(YSR_YSNR_20251202, resolution = 2)
DimPlot(YSR_YSNR_20251202, reduction = "umap", label = TRUE)
DimPlot(YSR_YSNR_20251202, reduction = "umap", group.by = "celltype")

# Feature plot
FeaturePlot(YSR_YSNR_20251202, "LOC107386542", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(YSR_YSNR_20251202, features = c("sall1"))

# QC metrics
VlnPlot(YSR_YSNR_20251202, features = c("percent.mt")) 
VlnPlot(YSR_YSNR_20251202, features = c("nFeature_RNA"))
VlnPlot(YSR_YSNR_20251202, features = c("nCount_RNA"))
YSR_YSNR_20251202[["percent.ribo"]] <- PercentageFeatureSet(YSR_YSNR_20251202, pattern = "^rps|^rpl")
VlnPlot(YSR_YSNR_20251202, features = c("percent.ribo"))

# Cluster markers
YSR_YSNR_20251202_clustermarkers_df <- as.data.frame(FindAllMarkers(YSR_YSNR_20251202)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Cluster 32 (res = 2) markers
markers_cluster32 <- FindMarkers(
  object = YSR_YSNR_20251202,
  ident.1 = 32,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.5
)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(YSR_YSNR_20251202_clustermarkers_df)) {
  killifish_NCBI_gene_name = YSR_YSNR_20251202_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {YSR_YSNR_20251202_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {YSR_YSNR_20251202_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_YSNR_20251202_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {YSR_YSNR_20251202_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_YSNR_20251202_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {YSR_YSNR_20251202_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Save and reload RDS
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(YSR_YSNR_20251202, file = "YSR_YSNR_20251202.rds")

# ----------------------------------------------------------------------------------------------------
# YSR_OSR
# ----------------------------------------------------------------------------------------------------

# Read in Seurat objects
YSR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/YSR_20251202.rds")
OSR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/OSR_20251202.rds")

# Merge
YSR_OSR_20251202 <- merge(YSR_20251202, y = OSR_20251202)
table(YSR_OSR_20251202@meta.data$orig.ident)

# Join layers
YSR_OSR_20251202 <- JoinLayers(YSR_OSR_20251202)

# Normalization
YSR_OSR_20251202 <- NormalizeData(YSR_OSR_20251202, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
YSR_OSR_20251202 <- FindVariableFeatures(YSR_OSR_20251202, selection.method = "vst", nfeatures = 2000)

# Scaling 
all.genes <- rownames(YSR_OSR_20251202)
YSR_OSR_20251202 <- ScaleData(YSR_OSR_20251202, features = all.genes)

# PCA
YSR_OSR_20251202 <- RunPCA(YSR_OSR_20251202, features = VariableFeatures(object = YSR_OSR_20251202))

# UMAP
YSR_OSR_20251202 <- RunUMAP(YSR_OSR_20251202, dims = 1:30)
YSR_OSR_20251202 <- FindNeighbors(YSR_OSR_20251202, dims = 1:30)
YSR_OSR_20251202 <- FindClusters(YSR_OSR_20251202, resolution = 0.25)
DimPlot(YSR_OSR_20251202, reduction = "umap", label = TRUE)
DimPlot(YSR_OSR_20251202, reduction = "umap", group.by = "orig.ident")

# Feature plot
FeaturePlot(YSR_OSR_20251202, "marco", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(YSR_OSR_20251202, features = c("marco")) +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95)
VlnPlot(YSR_OSR_20251202, features = c("atp6v0a1"), group.by = "orig.ident") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95)
VlnPlot(YSR_OSR_20251202, features = c("sall3"), group.by = "seurat_clusters", split.by = "orig.ident")

# QC metrics
VlnPlot(YSR_OSR_20251202, features = c("percent.mt")) 
VlnPlot(YSR_OSR_20251202, features = c("nFeature_RNA"))
VlnPlot(YSR_OSR_20251202, features = c("nCount_RNA"))
YSR_OSR_20251202[["percent.ribo"]] <- PercentageFeatureSet(YSR_OSR_20251202, pattern = "^rps|^rpl")
VlnPlot(YSR_OSR_20251202, features = c("percent.ribo"))

# Cluster markers
YSR_OSR_20251202_clustermarkers_df <- as.data.frame(FindAllMarkers(YSR_OSR_20251202)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(YSR_OSR_20251202_clustermarkers_df)) {
  killifish_NCBI_gene_name = YSR_OSR_20251202_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {YSR_OSR_20251202_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {YSR_OSR_20251202_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_OSR_20251202_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {YSR_OSR_20251202_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_OSR_20251202_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {YSR_OSR_20251202_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Save and reload RDS
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(YSR_OSR_20251202, file = "YSR_OSR_20251202.rds")
YSR_OSR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/YSR_OSR_20251202.rds")

# ----------------------------------------------------------------------------------------------------
# SR_WL
# ----------------------------------------------------------------------------------------------------

# Read in Seurat objects
SR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/SR_20251202.rds")
WL_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/WL_20251202.rds")

# Merge
SR_WL_20251202 <- merge(SR_20251202, y = WL_20251202)
table(SR_WL_20251202@meta.data$orig.ident)

# Join layers
SR_WL_20251202 <- JoinLayers(SR_WL_20251202)

# Normalization
SR_WL_20251202 <- NormalizeData(SR_WL_20251202, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
SR_WL_20251202 <- FindVariableFeatures(SR_WL_20251202, selection.method = "vst", nfeatures = 2000)

# Scaling 
all.genes <- rownames(SR_WL_20251202)
SR_WL_20251202 <- ScaleData(SR_WL_20251202, features = all.genes)

# PCA
SR_WL_20251202 <- RunPCA(SR_WL_20251202, features = VariableFeatures(object = SR_WL_20251202))

# UMAP
SR_WL_20251202 <- RunUMAP(SR_WL_20251202, dims = 1:30)
SR_WL_20251202 <- FindNeighbors(SR_WL_20251202, dims = 1:30)
SR_WL_20251202 <- FindClusters(SR_WL_20251202, resolution = 0.5)
DimPlot(SR_WL_20251202, reduction = "umap", label = TRUE)
DimPlot(SR_WL_20251202, reduction = "umap", group.by = "orig.ident")

# Feature plot
FeaturePlot(SR_WL_20251202, "osr2", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(SR_WL_20251202, features = c("osr2"))

# QC metrics
VlnPlot(SR_WL_20251202, features = c("percent.mt")) 
VlnPlot(SR_WL_20251202, features = c("nFeature_RNA"))
VlnPlot(SR_WL_20251202, features = c("nCount_RNA"))
SR_WL_20251202[["percent.ribo"]] <- PercentageFeatureSet(YSR_OSR_20251202, pattern = "^rps|^rpl")
VlnPlot(SR_WL_20251202, features = c("percent.ribo"))

# Cluster markers
SR_WL_20251202_clustermarkers_df <- as.data.frame(FindAllMarkers(SR_WL_20251202)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(SR_WL_20251202_clustermarkers_df)) {
  killifish_NCBI_gene_name = SR_WL_20251202_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {SR_WL_20251202_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {SR_WL_20251202_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {SR_WL_20251202_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {SR_WL_20251202_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {SR_WL_20251202_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {SR_WL_20251202_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Save and reload RDS
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(SR_WL_20251202, file = "SR_WL_20251202.rds")
#SR_WL_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/SR_WL_20251202.rds")