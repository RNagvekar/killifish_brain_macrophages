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
# YSR (red cells from young SP-oScarlet fish, August 2025 experiment)
# ----------------------------------------------------------------------------------------------------

# Create Seurat object
YSR_20251202.data <- Read10X(data.dir = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/YSR/raw_feature_bc_matrix")
YSR_20251202 <- CreateSeuratObject(counts = YSR_20251202.data, project = "20251202_10x", min.cells = 3, min.features = 200)

# Orig.ident
YSR_20251202$orig.ident <- "YSR"

# Percent mitochondrial reads
YSR_20251202[["percent.mt"]] <- PercentageFeatureSet(YSR_20251202, pattern = "^MT-")

# QC plots
VlnPlot(YSR_20251202, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
FeatureScatter(YSR_20251202, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtering
YSR_20251202 <- subset(YSR_20251202, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

# Normalization
YSR_20251202 <- NormalizeData(YSR_20251202, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
YSR_20251202 <- FindVariableFeatures(YSR_20251202, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(YSR_20251202)
YSR_20251202 <- ScaleData(YSR_20251202, features = all.genes)

# PCA
YSR_20251202 <- RunPCA(YSR_20251202, features = VariableFeatures(object = YSR_20251202))

# UMAP
YSR_20251202 <- RunUMAP(YSR_20251202, dims = 1:30)
YSR_20251202 <- FindNeighbors(YSR_20251202, dims = 1:30)
YSR_20251202 <- FindClusters(YSR_20251202, resolution = 0.2)
DimPlot(YSR_20251202, reduction = "umap", label = TRUE)

# Feature plot
FeaturePlot(YSR_20251202, "LOC107383266", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(YSR_20251202, features = c("sall3"))

# QC metrics
VlnPlot(YSR_20251202, features = c("percent.mt")) 
VlnPlot(YSR_20251202, features = c("nFeature_RNA"))
VlnPlot(YSR_20251202, features = c("nCount_RNA"))
YSR_20251202[["percent.ribo"]] <- PercentageFeatureSet(YSR_20251202, pattern = "^rps|^rpl")
VlnPlot(YSR_20251202, features = c("percent.ribo"))

# Cluster markers
YSR_20251202_clustermarkers_df <- as.data.frame(FindAllMarkers(YSR_20251202)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(YSR_20251202_clustermarkers_df)) {
  killifish_NCBI_gene_name = YSR_20251202_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {YSR_20251202_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {YSR_20251202_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_20251202_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {YSR_20251202_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_20251202_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {YSR_20251202_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Save and reload RDS
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(YSR_20251202, file = "YSR_20251202.rds")
YSR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/YSR_20251202.rds")

# ----------------------------------------------------------------------------------------------------
# YSNR (non-red cells from young SP-oScarlet fish, August 2025 experiment)
# ----------------------------------------------------------------------------------------------------

# Create Seurat object
YSNR_20251202.data <- Read10X(data.dir = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/YSNR/raw_feature_bc_matrix")
YSNR_20251202 <- CreateSeuratObject(counts = YSNR_20251202.data, project = "20251202_10x", min.cells = 3, min.features = 200)

# Orig.ident
YSNR_20251202$orig.ident <- "YSNR"

# Percent mitochondrial reads
YSNR_20251202[["percent.mt"]] <- PercentageFeatureSet(YSNR_20251202, pattern = "^MT-")

# QC plots
VlnPlot(YSNR_20251202, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
FeatureScatter(YSNR_20251202, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtering
YSNR_20251202 <- subset(YSNR_20251202, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

# Normalization
YSNR_20251202 <- NormalizeData(YSNR_20251202, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
YSNR_20251202 <- FindVariableFeatures(YSNR_20251202, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(YSNR_20251202)
YSNR_20251202 <- ScaleData(YSNR_20251202, features = all.genes)

# PCA
YSNR_20251202 <- RunPCA(YSNR_20251202, features = VariableFeatures(object = YSNR_20251202))

# UMAP
YSNR_20251202 <- RunUMAP(YSNR_20251202, dims = 1:30)
YSNR_20251202 <- FindNeighbors(YSNR_20251202, dims = 1:30)
YSNR_20251202 <- FindClusters(YSNR_20251202, resolution = 0.5)
DimPlot(YSNR_20251202, reduction = "umap", label = TRUE)

# Feature plot
FeaturePlot(YSNR_20251202, "eSPOwrong-insertion", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(YSNR_20251202, features = c("csf1r"))

# QC metrics
VlnPlot(YSNR_20251202, features = c("percent.mt")) 
VlnPlot(YSNR_20251202, features = c("nFeature_RNA"))
VlnPlot(YSNR_20251202, features = c("nCount_RNA"))
YSNR_20251202[["percent.ribo"]] <- PercentageFeatureSet(YSNR_20251202, pattern = "^rps|^rpl")
VlnPlot(YSNR_20251202, features = c("percent.ribo"))

# Cluster markers
YSNR_20251202_clustermarkers_df <- as.data.frame(FindAllMarkers(YSNR_20251202)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(YSNR_20251202_clustermarkers_df)) {
  killifish_NCBI_gene_name = YSNR_20251202_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {YSNR_20251202_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {YSNR_20251202_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSNR_20251202_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {YSNR_20251202_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSNR_20251202_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {YSNR_20251202_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Save and reload RDS
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(YSNR_20251202, file = "YSNR_20251202.rds")
YSNR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/YSNR_20251202.rds")

# ----------------------------------------------------------------------------------------------------
# SRR24058857 (young wildtype telencephalon sample from Ayana et al., 2024)
# N.B. in NCBI GEO/SRA there are two young wildtype telecephalon samples from this paper
# SRR24058857 is the one with more high-quality cells; the other is SRR14110665
# ----------------------------------------------------------------------------------------------------

# Create Seurat object
SRR24058857_20260315.data <- Read10X(data.dir = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/SRR24058857/raw_feature_bc_matrix")
SRR24058857_20260315 <- CreateSeuratObject(counts = SRR24058857_20260315.data, project = "20260315_10x", min.cells = 3, min.features = 200)

# Orig.ident
SRR24058857_20260315$orig.ident <- "SRR24058857"

# Percent mitochondrial reads
SRR24058857_20260315[["percent.mt"]] <- PercentageFeatureSet(SRR24058857_20260315, pattern = "^MT-")

# QC plots
VlnPlot(SRR24058857_20260315, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
FeatureScatter(SRR24058857_20260315, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtering
SRR24058857_20260315 <- subset(SRR24058857_20260315, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

# Normalization
SRR24058857_20260315 <- NormalizeData(SRR24058857_20260315, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
SRR24058857_20260315 <- FindVariableFeatures(SRR24058857_20260315, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(SRR24058857_20260315)
SRR24058857_20260315 <- ScaleData(SRR24058857_20260315, features = all.genes)

# PCA
SRR24058857_20260315 <- RunPCA(SRR24058857_20260315, features = VariableFeatures(object = SRR24058857_20260315))

# UMAP
SRR24058857_20260315 <- RunUMAP(SRR24058857_20260315, dims = 1:30)
SRR24058857_20260315 <- FindNeighbors(SRR24058857_20260315, dims = 1:30)
SRR24058857_20260315 <- FindClusters(SRR24058857_20260315, resolution = 0.5)
DimPlot(SRR24058857_20260315, reduction = "umap", label = TRUE)

# Feature plot
FeaturePlot(SRR24058857_20260315, "lyve1")

# Violin plot
VlnPlot(SRR24058857_20260315, features = c("csf1r"))

# QC metrics
VlnPlot(SRR24058857_20260315, features = c("percent.mt")) 
VlnPlot(SRR24058857_20260315, features = c("nFeature_RNA"))
VlnPlot(SRR24058857_20260315, features = c("nCount_RNA"))
SRR24058857_20260315[["percent.ribo"]] <- PercentageFeatureSet(SRR24058857_20260315, pattern = "^rps|^rpl")
VlnPlot(SRR24058857_20260315, features = c("percent.ribo"))

# Cluster markers
SRR24058857_20260315_clustermarkers_df <- as.data.frame(FindAllMarkers(SRR24058857_20260315)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(SRR24058857_20260315_clustermarkers_df)) {
  killifish_NCBI_gene_name = SRR24058857_20260315_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {SRR24058857_20260315_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {SRR24058857_20260315_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {SRR24058857_20260315_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {SRR24058857_20260315_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {SRR24058857_20260315_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {SRR24058857_20260315_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Save and reload RDS
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(SRR24058857_20260315, file = "SRR24058857_20260315.rds")
SRR24058857_20260315 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/SRR24058857_20260315.rds")

# ----------------------------------------------------------------------------------------------------
# OSR (red cells from old SP-oScarlet fish, August 2025 experiment)
# ----------------------------------------------------------------------------------------------------

# Create Seurat object
OSR_20251202.data <- Read10X(data.dir = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/OSR/raw_feature_bc_matrix")
OSR_20251202 <- CreateSeuratObject(counts = OSR_20251202.data, project = "20251202_10x", min.cells = 3, min.features = 200)

# Orig.ident
OSR_20251202$orig.ident <- "OSR"

# Percent mitochondrial reads
OSR_20251202[["percent.mt"]] <- PercentageFeatureSet(OSR_20251202, pattern = "^MT-")

# QC plots
VlnPlot(OSR_20251202, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
FeatureScatter(OSR_20251202, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtering
OSR_20251202 <- subset(OSR_20251202, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

# Normalization
OSR_20251202 <- NormalizeData(OSR_20251202, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
OSR_20251202 <- FindVariableFeatures(OSR_20251202, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(OSR_20251202)
OSR_20251202 <- ScaleData(OSR_20251202, features = all.genes)

# PCA
OSR_20251202 <- RunPCA(OSR_20251202, features = VariableFeatures(object = OSR_20251202))

# UMAP
OSR_20251202 <- RunUMAP(OSR_20251202, dims = 1:30)
OSR_20251202 <- FindNeighbors(OSR_20251202, dims = 1:30)
OSR_20251202 <- FindClusters(OSR_20251202, resolution = 0.5)
DimPlot(OSR_20251202, reduction = "umap", label = TRUE)

# Feature plot
FeaturePlot(OSR_20251202, "LOC107375757", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(OSR_20251202, features = c("LOC107375757"))

# QC metrics
VlnPlot(OSR_20251202, features = c("percent.mt")) 
VlnPlot(OSR_20251202, features = c("nFeature_RNA"))
VlnPlot(OSR_20251202, features = c("nCount_RNA"))
OSR_20251202[["percent.ribo"]] <- PercentageFeatureSet(OSR_20251202, pattern = "^rps|^rpl")
VlnPlot(OSR_20251202, features = c("percent.ribo"))

# Cluster markers
OSR_20251202_clustermarkers_df <- as.data.frame(FindAllMarkers(OSR_20251202)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(OSR_20251202_clustermarkers_df)) {
  killifish_NCBI_gene_name = OSR_20251202_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {OSR_20251202_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {OSR_20251202_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {OSR_20251202_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {OSR_20251202_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {OSR_20251202_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {OSR_20251202_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Save and reload RDS
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(OSR_20251202, file = "OSR_20251202.rds")
OSR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/OSR_20251202.rds")

# ----------------------------------------------------------------------------------------------------
# SR (red cells from old SP-oScarlet fish, October 2024 experiment)
# ----------------------------------------------------------------------------------------------------

# Create Seurat object
SR_20251202.data <- Read10X(data.dir = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/SR/raw_feature_bc_matrix")
SR_20251202 <- CreateSeuratObject(counts = SR_20251202.data, project = "20251202_10x", min.cells = 3, min.features = 200)

# Orig.ident
SR_20251202$orig.ident <- "SR"

# Percent mitochondrial reads
SR_20251202[["percent.mt"]] <- PercentageFeatureSet(SR_20251202, pattern = "^MT-")

# QC plots
VlnPlot(SR_20251202, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
FeatureScatter(SR_20251202, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtering
SR_20251202 <- subset(SR_20251202, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

# Normalization
SR_20251202 <- NormalizeData(SR_20251202, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
SR_20251202 <- FindVariableFeatures(SR_20251202, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(SR_20251202)
SR_20251202 <- ScaleData(SR_20251202, features = all.genes)

# PCA
SR_20251202 <- RunPCA(SR_20251202, features = VariableFeatures(object = SR_20251202))

# UMAP
SR_20251202 <- RunUMAP(SR_20251202, dims = 1:30)
SR_20251202 <- FindNeighbors(SR_20251202, dims = 1:30)
SR_20251202 <- FindClusters(SR_20251202, resolution = 0.5)
DimPlot(SR_20251202, reduction = "umap", label = TRUE)

# Feature plot
FeaturePlot(SR_20251202, "sall3", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(SR_20251202, features = c("LOC107386542"))

# QC metrics
VlnPlot(SR_20251202, features = c("percent.mt")) 
VlnPlot(SR_20251202, features = c("nFeature_RNA"))
VlnPlot(SR_20251202, features = c("nCount_RNA"))
SR_20251202[["percent.ribo"]] <- PercentageFeatureSet(SR_20251202, pattern = "^rps|^rpl")
VlnPlot(SR_20251202, features = c("percent.ribo"))

# Cluster markers
SR_20251202_clustermarkers_df <- as.data.frame(FindAllMarkers(SR_20251202)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(SR_20251202_clustermarkers_df)) {
  killifish_NCBI_gene_name = SR_20251202_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {SR_20251202_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {SR_20251202_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {SR_20251202_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {SR_20251202_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {SR_20251202_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {SR_20251202_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Save and reload RDS
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(SR_20251202, file = "SR_20251202.rds")
SR_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/SR_20251202.rds")

# ----------------------------------------------------------------------------------------------------
# WL (live cells from old wildtype fish, October 2024 experiment)
# ----------------------------------------------------------------------------------------------------

# Create Seurat object
WL_20251202.data <- Read10X(data.dir = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/WL/raw_feature_bc_matrix")
WL_20251202 <- CreateSeuratObject(counts = WL_20251202.data, project = "20251202_10x", min.cells = 3, min.features = 200)

# Orig.ident
WL_20251202$orig.ident <- "WL"

# Percent mitochondrial reads
WL_20251202[["percent.mt"]] <- PercentageFeatureSet(WL_20251202, pattern = "^MT-")

# QC plots
VlnPlot(WL_20251202, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
FeatureScatter(WL_20251202, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtering
WL_20251202 <- subset(WL_20251202, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

# Normalization
WL_20251202 <- NormalizeData(WL_20251202, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
WL_20251202 <- FindVariableFeatures(WL_20251202, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(WL_20251202)
WL_20251202 <- ScaleData(WL_20251202, features = all.genes)

# PCA
WL_20251202 <- RunPCA(WL_20251202, features = VariableFeatures(object = WL_20251202))

# UMAP
WL_20251202 <- RunUMAP(WL_20251202, dims = 1:30)
WL_20251202 <- FindNeighbors(WL_20251202, dims = 1:30)
WL_20251202 <- FindClusters(WL_20251202, resolution = 0.5)
DimPlot(WL_20251202, reduction = "umap", label = TRUE)

# Feature plot
FeaturePlot(WL_20251202, "cldn5", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(WL_20251202, features = c("csf1r"))

# QC metrics
VlnPlot(WL_20251202, features = c("percent.mt")) 
VlnPlot(WL_20251202, features = c("nFeature_RNA"))
VlnPlot(WL_20251202, features = c("nCount_RNA"))
WL_20251202[["percent.ribo"]] <- PercentageFeatureSet(WL_20251202, pattern = "^rps|^rpl")
VlnPlot(WL_20251202, features = c("elavl3"))

# Cluster markers
WL_20251202_clustermarkers_df <- as.data.frame(FindAllMarkers(WL_20251202)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(WL_20251202_clustermarkers_df)) {
  killifish_NCBI_gene_name = WL_20251202_clustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {WL_20251202_clustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {WL_20251202_clustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {WL_20251202_clustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {WL_20251202_clustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {WL_20251202_clustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {WL_20251202_clustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Doublet finder (heuristic - no doublet removal used in paper)
optimal_pK <- 0.10
annotations <- WL_20251202@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(optimal_pK * nrow(WL_20251202@meta.data))
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
WL_20251202 <- doubletFinder(seu = WL_20251202, PCs = 1:30, pK = optimal_pK, nExp = nExp.poi.adj)

# Singlets vs. doublets plot
DimPlot(WL_20251202, reduction = "umap", group.by = "DF.classifications_0.25_0.1_1654")

# Save and reload RDS
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(WL_20251202, file = "WL_20251202.rds")
WL_20251202 <- readRDS(file = "C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/WL_20251202.rds")