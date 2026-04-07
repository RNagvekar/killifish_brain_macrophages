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

# Read Seurat object (YSR_YSNR)
YSR_YSNR_20251202 <- readRDS(file = "YSR_YSNR_20251202.rds")

# Extract UMAP coordinates from the Seurat object
YSR_YSNR_umap_coords <- Embeddings(YSR_YSNR_20251202, "umap")
YSR_YSNR_xlims <- range(YSR_YSNR_umap_coords[,1])
YSR_YSNR_ylims <- range(YSR_YSNR_umap_coords[,2])

# Fig 1f
# UMAP by orig.ident
# YSNR = grey, YSR = red
# Save at 500W x 500H
DimPlot(YSR_YSNR_20251202, reduction = "umap", group.by = "orig.ident", cols = c("grey", "red")) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# UMAP with clustering at res = 0.1
YSR_YSNR_20251202 <- FindClusters(YSR_YSNR_20251202, resolution = 0.1)
DimPlot(YSR_YSNR_20251202, reduction = "umap", label = TRUE, label.size = 10) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Cluster markers, res = 0.1
YSR_YSNR_20251202_lowresclustermarkers_df <- as.data.frame(FindAllMarkers(YSR_YSNR_20251202)) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(p_val_adj < 0.05)

# Add human, mouse, and zebrafish gene names for cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(YSR_YSNR_20251202_lowresclustermarkers_df)) {
  killifish_NCBI_gene_name = YSR_YSNR_20251202_lowresclustermarkers_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {YSR_YSNR_20251202_lowresclustermarkers_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {YSR_YSNR_20251202_lowresclustermarkers_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_YSNR_20251202_lowresclustermarkers_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {YSR_YSNR_20251202_lowresclustermarkers_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_YSNR_20251202_lowresclustermarkers_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {YSR_YSNR_20251202_lowresclustermarkers_df$zebrafish_gene_name[i] = ""}
}

# Define cell types (manual - based on markers and visual inspection of UMAP)
YSR_YSNR_20251202$celltype <- case_when(
  Idents(YSR_YSNR_20251202) %in% c(0) ~ "Myeloid cells",
  Idents(YSR_YSNR_20251202) %in% c(1) ~ "Radial glia/progenitors/other",
  Idents(YSR_YSNR_20251202) %in% c(2, 3, 5) ~ "Neurons",
  Idents(YSR_YSNR_20251202) %in% c(4) ~ "Oligodendrocyte precursor cells",
  Idents(YSR_YSNR_20251202) %in% c(6) ~ "Oligodendrocytes",
  Idents(YSR_YSNR_20251202) %in% c(7) ~ "Pericytes",
  Idents(YSR_YSNR_20251202) %in% c(8) ~ "T cells")
table(YSR_YSNR_20251202$orig.ident, YSR_YSNR_20251202$celltype)

# Fig 1S2a
# Heatmap of celltype markers
# Save at 1000W x 600H
DoHeatmap(YSR_YSNR_20251202,
          features = top5markers_df$feature,
          lines.width = 10,
          group.by = "celltype",
          label = FALSE) +
  scale_y_discrete(labels = function(x) toupper(x)) +
  theme(axis.text.y = element_text(size = 10, face = "italic"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

# Fig 1g
# UMAP by celltype
# Save at 500W x 500H
DimPlot(YSR_YSNR_20251202, group.by = "celltype", reduction = "umap", label = FALSE) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 1g
# UMAP by celltype with legend
# Save at 700W x 500H
DimPlot(YSR_YSNR_20251202, group.by = "celltype", reduction = "umap", label = FALSE) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Cell type markers by AUC
YSR_YSNR_20251202_data <- GetAssayData(YSR_YSNR_20251202, layer = "data")
YSR_YSNR_20251202_celltypes <- YSR_YSNR_20251202$celltype
YSR_YSNR_20251202_celltypemarkers_AUC_df <- as.data.frame(wilcoxauc(YSR_YSNR_20251202_data, YSR_YSNR_20251202_celltypes))

top5markers_df <- YSR_YSNR_20251202_celltypemarkers_AUC_df %>% 
  group_by(group) %>%
  top_n(n = 5, wt = auc)

# Violin plots for myeloid markers

# Fig 1S2b
# apoeb
# Save at 700W x 500H
VlnPlot(YSR_YSNR_20251202, features = c("LOC107379395"), group.by = "celltype", pt.size = 0) +
  ggtitle(NULL) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$celltype),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Fig 1S2b
# iba1
# Save at 700W x 500H
VlnPlot(YSR_YSNR_20251202, features = c("LOC107378674"), group.by = "celltype", pt.size = 0) +
  ggtitle(NULL) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$celltype),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Fig 1S2b
# lcp1
# Save at 700W x 500H
VlnPlot(YSR_YSNR_20251202, features = c("lcp1"), group.by = "celltype", pt.size = 0) +
  ggtitle(NULL) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$celltype),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Fig 1S2b
# csf1r
# Save at 700W x 500H
VlnPlot(YSR_YSNR_20251202, features = c("csf1r"), group.by = "celltype", pt.size = 0) +
  ggtitle(NULL) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$celltype),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Fig 1S2b
# csf1r with labeled x-axis
# Save at 700W x 1000H
VlnPlot(YSR_YSNR_20251202, features = c("csf1r"), group.by = "celltype", pt.size = 0) +
  ggtitle(NULL) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$celltype),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 25, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Percent myeloid
percent_myeloid_YSR_YSNR_20251202_df <- YSR_YSNR_20251202@meta.data %>%
  group_by(orig.ident) %>%
  summarise(percent_myeloid = mean(celltype == "Myeloid cells") * 100)
ggplot(percent_myeloid_YSR_YSNR_20251202_df,
       aes(x = orig.ident, y = percent_myeloid, color = orig.ident)) +
  ylim(c(0, 100)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("YSNR" = "grey", "YSR" = "red")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line.x = element_line(color = "black", size = 0.25),
        axis.line.y = element_line(color = "black", size = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        legend.position = "none")

# Fig 1h
# csf1r feature plot
# Save at 500W x 500H
FeaturePlot(YSR_YSNR_20251202, "csf1r", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 1h
# csf1r feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_YSNR_20251202, "csf1r", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# apoeb feature plot
# Save at 500W x 500H
FeaturePlot(YSR_YSNR_20251202, "LOC107379395",  cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# apoeb feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_YSNR_20251202, "LOC107379395",  cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Fig 1i
# elavl3 feature plot
# Save at 500W x 500H
FeaturePlot(YSR_YSNR_20251202, "elavl3", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 1i
# elavl3 feature plot with legend
# Save at 600W x 500H
FeaturePlot(YSR_YSNR_20251202, "elavl3", cols = c("yellow", "purple")) +
  coord_cartesian(xlim = YSR_YSNR_xlims, ylim = YSR_YSNR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

# Save supplemental tables

# YSR_YSNR cluster markers (all, resolution = 0.1)
YSR_YSNR_20251202 <- FindClusters(YSR_YSNR_20251202, resolution = 0.1)
YSR_YSNR_20251202_lowresclustermarkers_all_df <- as.data.frame(FindAllMarkers(YSR_YSNR_20251202))

# Add human, mouse, and zebrafish gene names for all cluster markers
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(YSR_YSNR_20251202_lowresclustermarkers_all_df)) {
  killifish_NCBI_gene_name = YSR_YSNR_20251202_lowresclustermarkers_all_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {YSR_YSNR_20251202_lowresclustermarkers_all_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {YSR_YSNR_20251202_lowresclustermarkers_all_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_YSNR_20251202_lowresclustermarkers_all_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {YSR_YSNR_20251202_lowresclustermarkers_all_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {YSR_YSNR_20251202_lowresclustermarkers_all_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {YSR_YSNR_20251202_lowresclustermarkers_all_df$zebrafish_gene_name[i] = ""}
}

# Add celltypes
YSR_YSNR_20251202_lowresclustermarkers_all_df$celltype <- case_when(
  YSR_YSNR_20251202_lowresclustermarkers_all_df$cluster %in% c(0) ~ "Myeloid cells",
  YSR_YSNR_20251202_lowresclustermarkers_all_df$cluster %in% c(1) ~ "Radial glia/progenitors/other",
  YSR_YSNR_20251202_lowresclustermarkers_all_df$cluster %in% c(2, 3, 5) ~ "Neurons",
  YSR_YSNR_20251202_lowresclustermarkers_all_df$cluster %in% c(4) ~ "Oligodendrocyte precursor cells",
  YSR_YSNR_20251202_lowresclustermarkers_all_df$cluster %in% c(6) ~ "Oligodendrocytes",
  YSR_YSNR_20251202_lowresclustermarkers_all_df$cluster %in% c(7) ~ "Pericytes",
  YSR_YSNR_20251202_lowresclustermarkers_all_df$cluster %in% c(8) ~ "T cells")

# Save
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/")
write.table(YSR_YSNR_20251202_lowresclustermarkers_all_df, "SupplementaryTable7.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# YSR_YSNR cell type markers (by AUC)
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/")
write.table(YSR_YSNR_20251202_celltypemarkers_AUC_df, "SupplementaryTable8.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Save Seurat objects
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(YSR_YSNR_20251202, file = "YSR_YSNR_20251202.rds")