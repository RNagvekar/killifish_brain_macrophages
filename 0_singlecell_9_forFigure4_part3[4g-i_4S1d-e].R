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

# Read Seurat object (SR_WL)
SR_WL_20251202 <- readRDS(file = "SR_WL_20251202.rds")

# Extract UMAP coordinates from the Seurat object
SR_WL_umap_coords <- Embeddings(SR_WL_20251202, "umap")
SR_WL_xlims <- range(SR_WL_umap_coords[,1])
SR_WL_ylims <- range(SR_WL_umap_coords[,2])

# UMAP by orig.ident
# SR = red, WL = grey
# Save at 500W x 500H
DimPlot(SR_WL_20251202, reduction = "umap", group.by = "orig.ident", cols = c("red", "grey")) +
  coord_cartesian(xlim = SR_WL_xlims, ylim = SR_WL_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# UMAP with clustering at res = 0.2
# Save at 500W x 500H
SR_WL_20251202 <- FindClusters(SR_WL_20251202, resolution = 0.2)
DimPlot(SR_WL_20251202, reduction = "umap", label = TRUE, label.size = 8) +
  coord_cartesian(xlim = SR_WL_xlims, ylim = SR_WL_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Violin plot
VlnPlot(SR_WL_20251202, features = c("lcp1")) +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95)

# Percent red by cluster
# Save plot at 1000W x 500H
percent_red_df <- SR_WL_20251202@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    total = n(),
    n_red = sum(orig.ident == "SR"),
    pct_red = 100 * n_red / total
  )
ggplot(percent_red_df, aes(x = seurat_clusters, y = pct_red, color = "red")) +
  ylim(c(0, 100)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none")

# csf1r violin plot
# Save at 1000W x 500H
VlnPlot(SR_WL_20251202, features = c("csf1r")) +
  ggtitle(NULL) +
  theme(
    axis.text = element_text(size = 20),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# lyve1 violin plot
# Save at 1000W x 500H
VlnPlot(SR_WL_20251202, features = c("lyve1")) +
  ggtitle(NULL) +
  theme(
    axis.text = element_text(size = 20),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Subset to myeloid (csf1r+) clusters
SR_WL_myeloid_20251202 <- subset(SR_WL_20251202, idents = c(0, 1, 12))

# Redo-processing on subset object (myeloid only)

# Normalization
SR_WL_myeloid_20251202 <- NormalizeData(SR_WL_myeloid_20251202, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
SR_WL_myeloid_20251202 <- FindVariableFeatures(SR_WL_myeloid_20251202, selection.method = "vst", nfeatures = 2000)

# Scaling 
all.genes <- rownames(SR_WL_myeloid_20251202)
SR_WL_myeloid_20251202 <- ScaleData(SR_WL_myeloid_20251202, features = all.genes)

# PCA
SR_WL_myeloid_20251202 <- RunPCA(SR_WL_myeloid_20251202, features = VariableFeatures(object = SR_WL_myeloid_20251202))

# UMAP
SR_WL_myeloid_20251202 <- RunUMAP(SR_WL_myeloid_20251202, dims = 1:30)
SR_WL_myeloid_20251202 <- FindNeighbors(SR_WL_myeloid_20251202, dims = 1:30)
SR_WL_myeloid_20251202 <- FindClusters(SR_WL_myeloid_20251202, resolution = 0.5)
DimPlot(SR_WL_myeloid_20251202, reduction = "umap", label = TRUE)
DimPlot(SR_WL_myeloid_20251202, reduction = "umap", group.by = "orig.ident")

# Feature plot
FeaturePlot(SR_WL_myeloid_20251202, "LOC107386542", cols = c("yellow", "purple"))

# Violin plot
VlnPlot(SR_WL_myeloid_20251202, features = c("LOC107373603"))
VlnPlot(SR_WL_myeloid_20251202, features = c("LOC107383768"), group.by = "orig.ident")

# QC metrics
VlnPlot(SR_WL_myeloid_20251202, features = c("percent.mt")) 
VlnPlot(SR_WL_myeloid_20251202, features = c("nFeature_RNA"))
VlnPlot(SR_WL_myeloid_20251202, features = c("nCount_RNA"))
SR_WL_myeloid_20251202[["percent.ribo"]] <- PercentageFeatureSet(YSR_OSR_20251202, pattern = "^rps|^rpl")
VlnPlot(SR_WL_myeloid_20251202, features = c("percent.ribo"))

# Cells by cluster by orig.ident
table(Idents(SR_WL_myeloid_20251202), SR_WL_myeloid_20251202$orig.ident)

# Read Seurat object (SR_WL_myeloid)
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202")
SR_WL_myeloid_20251202 <- readRDS(file = "SR_WL_myeloid_20251202.rds")

# Extract UMAP coordinates from the Seurat object
SR_WL_myeloid_umap_coords <- Embeddings(SR_WL_myeloid_20251202, "umap")
SR_WL_myeloid_xlims <- range(SR_WL_myeloid_umap_coords[,1])
SR_WL_myeloid_ylims <- range(SR_WL_myeloid_umap_coords[,2])

# N.B. For Fig 4f please refer to script 0_singlecell_5_forFigure2_part2[2c_2S1b]_forFigure4_part2[4f_4S1c].R

# Fig 4g
# UMAP by orig.ident
# SR = purple, WL = grey
# Save at 500W x 500H
DimPlot(SR_WL_myeloid_20251202, reduction = "umap", group.by = "orig.ident", cols = c("purple", "darkgrey")) +
  coord_cartesian(xlim = SR_WL_myeloid_xlims, ylim = SR_WL_myeloid_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 4h
# Dotplot of marker expression in SR_myeloid (above) and WL_myeloid (below)
# Save at 500W x 400H
# LOC107379395 = apoeb
# LOC107372638 = HLA-DPB1 (MHCII)
# LOC107385068 = c1qc
# LOC107383768 = mrc1
markers <- c("csf1r",
             "LOC107379395",
             "LOC107372638",
             "LOC107385068",
             "LOC107383768")
dotplot_SR_WL_myeloid <- DotPlot(object = SR_WL_myeloid_20251202,
                                 features = markers,
                                 group.by = "orig.ident",
                                 scale = FALSE)
ggplot(dotplot_SR_WL_myeloid$data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_y_discrete(limits = rev) +
  scale_color_gradient(
    low = "yellow",
    high = "purple",
    limits = c(0, 7),
    trans = "sqrt",
    oob = scales::squish,
    na.value = "white",
    name = "Average\nlog-normalized expression") +
  scale_size_continuous(
    limits = c(0, 100),
    range  = c(1, 10),
    name   = "Percentage of\ncells with expression") +
  guides(
    color = guide_colorbar(order = 1, direction = "vertical"),
    size = guide_legend(order = 2)) +
  theme(legend.position = "right") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line.x = element_line(color = "black", size = 0.25),
        axis.line.y = element_line(color = "black", size = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.text = element_text(size = 12))

# DEG between SRmyeloid and WLmyeloid
DEG_SR_WL_myeloid_20251202 <- FindMarkers(
  object = SR_WL_myeloid_20251202,
  ident.1 = "SR",
  ident.2 = "WL",
  group.by = "orig.ident",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.5
)
DEG_SR_WL_myeloid_20251202_df <- as.data.frame(DEG_SR_WL_myeloid_20251202)
DEG_SR_WL_myeloid_20251202_df$gene <- rownames(DEG_SR_WL_myeloid_20251202_df)

# Define background (human genes)
background_SR_WL_myeloid_killifish <- rownames(SR_WL_myeloid_20251202)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
background_SR_WL_myeloid_human <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$`N. furzeri (NCBI)`[i] %in% background_SR_WL_myeloid_killifish) {background_SR_WL_myeloid_human <- append(background_SR_WL_myeloid_human, gene_names_df$Human[i])}
}

# Add human, mouse, and zebrafish gene names for DEGs
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(DEG_SR_WL_myeloid_20251202_df)) {
  killifish_NCBI_gene_name = DEG_SR_WL_myeloid_20251202_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {DEG_SR_WL_myeloid_20251202_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {DEG_SR_WL_myeloid_20251202_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {DEG_SR_WL_myeloid_20251202_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {DEG_SR_WL_myeloid_20251202_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {DEG_SR_WL_myeloid_20251202_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {DEG_SR_WL_myeloid_20251202_df$zebrafish_gene_name[i] = ""}
}

# Enriched in each sample
SRmyeloid_enriched_vsWLmyeloid_DEG_df <- DEG_SR_WL_myeloid_20251202_df %>% dplyr::filter(avg_log2FC > 0.5) %>% dplyr::filter(pct.1 > 0.1) %>% dplyr::filter(p_val_adj < 0.05)
WLmyeloid_enriched_vsSRmyeloid_DEG_df <- DEG_SR_WL_myeloid_20251202_df %>% dplyr::filter(avg_log2FC < -0.5) %>% dplyr::filter(pct.2 > 0.1) %>% dplyr::filter(p_val_adj < 0.05)

# GO terms on genes enriched in SRmyeloid (relative to WLmyeloid)
SRmyeloid_enriched_vsWLmyeloid_GOterms_BP <- enrichr(SRmyeloid_enriched_vsWLmyeloid_DEG_df$human_gene_name, "GO_Biological_Process_2025", background = background_SR_WL_myeloid_human)
SRmyeloid_enriched_vsWLmyeloid_GOterms_CC <- enrichr(SRmyeloid_enriched_vsWLmyeloid_DEG_df$human_gene_name, "GO_Cellular_Component_2025", background = background_SR_WL_myeloid_human)
SRmyeloid_enriched_vsWLmyeloid_GOterms_MF <- enrichr(SRmyeloid_enriched_vsWLmyeloid_DEG_df$human_gene_name, "GO_Molecular_Function_2025", background = background_SR_WL_myeloid_human)

SRmyeloid_enriched_vsWLmyeloid_GOterms_BP_df <- as.data.frame(SRmyeloid_enriched_vsWLmyeloid_GOterms_BP)
SRmyeloid_enriched_vsWLmyeloid_GOterms_BP_df$Term <- gsub(" \\(GO:[0-9]+\\)", "", SRmyeloid_enriched_vsWLmyeloid_GOterms_BP_df$GO_Biological_Process_2025.Term)
SRmyeloid_enriched_vsWLmyeloid_GOterms_BP_df$n_genes <- lengths(strsplit(as.character(SRmyeloid_enriched_vsWLmyeloid_GOterms_BP_df$GO_Biological_Process_2025.Genes), ";"))

# GO terms on genes enriched in WLmyeloid (relative to SRmyeloid)
WLmyeloid_enriched_vsSRmyeloid_GOterms_BP <- enrichr(WLmyeloid_enriched_vsSRmyeloid_DEG_df$human_gene_name, "GO_Biological_Process_2025", background = background_SR_WL_myeloid_human)
WLmyeloid_enriched_vsSRmyeloid_GOterms_CC <- enrichr(WLmyeloid_enriched_vsSRmyeloid_DEG_df$human_gene_name, "GO_Cellular_Component_2025", background = background_SR_WL_myeloid_human)
WLmyeloid_enriched_vsSRmyeloid_GOterms_MF <- enrichr(WLmyeloid_enriched_vsSRmyeloid_DEG_df$human_gene_name, "GO_Molecular_Function_2025", background = background_SR_WL_myeloid_human)

WLmyeloid_enriched_vsSRmyeloid_GOterms_BP_df <- as.data.frame(WLmyeloid_enriched_vsSRmyeloid_GOterms_BP)
WLmyeloid_enriched_vsSRmyeloid_GOterms_BP_df$Term <- gsub(" \\(GO:[0-9]+\\)", "", WLmyeloid_enriched_vsSRmyeloid_GOterms_BP_df$GO_Biological_Process_2025.Term)
WLmyeloid_enriched_vsSRmyeloid_GOterms_BP_df$n_genes <- lengths(strsplit(as.character(WLmyeloid_enriched_vsSRmyeloid_GOterms_BP_df$GO_Biological_Process_2025.Genes), ";"))

# Fig 4S1c
# Plot for GO terms on genes enriched in SRmyeloid (relative to WLmyeloid)
# Save at 1000W x 500H
ggplot(data = subset(SRmyeloid_enriched_vsWLmyeloid_GOterms_BP_df, GO_Biological_Process_2025.Adjusted.P.value < 1E-8),
       aes(x = -log10(GO_Biological_Process_2025.Adjusted.P.value),
           y = reorder(Term, -log10(GO_Biological_Process_2025.Adjusted.P.value)),
           size = n_genes)) +
  geom_line(linetype = "blank") +
  geom_point(color = "purple") +
  scale_size(range = c(16, 17)) +
  xlim(0, 17) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line.x = element_line(color = "black", size = 0.25),
        axis.line.y = element_line(color = "black", size = 0.25),
        axis.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        legend.position = "none")

# N.B. For Fig 4S1c please refer to script 0_singlecell_5_forFigure2_part2[2c_2S1b]_forFigure4_part2[4f_4S1c].R

# Fig 4S1d
# Plot for GO terms on genes enriched in WLmyeloid (relative to SRmyeloid)
# Save at 1200W x 500H
ggplot(data = subset(WLmyeloid_enriched_vsSRmyeloid_GOterms_BP_df, GO_Biological_Process_2025.Adjusted.P.value < 2E-3),
       aes(x = -log10(GO_Biological_Process_2025.Adjusted.P.value),
           y = reorder(Term, -log10(GO_Biological_Process_2025.Adjusted.P.value)),
           size = n_genes)) +
  geom_line(linetype = "blank") +
  geom_point(color = "grey") +
  scale_size(range = c(6, 10)) +
  xlim(0, 4) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line.x = element_line(color = "black", size = 0.25),
        axis.line.y = element_line(color = "black", size = 0.25),
        axis.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        legend.position = "none")

# Fig 4i
# Violin plot for GO term "cytoplasmic translation"
# SRmyeloid = purple, WLmyeloid = grey
# Save at 300W x 500H
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
cytoplasmictranslation_GO <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "go_id"), filters = "go", values = "GO:0002181", mart = ensembl)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
cytoplasmictranslation_killifish_GO <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Human[i] %in% cytoplasmictranslation_GO$hgnc_symbol) {
    cytoplasmictranslation_killifish_GO <- append(cytoplasmictranslation_killifish_GO, gene_names_df$`N. furzeri (NCBI)`[i])
  }
}
SR_WL_myeloid_20251202_cytoplasmictranslation_sum_expr <- rowSums(FetchData(object = SR_WL_myeloid_20251202, vars = cytoplasmictranslation_killifish_GO, layer = "data"))
SR_WL_myeloid_20251202$cytoplasmictranslation_sum_GO_expr <- SR_WL_myeloid_20251202_cytoplasmictranslation_sum_expr
VlnPlot(SR_WL_myeloid_20251202, features = "cytoplasmictranslation_sum_GO_expr", group.by = "orig.ident", cols = c("#b660cd", "grey"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = SR_WL_myeloid_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  ylim(c(-1, 400)) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_SRmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "SR")
cells_WLmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "WL")

# Extract the metadata column values
SRmyeloid_cytoplasmictranslation <- SR_WL_myeloid_20251202@meta.data[cells_SRmyeloid, "cytoplasmictranslation_sum_GO_expr"]
WLmyeloid_cytoplasmictranslation <- SR_WL_myeloid_20251202@meta.data[cells_WLmyeloid, "cytoplasmictranslation_sum_GO_expr"]

# Test
wilcox.test(SRmyeloid_cytoplasmictranslation, WLmyeloid_cytoplasmictranslation)

# Fig 4S1e
# Violin plot for GO term "endocytosis"
# SRmyeloid = purple, WLmyeloid = grey
# Save at 300W x 500H
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
endocytosis_GO <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "go_id"), filters = "go", values = "GO:0006897", mart = ensembl)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
endocytosis_killifish_GO <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Human[i] %in% endocytosis_GO$hgnc_symbol) {
    endocytosis_killifish_GO <- append(endocytosis_killifish_GO, gene_names_df$`N. furzeri (NCBI)`[i])
  }
}
SR_WL_myeloid_20251202_endocytosis_sum_expr <- rowSums(FetchData(object = SR_WL_myeloid_20251202, vars = endocytosis_killifish_GO, layer = "data"))
SR_WL_myeloid_20251202$endocytosis_sum_GO_expr <- SR_WL_myeloid_20251202_endocytosis_sum_expr
VlnPlot(SR_WL_myeloid_20251202, features = "endocytosis_sum_GO_expr", group.by = "orig.ident", cols = c("#b660cd", "grey"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = SR_WL_myeloid_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  ylim(c(-1, 100)) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_SRmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "SR")
cells_WLmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "WL")

# Extract the metadata column values
SRmyeloid_endocytosis <- SR_WL_myeloid_20251202@meta.data[cells_SRmyeloid, "endocytosis_sum_GO_expr"]
WLmyeloid_endocytosis <- SR_WL_myeloid_20251202@meta.data[cells_WLmyeloid, "endocytosis_sum_GO_expr"]

# Test
wilcox.test(SRmyeloid_endocytosis, WLmyeloid_endocytosis)

# Fig 4S1e
# Violin plot for GO term "positive regulation of endocytosis"
# SRmyeloid = purple, WLmyeloid = grey
# Save at 300W x 500H
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
PRendocytosis_GO <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "go_id"), filters = "go", values = "GO:0045807", mart = ensembl)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
PRendocytosis_killifish_GO <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Human[i] %in% PRendocytosis_GO$hgnc_symbol) {
    PRendocytosis_killifish_GO <- append(PRendocytosis_killifish_GO, gene_names_df$`N. furzeri (NCBI)`[i])
  }
}
SR_WL_myeloid_20251202_PRendocytosis_sum_expr <- rowSums(FetchData(object = SR_WL_myeloid_20251202, vars = PRendocytosis_killifish_GO, layer = "data"))
SR_WL_myeloid_20251202$PRendocytosis_sum_GO_expr <- SR_WL_myeloid_20251202_PRendocytosis_sum_expr
VlnPlot(SR_WL_myeloid_20251202, features = "PRendocytosis_sum_GO_expr", group.by = "orig.ident", cols = c("#b660cd", "grey"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = SR_WL_myeloid_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  ylim(c(-1, 25)) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_SRmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "SR")
cells_WLmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "WL")

# Extract the metadata column values
SRmyeloid_PRendocytosis <- SR_WL_myeloid_20251202@meta.data[cells_SRmyeloid, "PRendocytosis_sum_GO_expr"]
WLmyeloid_PRendocytosis <- SR_WL_myeloid_20251202@meta.data[cells_WLmyeloid, "PRendocytosis_sum_GO_expr"]

# Test
wilcox.test(SRmyeloid_PRendocytosis, WLmyeloid_PRendocytosis)

# Fig 4S1e
# Violin plot for GO term "lysosome"
# SRmyeloid = purple, WLmyeloid = grey
# Save at 300W x 500H
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
lysosome_GO <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "go_id"), filters = "go", values = "GO:0005764", mart = ensembl)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
lysosome_killifish_GO <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Human[i] %in% lysosome_GO$hgnc_symbol) {
    lysosome_killifish_GO <- append(lysosome_killifish_GO, gene_names_df$`N. furzeri (NCBI)`[i])
  }
}
SR_WL_myeloid_20251202_lysosome_sum_expr <- rowSums(FetchData(object = SR_WL_myeloid_20251202, vars = lysosome_killifish_GO, layer = "data"))
SR_WL_myeloid_20251202$lysosome_sum_GO_expr <- SR_WL_myeloid_20251202_lysosome_sum_expr
VlnPlot(SR_WL_myeloid_20251202, features = "lysosome_sum_GO_expr", group.by = "orig.ident", cols = c("#b660cd", "grey"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = SR_WL_myeloid_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  ylim(c(-1, 250)) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_SRmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "SR")
cells_WLmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "WL")

# Extract the metadata column values
SRmyeloid_lysosome <- SR_WL_myeloid_20251202@meta.data[cells_SRmyeloid, "lysosome_sum_GO_expr"]
WLmyeloid_lysosome <- SR_WL_myeloid_20251202@meta.data[cells_WLmyeloid, "lysosome_sum_GO_expr"]

# Test
wilcox.test(SRmyeloid_lysosome, WLmyeloid_lysosome)

# Fig 4S1e
# Violin plot for GO term "positive regulation of lysosome organization"
# SRmyeloid = purple, WLmyeloid = grey
# Save at 300W x 500H
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
PRlysosome_GO <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "go_id"), filters = "go", values = "GO:1905673", mart = ensembl)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
PRlysosome_killifish_GO <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Human[i] %in% PRlysosome_GO$hgnc_symbol) {
    PRlysosome_killifish_GO <- append(PRlysosome_killifish_GO, gene_names_df$`N. furzeri (NCBI)`[i])
  }
}
SR_WL_myeloid_20251202_PRlysosome_sum_expr <- rowSums(FetchData(object = SR_WL_myeloid_20251202, vars = PRlysosome_killifish_GO, layer = "data"))
SR_WL_myeloid_20251202$PRlysosome_sum_GO_expr <- SR_WL_myeloid_20251202_PRlysosome_sum_expr
VlnPlot(SR_WL_myeloid_20251202, features = "PRlysosome_sum_GO_expr", group.by = "orig.ident", cols = c("#b660cd", "grey"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = SR_WL_myeloid_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  ylim(c(-1, 8.5)) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_SRmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "SR")
cells_WLmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "WL")

# Extract the metadata column values
SRmyeloid_PRlysosome <- SR_WL_myeloid_20251202@meta.data[cells_SRmyeloid, "PRlysosome_sum_GO_expr"]
WLmyeloid_PRlysosome <- SR_WL_myeloid_20251202@meta.data[cells_WLmyeloid, "PRlysosome_sum_GO_expr"]

# Test
wilcox.test(SRmyeloid_PRlysosome, WLmyeloid_PRlysosome)

# Fig 4S1e
# Violin plot for GO term "vacuolar acidification"
# SRmyeloid = purple, WLmyeloid = grey
# Save at 300W x 500H
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
vacuolaracidification_GO <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "go_id"), filters = "go", values = "GO:0007035", mart = ensembl)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
vacuolaracidification_killifish_GO <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Human[i] %in% vacuolaracidification_GO$hgnc_symbol) {
    vacuolaracidification_killifish_GO <- append(vacuolaracidification_killifish_GO, gene_names_df$`N. furzeri (NCBI)`[i])
  }
}
SR_WL_myeloid_20251202_vacuolaracidification_sum_expr <- rowSums(FetchData(object = SR_WL_myeloid_20251202, vars = vacuolaracidification_killifish_GO, layer = "data"))
SR_WL_myeloid_20251202$vacuolaracidification_sum_GO_expr <- SR_WL_myeloid_20251202_vacuolaracidification_sum_expr
VlnPlot(SR_WL_myeloid_20251202, features = "vacuolaracidification_sum_GO_expr", group.by = "orig.ident", cols = c("#b660cd", "grey"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = SR_WL_myeloid_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  ylim(c(-1, 35)) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_SRmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "SR")
cells_WLmyeloid <- WhichCells(SR_WL_myeloid_20251202, expression = orig.ident == "WL")

# Extract the metadata column values
SRmyeloid_vacuolaracidification <- SR_WL_myeloid_20251202@meta.data[cells_SRmyeloid, "vacuolaracidification_sum_GO_expr"]
WLmyeloid_vacuolaracidification <- SR_WL_myeloid_20251202@meta.data[cells_WLmyeloid, "vacuolaracidification_sum_GO_expr"]

# Test
wilcox.test(SRmyeloid_vacuolaracidification, WLmyeloid_vacuolaracidification)

# Save supplemental tables

# DEG SR_WL_myeloid (all)
DEG_SR_WL_myeloid_20251202_all <- FindMarkers(
  object = SR_WL_myeloid_20251202,
  ident.1 = "SR",
  ident.2 = "WL",
  group.by = "orig.ident",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
)
DEG_SR_WL_myeloid_20251202_all_df <- as.data.frame(DEG_SR_WL_myeloid_20251202_all)
DEG_SR_WL_myeloid_20251202_all_df$gene <- rownames(DEG_SR_WL_myeloid_20251202_all_df)

# Add human, mouse, and zebrafish gene names for DEGs
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(DEG_SR_WL_myeloid_20251202_all_df)) {
  killifish_NCBI_gene_name = DEG_SR_WL_myeloid_20251202_all_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {DEG_SR_WL_myeloid_20251202_all_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {DEG_SR_WL_myeloid_20251202_all_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {DEG_SR_WL_myeloid_20251202_all_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {DEG_SR_WL_myeloid_20251202_all_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {DEG_SR_WL_myeloid_20251202_all_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {DEG_SR_WL_myeloid_20251202_all_df$zebrafish_gene_name[i] = ""}
}

# Save
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/")
write.table(DEG_SR_WL_myeloid_20251202_all_df, "SupplementaryTable12.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# GO:BP terms up in SRmyeloid
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/")
write.table(SRmyeloid_enriched_vsWLmyeloid_GOterms_BP_df, "SupplementaryTable13.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# GO:BP terms up in WLmyeloid
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/")
write.table(WLmyeloid_enriched_vsSRmyeloid_GOterms_BP_df, "SupplementaryTable14.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Save Seurat objects
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(SR_20251202, file = "SR_20251202.rds")
saveRDS(SR_WL_20251202, file = "SR_WL_20251202.rds")
saveRDS(SR_WL_myeloid_20251202, file = "SR_WL_myeloid_20251202.rds")