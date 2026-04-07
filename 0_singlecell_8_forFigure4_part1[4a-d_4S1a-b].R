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

# Read Seurat object (YSR_OSR)
YSR_OSR_20251202 <- readRDS(file = "YSR_OSR_20251202.rds")

# Extract UMAP coordinates from the Seurat object
YSR_OSR_umap_coords <- Embeddings(YSR_OSR_20251202, "umap")
YSR_OSR_xlims <- range(YSR_OSR_umap_coords[,1])
YSR_OSR_ylims <- range(YSR_OSR_umap_coords[,2])

# Reverse factor levels if needed (for young before old - check first otherwise you will un-reverse!)
# YSR_OSR_20251202$orig.ident <- factor(YSR_OSR_20251202$orig.ident, levels = rev(levels(factor(YSR_OSR_20251202$orig.ident))))

# Fig 4a
# UMAP by orig.ident
# OSR = purple, YSR = red
# Save at 500W x 500H
DimPlot(YSR_OSR_20251202, reduction = "umap", group.by = "orig.ident", cols = c("red", "#9b0fa5")) +
  coord_cartesian(xlim = YSR_OSR_xlims, ylim = YSR_OSR_ylims) +
  ggtitle(NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none")

# Fig 4b
# Dotplot of marker expression in OSR (above) and YSR (below)
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
dotplot_OSR_YSR <- DotPlot(object = YSR_OSR_20251202,
                            features = markers,
                            group.by = "orig.ident",
                            scale = FALSE)
ggplot(dotplot_OSR_YSR$data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(
    low = "yellow",
    high = "skyblue",
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

# Fig 4c
# Violin plot for GO term "cytoplasmic translation"
# YSR = red, OSR = purple
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
YSR_OSR_20251202_cytoplasmictranslation_sum_expr <- rowSums(FetchData(object = YSR_OSR_20251202, vars = cytoplasmictranslation_killifish_GO, layer = "data"))
YSR_OSR_20251202$cytoplasmictranslation_sum_GO_expr <- YSR_OSR_20251202_cytoplasmictranslation_sum_expr
VlnPlot(YSR_OSR_20251202, features = "cytoplasmictranslation_sum_GO_expr", group.by = "orig.ident", cols = c("#ff7f7f", "#e188ba"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_OSR_20251202$orig.ident),
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
cells_YSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "YSR")
cells_OSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "OSR")

# Extract the metadata column values
YSR_cytoplasmictranslation <- YSR_OSR_20251202@meta.data[cells_YSR, "cytoplasmictranslation_sum_GO_expr"]
OSR_cytoplasmictranslation <- YSR_OSR_20251202@meta.data[cells_OSR, "cytoplasmictranslation_sum_GO_expr"]

# Test
wilcox.test(YSR_cytoplasmictranslation, OSR_cytoplasmictranslation)

# Fig 4c
# Violin plot for GO term "vacuolar acidification"
# YSR = red, OSR = purple
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
YSR_OSR_20251202_vacuolaracidification_sum_expr <- rowSums(FetchData(object = YSR_OSR_20251202, vars = vacuolaracidification_killifish_GO, layer = "data"))
YSR_OSR_20251202$vacuolaracidification_sum_GO_expr <- YSR_OSR_20251202_vacuolaracidification_sum_expr
VlnPlot(YSR_OSR_20251202, features = "vacuolaracidification_sum_GO_expr", group.by = "orig.ident", cols = c("#ff7f7f", "#e188ba"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_OSR_20251202$orig.ident),
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
cells_YSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "YSR")
cells_OSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "OSR")

# Extract the metadata column values
YSR_vacuolaracidification <- YSR_OSR_20251202@meta.data[cells_YSR, "vacuolaracidification_sum_GO_expr"]
OSR_vacuolaracidification <- YSR_OSR_20251202@meta.data[cells_OSR, "vacuolaracidification_sum_GO_expr"]

# Test
wilcox.test(YSR_vacuolaracidification, OSR_vacuolaracidification)

# DEG between OSR and YSR
DEG_YSR_OSR_20251202 <- FindMarkers(
  object = YSR_OSR_20251202,
  ident.1 = "OSR",
  ident.2 = "YSR",
  group.by = "orig.ident",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.5
)
DEG_YSR_OSR_20251202_df <- as.data.frame(DEG_YSR_OSR_20251202)
DEG_YSR_OSR_20251202_df$gene <- rownames(DEG_YSR_OSR_20251202_df)

# Define background (human genes)
background_YSR_OSR_killifish <- rownames(YSR_OSR_20251202)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
background_YSR_OSR_human <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$`N. furzeri (NCBI)`[i] %in% background_YSR_OSR_killifish) {background_YSR_OSR_human <- append(background_YSR_OSR_human, gene_names_df$Human[i])}
}

# Add human, mouse, and zebrafish gene names for DEGs
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(DEG_YSR_OSR_20251202_df)) {
  killifish_NCBI_gene_name = DEG_YSR_OSR_20251202_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {DEG_YSR_OSR_20251202_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {DEG_YSR_OSR_20251202_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {DEG_YSR_OSR_20251202_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {DEG_YSR_OSR_20251202_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {DEG_YSR_OSR_20251202_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {DEG_YSR_OSR_20251202_df$zebrafish_gene_name[i] = ""}
}

# Enriched in each sample
OSR_enriched_vsYSR_DEG_df <- DEG_YSR_OSR_20251202_df %>% dplyr::filter(avg_log2FC > 0.5) %>% dplyr::filter(pct.1 > 0.1) %>% dplyr::filter(p_val_adj < 0.05)
YSR_enriched_vsOSR_DEG_df <- DEG_YSR_OSR_20251202_df %>% dplyr::filter(avg_log2FC < -0.5) %>% dplyr::filter(pct.2 > 0.1) %>% dplyr::filter(p_val_adj < 0.05)

# GO terms on genes enriched in OSR (relative to YSR)
OSR_enriched_vsYSR_GOterms_BP <- enrichr(OSR_enriched_vsYSR_DEG_df$human_gene_name, "GO_Biological_Process_2025", background = background_YSR_OSR_human)
OSR_enriched_vsYSR_GOterms_CC <- enrichr(OSR_enriched_vsYSR_DEG_df$human_gene_name, "GO_Cellular_Component_2025", background = background_YSR_OSR_human)
OSR_enriched_vsYSR_GOterms_MF <- enrichr(OSR_enriched_vsYSR_DEG_df$human_gene_name, "GO_Molecular_Function_2025", background = background_YSR_OSR_human)

OSR_enriched_vsYSR_GOterms_BP_df <- as.data.frame(OSR_enriched_vsYSR_GOterms_BP)
OSR_enriched_vsYSR_GOterms_BP_df$Term <- gsub(" \\(GO:[0-9]+\\)", "", OSR_enriched_vsYSR_GOterms_BP_df$GO_Biological_Process_2025.Term)
OSR_enriched_vsYSR_GOterms_BP_df$Term <- sub("^([a-z])", "\\U\\1", OSR_enriched_vsYSR_GOterms_BP_df$Term, perl = TRUE)
OSR_enriched_vsYSR_GOterms_BP_df$n_genes <- lengths(strsplit(as.character(OSR_enriched_vsYSR_GOterms_BP_df$GO_Biological_Process_2025.Genes), ";"))

# GO terms on genes enriched in YSR (relative to OSR)
YSR_enriched_vsOSR_GOterms_BP <- enrichr(YSR_enriched_vsOSR_DEG_df$human_gene_name, "GO_Biological_Process_2025", background = background_YSR_OSR_human)
YSR_enriched_vsOSR_GOterms_CC <- enrichr(YSR_enriched_vsOSR_DEG_df$human_gene_name, "GO_Cellular_Component_2025", background = background_YSR_OSR_human)
YSR_enriched_vsOSR_GOterms_MF <- enrichr(YSR_enriched_vsOSR_DEG_df$human_gene_name, "GO_Molecular_Function_2025", background = background_YSR_OSR_human)

YSR_enriched_vsOSR_GOterms_BP_df <- as.data.frame(YSR_enriched_vsOSR_GOterms_BP)
YSR_enriched_vsOSR_GOterms_BP_df$Term <- gsub(" \\(GO:[0-9]+\\)", "", YSR_enriched_vsOSR_GOterms_BP_df$GO_Biological_Process_2025.Term)
YSR_enriched_vsOSR_GOterms_BP_df$Term <- sub("^([a-z])", "\\U\\1", YSR_enriched_vsOSR_GOterms_BP_df$Term, perl = TRUE)
YSR_enriched_vsOSR_GOterms_BP_df$n_genes <- lengths(strsplit(as.character(YSR_enriched_vsOSR_GOterms_BP_df$GO_Biological_Process_2025.Genes), ";"))

# Fig 4d
# Plot for GO terms on genes enriched in YSR (relative to OSR)
# Save at 800W x 500H
ggplot(data = subset(YSR_enriched_vsOSR_GOterms_BP_df, GO_Biological_Process_2025.Adjusted.P.value < 0.05),
       aes(x = -log10(GO_Biological_Process_2025.Adjusted.P.value),
           y = reorder(Term, -log10(GO_Biological_Process_2025.Adjusted.P.value)),
           size = n_genes)) +
  geom_line(linetype = "blank") +
  geom_point(color = "red") +
  scale_size(range = c(6)) +
  xlim(0, 3) +
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

# Fig 4d
# Plot for GO terms on genes enriched in OSR (relative to YSR)
# Save at 1000W x 500H
ggplot(data = subset(OSR_enriched_vsYSR_GOterms_BP_df, GO_Biological_Process_2025.Adjusted.P.value < 1E-8),
       aes(x = -log10(GO_Biological_Process_2025.Adjusted.P.value),
           y = reorder(Term, -log10(GO_Biological_Process_2025.Adjusted.P.value)),
           size = n_genes)) +
  geom_line(linetype = "blank") +
  geom_point(color = "#9b0fa5") +
  scale_size(range = c(17, 23)) +
  xlim(0, 15) +
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

# Fig 4S1a
# Violin plot for GO term "endocytosis"
# YSR = red, OSR = purple
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
YSR_OSR_20251202_endocytosis_sum_expr <- rowSums(FetchData(object = YSR_OSR_20251202, vars = endocytosis_killifish_GO, layer = "data"))
YSR_OSR_20251202$endocytosis_sum_GO_expr <- YSR_OSR_20251202_endocytosis_sum_expr
VlnPlot(YSR_OSR_20251202, features = "endocytosis_sum_GO_expr", group.by = "orig.ident", cols = c("#ff7f7f", "#e188ba"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_OSR_20251202$orig.ident),
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
cells_YSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "YSR")
cells_OSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "OSR")

# Extract the metadata column values
YSR_endocytosis <- YSR_OSR_20251202@meta.data[cells_YSR, "endocytosis_sum_GO_expr"]
OSR_endocytosis <- YSR_OSR_20251202@meta.data[cells_OSR, "endocytosis_sum_GO_expr"]

# Test
wilcox.test(YSR_endocytosis, OSR_endocytosis)

# Fig 4S1a
# Violin plot for GO term "positive regulation of endocytosis"
# YSR = red, OSR = purple
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
YSR_OSR_20251202_PRendocytosis_sum_expr <- rowSums(FetchData(object = YSR_OSR_20251202, vars = PRendocytosis_killifish_GO, layer = "data"))
YSR_OSR_20251202$PRendocytosis_sum_GO_expr <- YSR_OSR_20251202_PRendocytosis_sum_expr
VlnPlot(YSR_OSR_20251202, features = "PRendocytosis_sum_GO_expr", group.by = "orig.ident", cols = c("#ff7f7f", "#e188ba"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_OSR_20251202$orig.ident),
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
cells_YSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "YSR")
cells_OSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "OSR")

# Extract the metadata column values
YSR_PRendocytosis <- YSR_OSR_20251202@meta.data[cells_YSR, "PRendocytosis_sum_GO_expr"]
OSR_PRendocytosis <- YSR_OSR_20251202@meta.data[cells_OSR, "PRendocytosis_sum_GO_expr"]

# Test
wilcox.test(YSR_PRendocytosis, OSR_PRendocytosis)

# Fig 4S1a
# Violin plot for GO term "lysosome"
# YSR = red, OSR = purple
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
YSR_OSR_20251202_lysosome_sum_expr <- rowSums(FetchData(object = YSR_OSR_20251202, vars = lysosome_killifish_GO, layer = "data"))
YSR_OSR_20251202$lysosome_sum_GO_expr <- YSR_OSR_20251202_lysosome_sum_expr
VlnPlot(YSR_OSR_20251202, features = "lysosome_sum_GO_expr", group.by = "orig.ident", cols = c("#ff7f7f", "#e188ba"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_OSR_20251202$orig.ident),
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
cells_YSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "YSR")
cells_OSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "OSR")

# Extract the metadata column values
YSR_lysosome <- YSR_OSR_20251202@meta.data[cells_YSR, "lysosome_sum_GO_expr"]
OSR_lysosome <- YSR_OSR_20251202@meta.data[cells_OSR, "lysosome_sum_GO_expr"]

# Test
wilcox.test(YSR_lysosome, OSR_lysosome)

# Fig 4S1a
# Violin plot for GO term "positive regulation of lysosome organization"
# YSR = red, OSR = purple
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
YSR_OSR_20251202_PRlysosome_sum_expr <- rowSums(FetchData(object = YSR_OSR_20251202, vars = PRlysosome_killifish_GO, layer = "data"))
YSR_OSR_20251202$PRlysosome_sum_GO_expr <- YSR_OSR_20251202_PRlysosome_sum_expr
VlnPlot(YSR_OSR_20251202, features = "PRlysosome_sum_GO_expr", group.by = "orig.ident", cols = c("#ff7f7f", "#e188ba"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_OSR_20251202$orig.ident),
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
cells_YSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "YSR")
cells_OSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "OSR")

# Extract the metadata column values
YSR_PRlysosome <- YSR_OSR_20251202@meta.data[cells_YSR, "PRlysosome_sum_GO_expr"]
OSR_PRlysosome <- YSR_OSR_20251202@meta.data[cells_OSR, "PRlysosome_sum_GO_expr"]

# Test
wilcox.test(YSR_PRlysosome, OSR_PRlysosome)

# Fig 4S1b
# Violin plot for GO term "inflammatory response"
# YSR = red, OSR = purple
# Save at 300W x 500H
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
inflammation_GO <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "go_id"), filters = "go", values = "GO:0006954", mart = ensembl)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
inflammation_killifish_GO <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Human[i] %in% inflammation_GO$hgnc_symbol) {
    inflammation_killifish_GO <- append(inflammation_killifish_GO, gene_names_df$`N. furzeri (NCBI)`[i])
  }
}
YSR_OSR_20251202_inflammation_sum_expr <- rowSums(FetchData(object = YSR_OSR_20251202, vars = inflammation_killifish_GO, layer = "data"))
YSR_OSR_20251202$inflammation_sum_GO_expr <- YSR_OSR_20251202_inflammation_sum_expr
VlnPlot(YSR_OSR_20251202, features = "inflammation_sum_GO_expr", group.by = "orig.ident", cols = c("#ff7f7f", "#e188ba"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_OSR_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  ylim(c(-1, 150)) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_YSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "YSR")
cells_OSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "OSR")

# Extract the metadata column values
YSR_inflammation <- YSR_OSR_20251202@meta.data[cells_YSR, "inflammation_sum_GO_expr"]
OSR_inflammation <- YSR_OSR_20251202@meta.data[cells_OSR, "inflammation_sum_GO_expr"]

# Test
wilcox.test(YSR_inflammation, OSR_inflammation)

# Fig 4S1b
# Violin plot for GO term "positive regulation of inflammatory response"
# YSR = red, OSR = purple
# Save at 300W x 500H
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
PRinflammation_GO <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "go_id"), filters = "go", values = "GO:0050729", mart = ensembl)
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
PRinflammation_killifish_GO <- c()
for (i in 1:nrow(gene_names_df)) {
  if (gene_names_df$Human[i] %in% PRinflammation_GO$hgnc_symbol) {
    PRinflammation_killifish_GO <- append(PRinflammation_killifish_GO, gene_names_df$`N. furzeri (NCBI)`[i])
  }
}
YSR_OSR_20251202_PRinflammation_sum_expr <- rowSums(FetchData(object = YSR_OSR_20251202, vars = PRinflammation_killifish_GO, layer = "data"))
YSR_OSR_20251202$PRinflammation_sum_GO_expr <- YSR_OSR_20251202_PRinflammation_sum_expr
VlnPlot(YSR_OSR_20251202, features = "PRinflammation_sum_GO_expr", group.by = "orig.ident", cols = c("#ff7f7f", "#e188ba"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_OSR_20251202$orig.ident),
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
cells_YSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "YSR")
cells_OSR <- WhichCells(YSR_OSR_20251202, expression = orig.ident == "OSR")

# Extract the metadata column values
YSR_PRinflammation <- YSR_OSR_20251202@meta.data[cells_YSR, "PRinflammation_sum_GO_expr"]
OSR_PRinflammation <- YSR_OSR_20251202@meta.data[cells_OSR, "PRinflammation_sum_GO_expr"]

# Test
wilcox.test(YSR_PRinflammation, OSR_PRinflammation)

# Save supplemental tables

# DEG YSR_OSR (all)
DEG_YSR_OSR_20251202_all <- FindMarkers(
  object = YSR_OSR_20251202,
  ident.1 = "OSR",
  ident.2 = "YSR",
  group.by = "orig.ident",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
)
DEG_YSR_OSR_20251202_all_df <- as.data.frame(DEG_YSR_OSR_20251202_all)
DEG_YSR_OSR_20251202_all_df$gene <- rownames(DEG_YSR_OSR_20251202_all_df)

# Add human, mouse, and zebrafish gene names for DEGs
gene_names_df <- read_excel("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/Symbols_processed_non-redundant_20170309.xlsx")
for (i in 1:nrow(DEG_YSR_OSR_20251202_all_df)) {
  killifish_NCBI_gene_name = DEG_YSR_OSR_20251202_all_df$gene[i]
  gene_name_row = which(gene_names_df[, "N. furzeri (NCBI)"] == killifish_NCBI_gene_name)
  if(length(gene_name_row) > 0) {DEG_YSR_OSR_20251202_all_df$human_gene_name[i] = gene_names_df$Human[gene_name_row]} else {DEG_YSR_OSR_20251202_all_df$human_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {DEG_YSR_OSR_20251202_all_df$mouse_gene_name[i] = gene_names_df$Mouse[gene_name_row]} else {DEG_YSR_OSR_20251202_all_df$mouse_gene_name[i] = ""}
  if(length(gene_name_row) > 0) {DEG_YSR_OSR_20251202_all_df$zebrafish_gene_name[i] = gene_names_df$Zebrafish[gene_name_row]} else {DEG_YSR_OSR_20251202_all_df$zebrafish_gene_name[i] = ""}
}

# Save
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/")
write.table(DEG_YSR_OSR_20251202_all_df, "SupplementaryTable9.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# GO:BP terms up in YSR
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/")
write.table(YSR_enriched_vsOSR_GOterms_BP_df, "SupplementaryTable10.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# GO:BP terms up in OSR
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/")
write.table(OSR_enriched_vsYSR_GOterms_BP_df, "SupplementaryTable11.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Save Seurat objects
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(OSR_20251202, file = "OSR_20251202.rds")
saveRDS(YSR_OSR_20251202, file = "YSR_OSR_20251202.rds")