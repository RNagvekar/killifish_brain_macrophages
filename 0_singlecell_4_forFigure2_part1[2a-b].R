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

# Fig 2a
# Dotplot of markers in YSNR (above) and YSR (below)
# Save at 700W x 400H
# LOC107379395 = apoeb
# LOC107372638 = HLA-DPB1 (MHCII)
# LOC107385068 = c1qc
# LOC107385067 = c1qa
# LOC107383768 = mrc1
markers <- c("csf1r",
             "LOC107379395",
             "LOC107372638",
             "cd74",
             "LOC107385068",
             "LOC107385067",
             "LOC107383768",
             "marco",
             "tmem119",
             "f13a1",
             "cldn5")
dotplot_YSNR_YSR <- DotPlot(object = YSR_YSNR_20251202,
                            features = markers,
                            group.by = "orig.ident",
                            scale = FALSE)
ggplot(dotplot_YSNR_YSR$data, aes(x = features.plot, y = id)) +
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

# Fig 2b
# Violin plot for GO term "endocytosis"
# YSNR = grey, YSR = red
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
YSR_YSNR_20251202_endocytosis_sum_expr <- rowSums(FetchData(object = YSR_YSNR_20251202, vars = endocytosis_killifish_GO, layer = "data"))
YSR_YSNR_20251202$endocytosis_sum_GO_expr <- YSR_YSNR_20251202_endocytosis_sum_expr
VlnPlot(YSR_YSNR_20251202, features = "endocytosis_sum_GO_expr", group.by = "orig.ident", cols = c("gray", "#ff7f7f"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_YSNR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSNR")
cells_YSR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSR")

# Extract the metadata column values
YSNR_endocytosis <- YSR_YSNR_20251202@meta.data[cells_YSNR, "endocytosis_sum_GO_expr"]
YSR_endocytosis <- YSR_YSNR_20251202@meta.data[cells_YSR, "endocytosis_sum_GO_expr"]

# Test
wilcox.test(YSNR_endocytosis, YSR_endocytosis)

# Fig 2b
# Violin plot for GO term "positive regulation of endocytosis"
# YSNR = grey, YSR = red
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
YSR_YSNR_20251202_PRendocytosis_sum_expr <- rowSums(FetchData(object = YSR_YSNR_20251202, vars = PRendocytosis_killifish_GO, layer = "data"))
YSR_YSNR_20251202$PRendocytosis_sum_GO_expr <- YSR_YSNR_20251202_PRendocytosis_sum_expr
VlnPlot(YSR_YSNR_20251202, features = "PRendocytosis_sum_GO_expr", group.by = "orig.ident", cols = c("gray", "#ff7f7f"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_YSNR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSNR")
cells_YSR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSR")

# Extract the metadata column values
YSNR_PRendocytosis <- YSR_YSNR_20251202@meta.data[cells_YSNR, "PRendocytosis_sum_GO_expr"]
YSR_PRendocytosis <- YSR_YSNR_20251202@meta.data[cells_YSR, "PRendocytosis_sum_GO_expr"]

# Test
wilcox.test(YSNR_PRendocytosis, YSR_PRendocytosis)

# Fig 2b
# Violin plot for GO term "lysosome"
# YSNR = grey, YSR = red
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
YSR_YSNR_20251202_lysosome_sum_expr <- rowSums(FetchData(object = YSR_YSNR_20251202, vars = lysosome_killifish_GO, layer = "data"))
YSR_YSNR_20251202$lysosome_sum_GO_expr <- YSR_YSNR_20251202_lysosome_sum_expr
VlnPlot(YSR_YSNR_20251202, features = "lysosome_sum_GO_expr", group.by = "orig.ident", cols = c("gray", "#ff7f7f"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_YSNR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSNR")
cells_YSR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSR")

# Extract the metadata column values
YSNR_lysosome <- YSR_YSNR_20251202@meta.data[cells_YSNR, "lysosome_sum_GO_expr"]
YSR_lysosome <- YSR_YSNR_20251202@meta.data[cells_YSR, "lysosome_sum_GO_expr"]

# Test
wilcox.test(YSNR_lysosome, YSR_lysosome)

# Fig 2b
# Violin plot for GO term "positive regulation of lysosome organization"
# YSNR = grey, YSR = red
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
YSR_YSNR_20251202_PRlysosome_sum_expr <- rowSums(FetchData(object = YSR_YSNR_20251202, vars = PRlysosome_killifish_GO, layer = "data"))
YSR_YSNR_20251202$PRlysosome_sum_GO_expr <- YSR_YSNR_20251202_PRlysosome_sum_expr
VlnPlot(YSR_YSNR_20251202, features = "PRlysosome_sum_GO_expr", group.by = "orig.ident", cols = c("gray", "#ff7f7f"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_YSNR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSNR")
cells_YSR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSR")

# Extract the metadata column values
YSNR_PRlysosome <- YSR_YSNR_20251202@meta.data[cells_YSNR, "PRlysosome_sum_GO_expr"]
YSR_PRlysosome <- YSR_YSNR_20251202@meta.data[cells_YSR, "PRlysosome_sum_GO_expr"]

# Test
wilcox.test(YSNR_PRlysosome, YSR_PRlysosome)

# Fig 2b
# Violin plot for GO term "vacuolar acidification"
# YSNR = grey, YSR = red
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
YSR_YSNR_20251202_vacuolaracidification_sum_expr <- rowSums(FetchData(object = YSR_YSNR_20251202, vars = vacuolaracidification_killifish_GO, layer = "data"))
YSR_YSNR_20251202$vacuolaracidification_sum_GO_expr <- YSR_YSNR_20251202_vacuolaracidification_sum_expr
VlnPlot(YSR_YSNR_20251202, features = "vacuolaracidification_sum_GO_expr", group.by = "orig.ident", cols = c("gray", "#ff7f7f"), pt.size = 0) +
  geom_point(size = 0.5, alpha = 0.4, color = "gray", position = position_jitter(w = 0.35, h = 0)) +
  geom_violin(aes(fill = YSR_YSNR_20251202$orig.ident),
              alpha = 0.8,
              trim = TRUE,
              scale = "width") +
  stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.line = element_line(color = "black", size = 0.5))

# Wilcoxon test
# Cells in each orig.ident group
cells_YSNR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSNR")
cells_YSR <- WhichCells(YSR_YSNR_20251202, expression = orig.ident == "YSR")

# Extract the metadata column values
YSNR_vacuolaracidification <- YSR_YSNR_20251202@meta.data[cells_YSNR, "vacuolaracidification_sum_GO_expr"]
YSR_vacuolaracidification <- YSR_YSNR_20251202@meta.data[cells_YSR, "vacuolaracidification_sum_GO_expr"]

# Test
wilcox.test(YSNR_vacuolaracidification, YSR_vacuolaracidification)

# Save Seurat objects
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_singlecell/20251202/")
saveRDS(YSR_YSNR_20251202, file = "YSR_YSNR_20251202.rds")