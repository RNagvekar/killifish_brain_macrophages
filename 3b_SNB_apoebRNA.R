rm(list=ls())
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(FNN)

# Set working directory
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/3b/3b_QuPath_proj/")

# Read .csv
measurements_SNB_apoebRNA_df <- read.csv("SNB_apoebRNA_measurements.csv")

# Add columns with image data
measurements_SNB_apoebRNA_df <- measurements_SNB_apoebRNA_df %>%
  separate(col = Image, into = c("Date", "nImage", "Slide", "Animal", "Section", "Region", "Age", "Sex", "Channel405", "Channel488", "Channel546", "Channel647", "Magnification", "Notes"), sep = "_", remove = FALSE) %>%
  mutate(Notes = sub("\\.[^.]+$", "", Notes))

# Add a logical column for SNBpositive based on Cell..SNB.mean
measurements_SNB_apoebRNA_df <- measurements_SNB_apoebRNA_df %>%
  group_by(Image) %>%
  mutate(
    SNBpositive_threshold = mean(Cell..SNB.mean) + 1.5 * sd(Cell..SNB.mean),
    SNBpositive = Cell..SNB.mean > SNBpositive_threshold) %>%
  ungroup()

# Add a logical column for apoebRNApositive based on Cell..apoebRNA.mean
measurements_SNB_apoebRNA_df <- measurements_SNB_apoebRNA_df %>%
  group_by(Image) %>%
  mutate(
    apoebRNApositive_threshold = mean(Cell..apoebRNA.mean) + 2.5 * sd(Cell..apoebRNA.mean),
    apoebRNApositive = Cell..apoebRNA.mean > apoebRNApositive_threshold) %>%
  ungroup()

# Add a column for Experiment
measurements_SNB_apoebRNA_df <- measurements_SNB_apoebRNA_df %>%
  mutate(Experiment = case_when(
    Animal %in% c("W121", "W128", "W133", "W151") ~ "Imaging4",
    Animal %in% c("W125", "W148") ~ "Imaging5",
    TRUE ~ NA_character_
  ))

# Distance to nearest SNB+ cell
measurements_SNB_apoebRNA_df <- measurements_SNB_apoebRNA_df %>% 
  group_by(Image) %>% 
  group_modify(~{
    A <- .x
    B <- A[A$SNBpositive == TRUE, ]
    
    if(nrow(B) == 0){
      A$dist_to_SNBpositive <- NA
      return(A)
    }
    
    # build matrices of centroids
    matA <- as.matrix(A[, c("Centroid.X.µm", "Centroid.Y.µm")])
    matB <- as.matrix(B[, c("Centroid.X.µm", "Centroid.Y.µm")])
    
    # k=1 nearest neighbor distances
    nn <- get.knnx(data = matB, query = matA, k = 1)
    
    A$dist_to_SNBpositive <- nn$nn.dist[, 1]
    A
  }) %>%
  ungroup()

# Percent SNBpositive by region, animal
# Not used for quantification, only supplementary table
pctSNBpositive_measurements_SNB_apoebRNA_df <- measurements_SNB_apoebRNA_df %>%
  group_by(Region, Animal, Experiment, Sex) %>%
  summarise(
    pct_SNBpositive = mean(SNBpositive) * 100,
  )

# Percent apoebRNApositive by region, animal
# Used for quantification
pctapoebRNApositive_measurements_SNB_apoebRNA_df <- measurements_SNB_apoebRNA_df %>%
  group_by(Region, Animal, Experiment, Sex) %>%
  summarise(
    pct_apoebRNApositive_within_1µm_of_SNBpositive = mean(apoebRNApositive[dist_to_SNBpositive < 1], na.rm = TRUE) * 100,
    pct_apoebRNApositive_within_5µm_of_SNBpositive = mean(apoebRNApositive[dist_to_SNBpositive < 5], na.rm = TRUE) * 100,
    pct_apoebRNApositive_within_10µm_of_SNBpositive = mean(apoebRNApositive[dist_to_SNBpositive < 10], na.rm = TRUE) * 100,
    pct_apoebRNApositive_ofallcells = mean(apoebRNApositive) * 100
  )

# Summarize by apoebRNA status, animal, and region
# Not used for quantification, only supplementary table
by_animal_region_and_apoebRNA_status_measurements_SNB_apoebRNA_df <- measurements_SNB_apoebRNA_df %>%
  group_by(apoebRNApositive, Animal, Region, Experiment, Sex) %>%
  summarise(
    n = n(),
    mean_dist_to_SNBpositive = mean(dist_to_SNBpositive),
    median_dist_to_SNBpositive = median(dist_to_SNBpositive),
    pct_within_1µm_of_SNBpositive = 100 * sum(dist_to_SNBpositive < 1) / n,
    pct_within_5µm_of_SNBpositive = 100 * sum(dist_to_SNBpositive < 5) / n,
    pct_within_10µm_of_SNBpositive = 100 * sum(dist_to_SNBpositive < 10) / n
  )

# Same-experiment all cells (sE_allcells) averages per region
sE_allcells_avg_pctapoebRNApositive_measurements_SNB_apoebRNA_df <- pctapoebRNApositive_measurements_SNB_apoebRNA_df %>%
  group_by(Experiment, Region) %>%
  summarise(
    sE_allcells_avg_pctapoebRNApositive = mean(pct_apoebRNApositive_ofallcells)
  )

# Create a norm df with normalized values per region
norm_pctapoebRNApositive_measurements_SNB_apoebRNA_df <- pctapoebRNApositive_measurements_SNB_apoebRNA_df %>%
  left_join(sE_allcells_avg_pctapoebRNApositive_measurements_SNB_apoebRNA_df,
            by = c("Region", "Experiment")) %>%
  mutate(across(starts_with("pct"),
                ~ .x / sE_allcells_avg_pctapoebRNApositive, 
                .names = "norm_{.col}")) %>%
  dplyr::select(Region, Animal, Experiment, Sex, starts_with("norm_"))

# Pivot to long
norm_forplot_df <- norm_pctapoebRNApositive_measurements_SNB_apoebRNA_df %>%
  pivot_longer(
    cols = c(norm_pct_apoebRNApositive_ofallcells, norm_pct_apoebRNApositive_within_5µm_of_SNBpositive),
    names_to = "xgroup",
    values_to = "yval"
  ) %>%
  dplyr::select(Region, Animal, Experiment, Sex, xgroup, yval)

# Summary statistics
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
summary_stats_norm_forplot_df <- norm_forplot_df %>%
  group_by(Region, xgroup) %>%
  summarise(across(c("yval"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))
summary_stats_norm_forplot_df <- summary_stats_norm_forplot_df %>%
  mutate(xgroup_factor = factor(xgroup))

# Mean colocalized (pct_within_5µm_of_SNBpositive) by region

# Fig 3b
# Forebrain (telencephalon; TEL)
# Save at 200W x 300H
ggplot(data = subset(norm_forplot_df, Region == "TEL"),
       aes(x = xgroup, y = yval)) +
  ylim(c(0, 4)) +
  geom_point(aes(shape = Sex), size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  scale_shape_manual(values = c("M" = 15, "F" = 17)) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_forplot_df, Region == "TEL"),
               aes(x = as.numeric(xgroup_factor) - 0.1, xend = as.numeric(xgroup_factor) + 0.1,
                   y = yval_mean, yend = yval_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_forplot_df, Region == "TEL"),
                aes(x = as.numeric(xgroup_factor), y = yval_mean,
                    ymin = yval_mean - yval_se,
                    ymax = yval_mean + yval_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5)+
  #ggtitle("norm pct apoebRNApositive") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line.x = element_line(color = "black", size = 0.25),
        axis.line.y = element_line(color = "black", size = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

# Wilcoxon test (paired)
norm_forplot_TEL_df <- subset(norm_forplot_df, Region == "TEL")
wilcox.test(subset(norm_forplot_TEL_df, xgroup == "norm_pct_apoebRNApositive_ofallcells")$yval,
            subset(norm_forplot_TEL_df, xgroup == "norm_pct_apoebRNApositive_within_5µm_of_SNBpositive")$yval,
            paired = TRUE)

# Fig 3b
# Midbrain (optic tectum; OT)
# Save at 200W x 300H
ggplot(data = subset(norm_forplot_df, Region == "OT"),
       aes(x = xgroup, y = yval)) +
  ylim(c(0, 4)) +
  geom_point(aes(shape = Sex), size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  scale_shape_manual(values = c("M" = 15, "F" = 17)) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_forplot_df, Region == "OT"),
               aes(x = as.numeric(xgroup_factor) - 0.1, xend = as.numeric(xgroup_factor) + 0.1,
                   y = yval_mean, yend = yval_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_forplot_df, Region == "OT"),
                aes(x = as.numeric(xgroup_factor), y = yval_mean,
                    ymin = yval_mean - yval_se,
                    ymax = yval_mean + yval_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5)+
  #ggtitle("norm pct apoebRNApositive") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line.x = element_line(color = "black", size = 0.25),
        axis.line.y = element_line(color = "black", size = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

# Wilcoxon test (paired)
norm_forplot_OT_df <- subset(norm_forplot_df, Region == "OT")
wilcox.test(subset(norm_forplot_OT_df, xgroup == "norm_pct_apoebRNApositive_ofallcells")$yval,
            subset(norm_forplot_OT_df, xgroup == "norm_pct_apoebRNApositive_within_5µm_of_SNBpositive")$yval,
            paired = TRUE)

# Fig 3b
# Hindbrain (HB)
# Save at 200W x 300H
ggplot(data = subset(norm_forplot_df, Region == "HB"),
       aes(x = xgroup, y = yval)) +
  ylim(c(0, 4)) +
  geom_point(aes(shape = Sex), size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  scale_shape_manual(values = c("M" = 15, "F" = 17)) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_forplot_df, Region == "HB"),
               aes(x = as.numeric(xgroup_factor) - 0.1, xend = as.numeric(xgroup_factor) + 0.1,
                   y = yval_mean, yend = yval_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_forplot_df, Region == "HB"),
                aes(x = as.numeric(xgroup_factor), y = yval_mean,
                    ymin = yval_mean - yval_se,
                    ymax = yval_mean + yval_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5)+
  #ggtitle("norm pct apoebRNApositive") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line.x = element_line(color = "black", size = 0.25),
        axis.line.y = element_line(color = "black", size = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

# Wilcoxon test (paired)
norm_forplot_HB_df <- subset(norm_forplot_df, Region == "HB")
wilcox.test(subset(norm_forplot_HB_df, xgroup == "norm_pct_apoebRNApositive_ofallcells")$yval,
            subset(norm_forplot_HB_df, xgroup == "norm_pct_apoebRNApositive_within_5µm_of_SNBpositive")$yval,
            paired = TRUE)

# Spatial maps
ggplot(data = subset(measurements_SNB_apoebRNA_df,
                     (Date == "20251125" & nImage == "06")),
       aes(x = Centroid.X.µm,
           y = -Centroid.Y.µm,
           color = SNBpositive)) +
  geom_point(size = 1,
             position = position_jitter(width = 0.0, height = 0.0)) +
  theme_minimal()

# Save source data tables
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/3b/")
write.csv(pctSNBpositive_measurements_SNB_apoebRNA_df, "SourceDataTable4.csv", fileEncoding = "latin1")
write.csv(by_animal_region_and_apoebRNA_status_measurements_SNB_apoebRNA_df, "SourceDataTable5.csv", fileEncoding = "latin1")
write.csv(pctapoebRNApositive_measurements_SNB_apoebRNA_df, "SourceDataTable6.csv", fileEncoding = "latin1")