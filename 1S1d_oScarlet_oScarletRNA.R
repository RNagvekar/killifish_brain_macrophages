rm(list=ls())
library(ggplot2)
library(patchwork)
library(tidyverse)
library(dplyr)

# Set working directory
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/1S1e/1S1e_QuPath_proj/")

# Read .csv
measurements_oScarlet_oScarletRNA_df <- read.csv("1S1e_oScarlet_oScarletRNA_measurements.csv")

# Add columns with image data
measurements_oScarlet_oScarletRNA_df <- measurements_oScarlet_oScarletRNA_df %>%
  separate(col = Image, into = c("Date", "nImage", "Slide", "Animal", "Section", "Region", "Age", "Sex", "Channel405", "Channel488", "Channel546", "Channel647", "Magnification", "Notes"), sep = "_", remove = FALSE) %>%
  mutate(Notes = sub("\\.[^.]+$", "", Notes))

# Add a logical column for oScarletpositivebymean based on Cell..oScarlet.mean
measurements_oScarlet_oScarletRNA_df <- measurements_oScarlet_oScarletRNA_df %>%
  group_by(Image) %>%
  mutate(
    oScarletpositivebymean_threshold = mean(Cell..oScarlet.mean) + 2.5 * sd(Cell..oScarlet.mean),
    oScarletpositivebymean = Cell..oScarlet.mean > oScarletpositivebymean_threshold) %>%
  ungroup()

# Add a logical column for oScarletpositivebymax based on Cell..oScarlet.max
measurements_oScarlet_oScarletRNA_df <- measurements_oScarlet_oScarletRNA_df %>%
  mutate(
    oScarletpositivebymax_threshold = 25000,
    oScarletpositivebymax = Cell..oScarlet.max > oScarletpositivebymax_threshold)

# Add a logical column for oScarletpositive
measurements_oScarlet_oScarletRNA_df <- measurements_oScarlet_oScarletRNA_df %>%
  mutate(
    oScarletpositive = oScarletpositivebymean & oScarletpositivebymax)

# Add a logical column for oScarlet protein positive with 0 oScarlet RNA transcripts
measurements_oScarlet_oScarletRNA_df <- measurements_oScarlet_oScarletRNA_df %>%
  mutate(
    oScarletpositive_withoutoScarletRNA = oScarletpositive & (Subcellular..Channel.4..Num.spots.estimated == 0))

# Histograms
hist(subset(measurements_oScarlet_oScarletRNA_df, oScarletpositive_withoutoScarletRNA)$Cell..oScarlet.mean)
hist(subset(measurements_oScarlet_oScarletRNA_df, !oScarletpositive_withoutoScarletRNA)$Cell..oScarlet.mean)

bins <- 50
ymax <- max(
  max(hist(subset(measurements_oScarlet_oScarletRNA_df, oScarletpositive_withoutoScarletRNA)$Cell..oScarlet.mean, breaks = bins, plot = FALSE)$counts),
  max(hist(subset(measurements_oScarlet_oScarletRNA_df, !oScarletpositive_withoutoScarletRNA)$Cell..oScarlet.mean, breaks = bins, plot = FALSE)$counts)
)

hist_oScarletpositive_withoutoScarletRNA <- ggplot(subset(measurements_oScarlet_oScarletRNA_df, oScarletpositive_withoutoScarletRNA), aes(x = Cell..oScarlet.mean)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins = bins, fill = "red", color = NA) +
  coord_cartesian(ylim = c(0, 20)) +
  ylab("Percent of cells") +
  theme_classic()
hist_other <- ggplot(subset(measurements_oScarlet_oScarletRNA_df, !oScarletpositive_withoutoScarletRNA), aes(x = Cell..oScarlet.mean)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins = bins, fill = "grey", color = NA) +
  coord_cartesian(ylim = c(0, 20)) +
  ylab("Percent of cells") +
  theme_classic()

hist_oScarletpositive_withoutoScarletRNA / hist_other

# oScarlet mean vs. max
ggplot(measurements_oScarlet_oScarletRNA_df,
       aes(x = Cell..oScarlet.mean,
           y = Cell..oScarlet.max,
           color = oScarletpositive)) +
  geom_point(size = 1,
             position = position_jitter(width = 0.0, height = 0.0)) +
  scale_x_continuous(limits = c(0, 32000), labels = scales::comma) +
  scale_color_manual(values = c(
    "TRUE" = "red",
    "FALSE" = "grey")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

# Fig 1S1e - oScarlet vs. oScarlet RNA
# oScarlet RNA transcripts on y-axis
# Save at 1000W x 400H
ggplot(measurements_oScarlet_oScarletRNA_df,
       aes(x = Cell..oScarlet.mean,
           y = Subcellular..Channel.4..Num.spots.estimated,
           color = oScarletpositive_withoutoScarletRNA,
           alpha = 0.5)) +
  geom_point(size = 2,
             position = position_jitter(width = 0.0, height = 0.0)) +
  scale_x_continuous(limits = c(0, 32000), labels = scales::comma) +
  scale_color_manual(values = c(
    "TRUE" = "red",
    "FALSE" = "grey")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line.x = element_line(color = "black", size = 0.25),
        axis.line.y = element_line(color = "black", size = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        legend.position = "none")