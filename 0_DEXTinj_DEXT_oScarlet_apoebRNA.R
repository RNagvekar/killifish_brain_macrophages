rm(list=ls())
library(ggplot2)
library(tidyverse)
library(dplyr)

# Set working directory
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_DEXTinj/0_DEXTinj_QuPath_proj/")

# Read .csv
measurements_DEXT_oScarlet_apoebRNA_df <- read.csv("DEXT_oScarlet_apoebRNA_measurements.csv")

# Add columns with image data
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  separate(col = Image, into = c("Date", "nImage", "Slide", "Animal", "Section", "Region", "Age", "Sex", "Channel405", "Channel488", "Channel546", "Channel647", "Magnification", "Notes"), sep = "_", remove = FALSE) %>%
  mutate(Notes = sub("\\.[^.]+$", "", Notes))

# Add notes for images potentially confounded by areas with abnormally high density of apoebRNA signal (possible injury site?), as determined by visual inspection
# Image 20251123_15... also has a region with a high density of apoebRNA signal, but this is clearly associated with a blood vessel rather than a potential injury site
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  mutate(Notes = case_when(
    (Date == "20251122" & nImage == "11") |
      (Date == "20251122" & nImage == "18") |
      (Date == "20251123" & nImage == "08") |
      (Date == "20251123" & nImage == "11") ~ "apoebRNAhighdensity",
    TRUE ~ Notes
  ))

# Add a logical column for DEXTpositivebymean based on Cell..DEXT.mean
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Image) %>%
  mutate(
    DEXTpositivebymean_threshold = mean(Cell..DEXT.mean) + 2.5 * sd(Cell..DEXT.mean),
    DEXTpositivebymean = Cell..DEXT.mean > DEXTpositivebymean_threshold) %>%
  ungroup()

# Add a logical column for DEXTpositivebymax based on Cell..DEXT.max
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  mutate(
    DEXTpositivebymax_threshold = 25000,
    DEXTpositivebymax = Cell..DEXT.max > DEXTpositivebymax_threshold)

# Add a logical column for DEXTpositive
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  mutate(
    DEXTpositive = DEXTpositivebymean & DEXTpositivebymax)

# Add a logical column for oScarletpositivebymean based on Cell..oScarlet.mean
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Image) %>%
  mutate(
    oScarletpositivebymean_threshold = mean(Cell..oScarlet.mean) + 2.5 * sd(Cell..oScarlet.mean),
    oScarletpositivebymean = Cell..oScarlet.mean > oScarletpositivebymean_threshold) %>%
  ungroup()

# Add a logical column for oScarletpositivebymax based on Cell..oScarlet.max
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  mutate(
    oScarletpositivebymax_threshold = 25000,
    oScarletpositivebymax = Cell..oScarlet.max > oScarletpositivebymax_threshold)

# Add a logical column for oScarletpositive
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  mutate(
    oScarletpositive = oScarletpositivebymean & oScarletpositivebymax)

# Add a logical column for apoebRNApositive based on Cell..apoebRNA.mean
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Image) %>%
  mutate(
    apoebRNApositive_threshold = mean(Cell..apoebRNA.mean) + 2.5 * sd(Cell..apoebRNA.mean),
    apoebRNApositive = Cell..apoebRNA.mean > apoebRNApositive_threshold) %>%
  ungroup()

# Add a column for Experiment
measurements_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  mutate(Experiment = case_when(
    Animal %in% c("S1-DEXT", "S2-DEXT", "S3-DEXT", "S4-DEXT") ~ "Imaging2",
    Animal %in% c("S10-PBS-DEXT", "S11-PBS-DEXT", "S12-PBS-DEXT") ~ "Imaging3",
     TRUE ~ NA_character_
  ))

# Percent DEXTpositive, oScarletpositive, apoebRNApositive by image
# Note the "apoebRNAhighdensity" images are not outliers in apoebRNA+ in this analysis
# So we will do quantification with these images and also excluding them
pctpositivecells_DEXT_oScarlet_apoebRNA_byimage_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Image, Animal, Experiment) %>%
    summarise(
      pct_DEXTpositive = mean(DEXTpositive) * 100,
      pct_oScarletpositive = mean(oScarletpositive) * 100,
      pct_apoebRNApositive = mean(apoebRNApositive) * 100
    )

# Percent DEXTpositive, oScarletpositive, apoebRNApositive by animal
pctpositivecells_DEXT_oScarlet_apoebRNA_df <- measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Animal, Experiment) %>%
  summarise(
    pct_DEXTpositive = mean(DEXTpositive) * 100,
    pct_oScarletpositive = mean(oScarletpositive) * 100,
    pct_apoebRNApositive = mean(apoebRNApositive) * 100
  )

# Same-experiment (sE) averages
sE_avg_pctpositivecells_DEXT_oScarlet_apoebRNA_df <- pctpositivecells_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Experiment) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(.x, na.rm = TRUE), 
                   .names = "sE_avg_{.col}"))

# Create a norm df with normalized values 
norm_pctpositivecells_DEXT_oScarlet_apoebRNA_df <- pctpositivecells_DEXT_oScarlet_apoebRNA_df %>%
  left_join(sE_avg_pctpositivecells_DEXT_oScarlet_apoebRNA_df, by = c("Experiment")) %>%
  mutate(across(starts_with("pct"), 
                ~ . / get(paste0("sE_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(Animal, Experiment, starts_with("norm_"))

# Summary statistics
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
summary_stats_norm_pctpositivecells_DEXT_oScarlet_apoebRNA_df <- norm_pctpositivecells_DEXT_oScarlet_apoebRNA_df %>%
  ungroup() %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))

# Fig 3S1b
# Percent DEXT+
# Save at 200W x 300H
ggplot(data = norm_pctpositivecells_DEXT_oScarlet_apoebRNA_df, aes(x = 1, y = norm_pct_DEXTpositive)) +
  ylim(c(0, 1.5)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0),
             color = ifelse(norm_pctpositivecells_DEXT_oScarlet_apoebRNA_df$Animal == "S4-DEXT",
                            "black", "#0b6623"), alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = summary_stats_norm_pctpositivecells_DEXT_oScarlet_apoebRNA_df,
               aes(x = 0.9, xend = 1.1,
                   y = norm_pct_DEXTpositive_mean, yend = norm_pct_DEXTpositive_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = summary_stats_norm_pctpositivecells_DEXT_oScarlet_apoebRNA_df,
                aes(x = 1, y = norm_pct_DEXTpositive_mean,
                    ymin = norm_pct_DEXTpositive_mean - norm_pct_DEXTpositive_se,
                    ymax = norm_pct_DEXTpositive_mean + norm_pct_DEXTpositive_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  scale_x_continuous(
    breaks = c(1),
    labels = c("All"),
    limits = c(0.5, 1.5),
    expand = c(0,0)) +
  #ggtitle("norm pct DEXT") +
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

# Exclude S4-DEXT (low outlier in pctDEXTpositive, visually confirmed)
excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- subset(measurements_DEXT_oScarlet_apoebRNA_df, Animal != "S4-DEXT")

# Summarize by animal and DEXT, oScarlet, and apoebRNA status

# (1) Summarize by animal and DEXT status
by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Animal, DEXTpositive) %>%
  summarise(
    mean_DEXT = mean(Cell..DEXT.mean, na.rm = TRUE),
    mean_oScarlet = mean(Cell..oScarlet.mean, na.rm = TRUE),
    mean_apoebRNA = mean(Cell..apoebRNA.mean, na.rm = TRUE),
    Experiment = first(Experiment)
  )

# Same-experiment DEXTnegative (sE_DEXTnegative) averages
sE_DEXTnegative_avg_by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Experiment) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(.x[DEXTpositive == FALSE], na.rm = TRUE), 
                   .names = "sE_DEXTnegative_avg_{.col}"))

# Create a norm df with normalized values 
norm_by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  left_join(sE_DEXTnegative_avg_by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, by = c("Experiment")) %>%
  mutate(across(starts_with("mean"), 
                ~ . / get(paste0("sE_DEXTnegative_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(DEXTpositive, Animal, Experiment, starts_with("norm_"))

# Summary statistics
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
summary_stats_norm_by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- norm_by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(DEXTpositive) %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))

# (2) Summarize by animal and oScarlet status
by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Animal, oScarletpositive) %>%
  summarise(
    mean_DEXT = mean(Cell..DEXT.mean, na.rm = TRUE),
    mean_oScarlet = mean(Cell..oScarlet.mean, na.rm = TRUE),
    mean_apoebRNA = mean(Cell..apoebRNA.mean, na.rm = TRUE),
    Experiment = first(Experiment)
  )

# Same-experiment oScarletnegative (sE_oScarletnegative) averages
sE_oScarletnegative_avg_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Experiment) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(.x[oScarletpositive == FALSE], na.rm = TRUE), 
                   .names = "sE_oScarletnegative_avg_{.col}"))

# Create a norm df with normalized values 
norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  left_join(sE_oScarletnegative_avg_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, by = c("Experiment")) %>%
  mutate(across(starts_with("mean"), 
                ~ . / get(paste0("sE_oScarletnegative_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(oScarletpositive, Animal, Experiment, starts_with("norm_"))

# Summary statistics
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
summary_stats_norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(oScarletpositive) %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))

# (3) Summarize by animal and apoebRNA status
by_animal_and_apoebRNA_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Animal, apoebRNApositive) %>%
  summarise(
    mean_DEXT = mean(Cell..DEXT.mean, na.rm = TRUE),
    mean_oScarlet = mean(Cell..oScarlet.mean, na.rm = TRUE),
    mean_apoebRNA = mean(Cell..apoebRNA.mean, na.rm = TRUE),
    Experiment = first(Experiment)
  )

# Same-experiment apoebRNAnegative (sE_apoebRNAnegative) averages
sE_apoebRNAnegative_avg_by_animal_and_apoebRNA_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_apoebRNA_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Experiment) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(.x[apoebRNApositive == FALSE], na.rm = TRUE), 
                   .names = "sE_apoebRNAnegative_avg_{.col}"))

# Create a norm df with normalized values 
norm_by_animal_and_apoebRNA_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_apoebRNA_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  left_join(sE_apoebRNAnegative_avg_by_animal_and_apoebRNA_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, by = c("Experiment")) %>%
  mutate(across(starts_with("mean"), 
                ~ . / get(paste0("sE_apoebRNAnegative_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(apoebRNApositive, Animal, Experiment, starts_with("norm_"))

# Summary statistics
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
summary_stats_norm_by_animal_and_apoebRNA_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- norm_by_animal_and_apoebRNA_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(apoebRNApositive) %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))

# Mean DEXT by oScarlet status
# Save at 200W x 300H
ggplot(data = norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, aes(x = oScarletpositive, y = norm_mean_DEXT)) +
  ylim(c(0, 7)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#0b6623", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = summary_stats_norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
               aes(x = as.numeric(oScarletpositive) - 0.1, xend = as.numeric(oScarletpositive) + 0.1,
                   y = norm_mean_DEXT_mean, yend = norm_mean_DEXT_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = summary_stats_norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
                aes(x = oScarletpositive, y = norm_mean_DEXT_mean,
                    ymin = norm_mean_DEXT_mean - norm_mean_DEXT_se,
                    ymax = norm_mean_DEXT_mean + norm_mean_DEXT_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("oScarletnegative", "oScarletpositive"),
    limits = c(-0.5, 1.5),
    expand = c(0,0)) +
  #ggtitle("norm mean DEXT of oScarletnegative vs. oScarletpositive") +
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
wilcox.test(subset(norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, oScarletpositive == FALSE)$norm_mean_DEXT,
            subset(norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, oScarletpositive == TRUE)$norm_mean_DEXT,
            paired = TRUE)

# Mean apoebRNA by oScarlet status
# Save at 200W x 300H
ggplot(data = norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, aes(x = oScarletpositive, y = norm_mean_apoebRNA)) +
  ylim(c(0, 7)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = summary_stats_norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
               aes(x = as.numeric(oScarletpositive) - 0.1, xend = as.numeric(oScarletpositive) + 0.1,
                   y = norm_mean_apoebRNA_mean, yend = norm_mean_apoebRNA_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = summary_stats_norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
                aes(x = oScarletpositive, y = norm_mean_apoebRNA_mean,
                    ymin = norm_mean_apoebRNA_mean - norm_mean_apoebRNA_se,
                    ymax = norm_mean_apoebRNA_mean + norm_mean_apoebRNA_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("oScarletnegative", "oScarletpositive"),
    limits = c(-0.5, 1.5),
    expand = c(0,0)) +
  #ggtitle("norm mean apoebRNA of oScarletnegative vs. oScarletpositive") +
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
wilcox.test(subset(norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, oScarletpositive == FALSE)$norm_mean_apoebRNA,
            subset(norm_by_animal_and_oScarlet_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, oScarletpositive == TRUE)$norm_mean_apoebRNA,
            paired = TRUE)

# Mean apoebRNA by DEXT status
# Save at 200W x 300H
ggplot(data = norm_by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, aes(x = DEXTpositive, y = norm_mean_apoebRNA)) +
  ylim(c(0, 4)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = summary_stats_norm_by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
               aes(x = as.numeric(DEXTpositive) - 0.1, xend = as.numeric(DEXTpositive) + 0.1,
                   y = norm_mean_apoebRNA_mean, yend = norm_mean_apoebRNA_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = summary_stats_norm_by_animal_and_DEXT_status_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
                aes(x = DEXTpositive, y = norm_mean_apoebRNA_mean,
                    ymin = norm_mean_apoebRNA_mean - norm_mean_apoebRNA_se,
                    ymax = norm_mean_apoebRNA_mean + norm_mean_apoebRNA_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("DEXTnegative", "DEXTpositive"),
    limits = c(-0.5, 1.5),
    expand = c(0,0)) +
  #ggtitle("norm mean apoebRNA of DEXTnegative vs. DEXTpositive") +
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

# Exclude images potentially confounded by injury sites (areas of high apoebRNA density)
excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- subset(excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, Notes != "apoebRNAhighdensity")

# Summarize by animal and DEXT, oScarlet, and apoebRNA status

# (1) Summarize by animal and DEXT status
by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Animal, DEXTpositive) %>%
  summarise(
    mean_DEXT = mean(Cell..DEXT.mean, na.rm = TRUE),
    mean_oScarlet = mean(Cell..oScarlet.mean, na.rm = TRUE),
    mean_apoebRNA = mean(Cell..apoebRNA.mean, na.rm = TRUE),
    Experiment = first(Experiment)
  )

# Same-experiment DEXTnegative (sE_DEXTnegative) averages
sE_DEXTnegative_avg_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Experiment) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(.x[DEXTpositive == FALSE], na.rm = TRUE), 
                   .names = "sE_DEXTnegative_avg_{.col}"))

# Create a norm df with normalized values 
norm_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  left_join(sE_DEXTnegative_avg_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, by = c("Experiment")) %>%
  mutate(across(starts_with("mean"), 
                ~ . / get(paste0("sE_DEXTnegative_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(DEXTpositive, Animal, Experiment, starts_with("norm_"))

# Summary statistics
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
summary_stats_norm_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- norm_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(DEXTpositive) %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))

# (2) Summarize by animal and oScarlet status
by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Animal, oScarletpositive) %>%
  summarise(
    mean_DEXT = mean(Cell..DEXT.mean, na.rm = TRUE),
    mean_oScarlet = mean(Cell..oScarlet.mean, na.rm = TRUE),
    mean_apoebRNA = mean(Cell..apoebRNA.mean, na.rm = TRUE),
    Experiment = first(Experiment)
  )

# Same-experiment oScarletnegative (sE_oScarletnegative) averages
sE_oScarletnegative_avg_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Experiment) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(.x[oScarletpositive == FALSE], na.rm = TRUE), 
                   .names = "sE_oScarletnegative_avg_{.col}"))

# Create a norm df with normalized values 
norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  left_join(sE_oScarletnegative_avg_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, by = c("Experiment")) %>%
  mutate(across(starts_with("mean"), 
                ~ . / get(paste0("sE_oScarletnegative_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(oScarletpositive, Animal, Experiment, starts_with("norm_"))

# Summary statistics
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
summary_stats_norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(oScarletpositive) %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))

# (3) Summarize by animal and apoebRNA status
by_animal_and_apoebRNA_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Animal, apoebRNApositive) %>%
  summarise(
    mean_DEXT = mean(Cell..DEXT.mean, na.rm = TRUE),
    mean_oScarlet = mean(Cell..oScarlet.mean, na.rm = TRUE),
    mean_apoebRNA = mean(Cell..apoebRNA.mean, na.rm = TRUE),
    Experiment = first(Experiment)
  )

# Same-experiment apoebRNAnegative (sE_apoebRNAnegative) averages
sE_apoebRNAnegative_avg_by_animal_and_apoebRNA_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_apoebRNA_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(Experiment) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(.x[apoebRNApositive == FALSE], na.rm = TRUE), 
                   .names = "sE_apoebRNAnegative_avg_{.col}"))

# Create a norm df with normalized values 
norm_by_animal_and_apoebRNA_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- by_animal_and_apoebRNA_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  left_join(sE_apoebRNAnegative_avg_by_animal_and_apoebRNA_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, by = c("Experiment")) %>%
  mutate(across(starts_with("mean"), 
                ~ . / get(paste0("sE_apoebRNAnegative_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(apoebRNApositive, Animal, Experiment, starts_with("norm_"))

# Summary statistics
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
summary_stats_norm_by_animal_and_apoebRNA_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df <- norm_by_animal_and_apoebRNA_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df %>%
  group_by(apoebRNApositive) %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))

# Fig 3a
# Mean DEXT by oScarlet status
# Save at 200W x 300H
ggplot(data = norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, aes(x = oScarletpositive, y = norm_mean_DEXT)) +
  ylim(c(0, 7)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#0b6623", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = summary_stats_norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
               aes(x = as.numeric(oScarletpositive) - 0.1, xend = as.numeric(oScarletpositive) + 0.1,
                   y = norm_mean_DEXT_mean, yend = norm_mean_DEXT_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = summary_stats_norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
                aes(x = oScarletpositive, y = norm_mean_DEXT_mean,
                    ymin = norm_mean_DEXT_mean - norm_mean_DEXT_se,
                    ymax = norm_mean_DEXT_mean + norm_mean_DEXT_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("oScarletnegative", "oScarletpositive"),
    limits = c(-0.5, 1.5),
    expand = c(0,0)) +
  #ggtitle("norm mean DEXT of oScarletnegative vs. oScarletpositive") +
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
wilcox.test(subset(norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, oScarletpositive == FALSE)$norm_mean_DEXT,
            subset(norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, oScarletpositive == TRUE)$norm_mean_DEXT,
            paired = TRUE)

# Fig 3a
# Mean apoebRNA by oScarlet status
# Save at 200W x 300H
ggplot(data = norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, aes(x = oScarletpositive, y = norm_mean_apoebRNA)) +
  ylim(c(0, 7)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = summary_stats_norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
               aes(x = as.numeric(oScarletpositive) - 0.1, xend = as.numeric(oScarletpositive) + 0.1,
                   y = norm_mean_apoebRNA_mean, yend = norm_mean_apoebRNA_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = summary_stats_norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
                aes(x = oScarletpositive, y = norm_mean_apoebRNA_mean,
                    ymin = norm_mean_apoebRNA_mean - norm_mean_apoebRNA_se,
                    ymax = norm_mean_apoebRNA_mean + norm_mean_apoebRNA_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("oScarletnegative", "oScarletpositive"),
    limits = c(-0.5, 1.5),
    expand = c(0,0)) +
  #ggtitle("norm mean apoebRNA of oScarletnegative vs. oScarletpositive") +
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
wilcox.test(subset(norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, oScarletpositive == FALSE)$norm_mean_apoebRNA,
            subset(norm_by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, oScarletpositive == TRUE)$norm_mean_apoebRNA,
            paired = TRUE)

# Fig 3S1a
# Mean oScarlet by DEXT status
# Save at 200W x 300H
ggplot(data = norm_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, aes(x = DEXTpositive, y = norm_mean_oScarlet)) +
  ylim(c(0, 7)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#ee2400", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = summary_stats_norm_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
               aes(x = as.numeric(DEXTpositive) - 0.1, xend = as.numeric(DEXTpositive) + 0.1,
                   y = norm_mean_oScarlet_mean, yend = norm_mean_oScarlet_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = summary_stats_norm_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
                aes(x = DEXTpositive, y = norm_mean_oScarlet_mean,
                    ymin = norm_mean_oScarlet_mean - norm_mean_oScarlet_se,
                    ymax = norm_mean_oScarlet_mean + norm_mean_oScarlet_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("DEXTnegative", "DEXTpositive"),
    limits = c(-0.5, 1.5),
    expand = c(0,0)) +
  #ggtitle("norm mean oScarlet of DEXTnegative vs. DEXTpositive") +
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

# Fig 3S1a
# Mean apoebRNA by DEXT status
# Save at 200W x 300H
ggplot(data = norm_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, aes(x = DEXTpositive, y = norm_mean_apoebRNA)) +
  ylim(c(0, 7)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = summary_stats_norm_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
               aes(x = as.numeric(DEXTpositive) - 0.1, xend = as.numeric(DEXTpositive) + 0.1,
                   y = norm_mean_apoebRNA_mean, yend = norm_mean_apoebRNA_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = summary_stats_norm_by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
                aes(x = DEXTpositive, y = norm_mean_apoebRNA_mean,
                    ymin = norm_mean_apoebRNA_mean - norm_mean_apoebRNA_se,
                    ymax = norm_mean_apoebRNA_mean + norm_mean_apoebRNA_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("DEXTnegative", "DEXTpositive"),
    limits = c(-0.5, 1.5),
    expand = c(0,0)) +
  #ggtitle("norm mean apoebRNA of DEXTnegative vs. DEXTpositive") +
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

# Spatial maps
ggplot(data = subset(measurements_DEXT_oScarlet_apoebRNA_df,
                     (Date == "20251122" & nImage == "07")),
       aes(x = Centroid.X.µm,
           y = -Centroid.Y.µm,
           color = oScarletpositive)) +
  geom_point(size = 1,
             position = position_jitter(width = 0.0, height = 0.0)) +
  theme_minimal()

# Scatterplots (exclude S4-DEXT but not images with apoebRNAhighdensity)

# oScarlet vs. DEXT
ggplot(measurements_DEXT_oScarlet_apoebRNA_df,
       aes(x = Cell..oScarlet.mean,
           y = Cell..DEXT.mean,
           color = oScarletpositive)) +
  geom_point(size = 1,
             position = position_jitter(width = 0.0, height = 0.0)) +
  scale_x_continuous(limits = c(0, 65535), labels = scales::comma) +
  scale_y_continuous(limits = c(0, 65535), labels = scales::comma) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "none")

# oScarlet vs. apoebRNA
ggplot(measurements_DEXT_oScarlet_apoebRNA_df,
       aes(x = Cell..oScarlet.mean,
           y = Cell..apoebRNA.mean,
           color = oScarletpositive)) +
  geom_point(size = 1,
             position = position_jitter(width = 0.0, height = 0.0)) +
  scale_x_continuous(limits = c(0, 65535), labels = scales::comma) +
  scale_y_continuous(limits = c(0, 65535), labels = scales::comma) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "none")

# DEXT vs. apoebRNA
ggplot(excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df,
       aes(x = Cell..DEXT.mean,
           y = Cell..apoebRNA.mean)) +
  geom_point(size = 1,
             position = position_jitter(width = 0.0, height = 0.0)) +
  scale_x_continuous(limits = c(0, 65535), labels = scales::comma) +
  scale_y_continuous(limits = c(0, 65535), labels = scales::comma) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "none")

# Save source data tables
write.csv(by_animal_and_oScarlet_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, "SourceDataTable1.csv")
write.csv(by_animal_and_DEXT_status_excludeapoebRNAhighdensity_excludeS4DEXT_measurements_DEXT_oScarlet_apoebRNA_df, "SourceDataTable2.csv")
write.csv(pctpositivecells_DEXT_oScarlet_apoebRNA_df, "SourceDataTable3.csv")