rm(list=ls())
library(ggplot2)
library(tidyverse)
library(dplyr)

# Set working directory
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/1S1c/1S1c_QuPath_proj/")

# Read .csv
measurements_ELAVL3_df <- read.csv("ELAVL3_measurements.csv")

# Add columns with image data
measurements_ELAVL3_df <- measurements_ELAVL3_df %>%
  separate(col = Image, into = c("Date", "nImage", "Slide", "Animal", "Section", "Region", "Age", "Sex", "Channel405", "Channel488", "Channel546", "Channel647", "Magnification", "Notes"), sep = "_", remove = FALSE) %>%
  mutate(Notes = sub("\\.[^.]+$", "", Notes))

# Add a logical column for ELAVL3positive based on Cytoplasm..ELAVL3.mean
measurements_ELAVL3_df <- measurements_ELAVL3_df %>%
  group_by(Image) %>%
  mutate(
    ELAVL3positive_threshold = quantile(Cell..ELAVL3.mean, 0.10),
    ELAVL3positive = Cell..ELAVL3.mean > ELAVL3positive_threshold) %>%
  ungroup()

# Add a column for Genotype
measurements_ELAVL3_df <- measurements_ELAVL3_df %>%
  mutate(Genotype = case_when(
    startsWith(Animal, "S") ~ "SP-oScarlet",
    startsWith(Animal, "W") ~ "Wildtype",
    TRUE ~ NA_character_
  ))

# Add a column for Experiment
measurements_ELAVL3_df <- measurements_ELAVL3_df %>%
  mutate(Experiment = "Imaging1")

# Summarize by animal and ELAVL3 status
by_animal_and_ELAVL3_status_measurements_ELAVL3_df <- measurements_ELAVL3_df %>%
  group_by(Animal, ELAVL3positive) %>%
  summarise(
    n = n(),
    mean_ELAVL3 = mean(Cell..ELAVL3.mean, na.rm = TRUE),
    Genotype = first(Genotype),
    Experiment = first(Experiment)
  )

# Same-experiment wildtype (sE_wildtype_ELAVL3positive_ELAVL3positive) averages
sE_wildtype_ELAVL3positive_avg_by_animal_and_ELAVL3_status_measurements_ELAVL3_df <- by_animal_and_ELAVL3_status_measurements_ELAVL3_df %>%
  group_by(Experiment) %>%
  summarise(across(starts_with("mean"), 
                   ~ mean(.x[Genotype == "Wildtype" & ELAVL3positive == TRUE], na.rm = TRUE), 
                   .names = "sE_wildtype_ELAVL3positive_avg_{.col}"))

# Create a norm df with normalized values 
norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df <- by_animal_and_ELAVL3_status_measurements_ELAVL3_df %>%
  left_join(sE_wildtype_ELAVL3positive_avg_by_animal_and_ELAVL3_status_measurements_ELAVL3_df, by = c("Experiment")) %>%
  mutate(across(starts_with("mean"), 
                ~ . / get(paste0("sE_wildtype_ELAVL3positive_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(ELAVL3positive, Genotype, Animal, Experiment, starts_with("norm_"))

# Summary statistics
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
summary_stats_norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df <- norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df %>%
  group_by(Genotype, ELAVL3positive) %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))
summary_stats_norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df <- summary_stats_norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df %>%
  mutate(Genotype_factor = factor(Genotype))

# Fig 1S1c
# Mean ELAVL3 of ELAVL3positive by genotype
# Save at 200W x 300H
ggplot(data = subset(norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df, ELAVL3positive),
       aes(x = Genotype, y = norm_mean_ELAVL3)) +
  ylim(c(0, 2)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#0b6623", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df, ELAVL3positive),
               aes(x = as.numeric(Genotype_factor) - 0.1, xend = as.numeric(Genotype_factor) + 0.1,
                   y = norm_mean_ELAVL3_mean, yend = norm_mean_ELAVL3_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df, ELAVL3positive),
                aes(x = Genotype, y = norm_mean_ELAVL3_mean,
                    ymin = norm_mean_ELAVL3_mean - norm_mean_ELAVL3_se,
                    ymax = norm_mean_ELAVL3_mean + norm_mean_ELAVL3_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean ELAVL3 of ELAVL3positive by genotype") +
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

# Wilcoxon test (unpaired)
wilcox.test(subset(norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df, (ELAVL3positive & Genotype == "SP-oScarlet"))$norm_mean_ELAVL3,
            subset(norm_by_animal_and_ELAVL3_status_measurements_ELAVL3_df, (ELAVL3positive & Genotype == "Wildtype"))$norm_mean_ELAVL3,
            paired = FALSE)

# Spatial maps
ggplot(data = subset(measurements_ELAVL3_df,
                     (Date == "20250901" & nImage == "06")),
       aes(x = Centroid.X.µm,
           y = -Centroid.Y.µm,
           color = ELAVL3positive)) +
  geom_point(size = 1,
             position = position_jitter(width = 0.0, height = 0.0)) +
  theme_minimal()

# Save source data table
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/1S1c/")
write.csv(by_animal_and_ELAVL3_status_measurements_ELAVL3_df, "SourceDataTable1.csv")