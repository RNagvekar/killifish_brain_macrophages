rm(list = ls())
library(dplyr)
library(ggplot2)
library(stringr)

# Set the working directory to the folder containing subfolders
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/0_OVAuptake/scaledata/")

# Get a list of all .csv files in the folder and subfolders
csv_files <- list.files(path = ".", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

# Initialize an empty list to store data frames
scaledata_list <- list()

# Loop through each file and read it into R
for (file in csv_files) {
  
  # Read the .csv file
  thisfile_scaledata_df <- read.csv(file)
  
  # Extract the file name and subfolder name
  file_name <- basename(file)  # Get the file name
  subfolder_name <- dirname(file)  # Get the subfolder (experiment) path
  subfolder_name <- gsub(paste0(getwd(), "/"), "", subfolder_name)
  subfolder_name <- gsub("^\\./", "", subfolder_name)
  
  # Add new columns for file name and subfolder name
  thisfile_scaledata_df <- thisfile_scaledata_df %>%
    mutate(source_file = as.character(file_name), source_subfolder = as.character(subfolder_name))
  
  # Append the data frame to the list
  scaledata_list[[file]] <- thisfile_scaledata_df
}

# Combine all data frames into one
scaledata_df <- bind_rows(scaledata_list)

# Change column names
colnames(scaledata_df) <- c("FSC-A", "OVA", "oScarlet", "source_file", "source_subfolder")

# Add a column for sex (specified in source_subfolder name)
scaledata_df$sex <- ifelse(grepl("_males_", scaledata_df$source_subfolder), "M", "F")

# Add a column for experiment date (first 8 characters of source_subfolder name)
scaledata_df$date <- substr(scaledata_df$source_subfolder, 1, 8)

# Add a column for animal (substring between first two underscores in source_file + sex + date)
scaledata_df$animal <- paste0(sapply(strsplit(scaledata_df$source_file, "_"), function(x) x[2]), "_", scaledata_df$sex, "_", scaledata_df$date)

# Add a column for age (first character of animal name)
scaledata_df$age <- ifelse(grepl("^Y", scaledata_df$animal), "young", "old")

# Add a column for oScarlet normalized to FSC-A (x1000)
scaledata_df$oScarlet_norm_to_FSCA <- (scaledata_df$oScarlet/scaledata_df$`FSC-A`) * 1000

# Add a column for OVA normalized to FSC-A (x1000)
scaledata_df$OVA_norm_to_FSCA <- (scaledata_df$OVA/scaledata_df$`FSC-A`) * 1000

# Add a column for OVA normalized to oScarlet (x1000)
scaledata_df$OVA_norm_to_oScarlet <- (scaledata_df$OVA/scaledata_df$oScarlet) * 1000

# Add a column for oScarlet status
scaledata_df$oScarlet_status <- case_when(
  scaledata_df$oScarlet < 500 ~ "oScarletLOW",
  scaledata_df$oScarlet > 500 ~ "oScarletHIGH",
  TRUE ~ "oScarletMIDDLE"  # Default case
)               

# Add a column for OVA status
scaledata_df$OVA_status <- case_when(
  (scaledata_df$date == "20250826" & scaledata_df$OVA > 1000) ~ "OVA+",
  (scaledata_df$date == "20250827" & scaledata_df$OVA > 1000) ~ "OVA+",
  (scaledata_df$date == "20250902" & scaledata_df$OVA > 500) ~ "OVA+",
  TRUE ~ "OVA-"  # Default case
)

# Summarize by animal and oScarlet_status
by_animal_and_oScarlet_status_scaledata_df <- scaledata_df %>%
  group_by(animal, oScarlet_status) %>%
  summarise(
    `percent_OVA+` = sum(OVA_status == "OVA+")/n() * 100,
    mean_FSCA = mean(`FSC-A`, na.rm = TRUE),
    mean_oScarlet = mean(oScarlet, na.rm = TRUE),
    mean_oScarlet_norm_to_FSCA = mean(oScarlet_norm_to_FSCA, na.rm = TRUE),
    mean_OVA = mean(OVA, na.rm = TRUE),
    mean_OVA_norm_to_FSCA = mean(OVA_norm_to_FSCA, na.rm = TRUE),
    mean_OVA_norm_to_oScarlet = mean(OVA_norm_to_oScarlet, na.rm = TRUE),
    `mean_FSCA_of_OVA+` = mean(`FSC-A`[OVA_status == "OVA+"], na.rm = TRUE),
    `mean_oScarlet_of_OVA+` = mean(oScarlet[OVA_status == "OVA+"], na.rm = TRUE),
    `mean_OVA_of_OVA+` = mean(OVA[OVA_status == "OVA+"], na.rm = TRUE),
    `mean_OVA_norm_to_FSCA_of_OVA+` = mean(OVA_norm_to_FSCA[OVA_status == "OVA+"], na.rm = TRUE),
    `mean_OVA_norm_to_oScarlet_of_OVA+` = mean(OVA_norm_to_oScarlet[OVA_status == "OVA+"], na.rm = TRUE),
    age = first(age),
    sex = first(sex),
    date = first(date),
)

# Young before old (for plots)
by_animal_and_oScarlet_status_scaledata_df$age <- factor(by_animal_and_oScarlet_status_scaledata_df$age, levels = rev(sort(unique(by_animal_and_oScarlet_status_scaledata_df$age))))
str(by_animal_and_oScarlet_status_scaledata_df)

# Plot for non-normalized mean OVA of oScarletHIGH (males, 20250826 experiment)
ggplot(data = subset(by_animal_and_oScarlet_status_scaledata_df, (by_animal_and_oScarlet_status_scaledata_df$date == "20250826" & by_animal_and_oScarlet_status_scaledata_df$oScarlet_status == "oScarletHIGH")), aes(x = age, y = mean_OVA, color = age)) +
  ylim(c(0, 2500)) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0)) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#ee2400", "#ffb09c")) +
  ggtitle("mean OVA of oScarletHIGH (males only)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none")

# Plot for non-normalized mean OVA of oScarletHIGH (females, 20250827 experiment)
ggplot(data = subset(by_animal_and_oScarlet_status_scaledata_df, (by_animal_and_oScarlet_status_scaledata_df$date == "20250827" & by_animal_and_oScarlet_status_scaledata_df$oScarlet_status == "oScarletHIGH")), aes(x = age, y = mean_OVA, color = age)) +
  ylim(c(0, 2500)) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0)) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#ee2400", "#ffb09c")) +
  ggtitle("mean OVA of oScarletHIGH (females only)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none")

# Plot for non-normalized mean OVA of oScarletHIGH (males, 20250902 experiment)
ggplot(data = subset(by_animal_and_oScarlet_status_scaledata_df, (by_animal_and_oScarlet_status_scaledata_df$sex == "M" & by_animal_and_oScarlet_status_scaledata_df$date == "20250902" & by_animal_and_oScarlet_status_scaledata_df$oScarlet_status == "oScarletHIGH")), aes(x = age, y = mean_OVA, color = age)) +
  ylim(c(0, 500)) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0)) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#ee2400", "#ffb09c")) +
  ggtitle("mean OVA of oScarletHIGH (males)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none")

# Plot for non-normalized mean OVA of oScarletHIGH (females, 20250902 experiment)
ggplot(data = subset(by_animal_and_oScarlet_status_scaledata_df, (by_animal_and_oScarlet_status_scaledata_df$sex == "F" & by_animal_and_oScarlet_status_scaledata_df$date == "20250902" & by_animal_and_oScarlet_status_scaledata_df$oScarlet_status == "oScarletHIGH")), aes(x = age, y = mean_OVA, color = age)) +
  ylim(c(0, 500)) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0)) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#ee2400", "#ffb09c")) +
  ggtitle("mean OVA of oScarletHIGH (females)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none")

# ------------------------------------------------------------------
# oScarletHIGH vs. oScarletLOW in young comparison
# ------------------------------------------------------------------

# For normalization to oScarletLOW, first subset to young oScarletHIGH or oScarletLOW only
by_animal_young_scaledata_df <- subset(by_animal_and_oScarlet_status_scaledata_df,
                                       (by_animal_and_oScarlet_status_scaledata_df$age == "young" & (by_animal_and_oScarlet_status_scaledata_df$oScarlet_status == "oScarletLOW" | by_animal_and_oScarlet_status_scaledata_df$oScarlet_status == "oScarletHIGH")))

# Same-sex same-day oScarletLOW (ssd_oScarletLOW) averages
ssd_oScarletLOW_avg_by_animal_young_scaledata_df <- by_animal_young_scaledata_df %>%
  group_by(sex, date) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(.x[oScarlet_status == "oScarletLOW"], na.rm = TRUE), 
                   .names = "ssd_oScarletLOW_avg_{.col}"))

# Create a norm df with normalized values 
norm_by_animal_young_scaledata_df <- by_animal_young_scaledata_df %>%
  left_join(ssd_oScarletLOW_avg_by_animal_young_scaledata_df, by = c("sex", "date")) %>%
  mutate(across(starts_with("percent") | starts_with("mean"), 
                ~ . / get(paste0("ssd_oScarletLOW_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(oScarlet_status, sex, date, starts_with("norm_"))

# Summary statistics for norm df by oScarlet_status and sex
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}

summary_stats_norm_by_animal_young_scaledata_df <- norm_by_animal_young_scaledata_df %>%
  group_by(oScarlet_status, sex) %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))

# oScarletLOW before oScarletHIGH in norm df (for plots)
norm_by_animal_young_scaledata_df$oScarlet_status <- factor(norm_by_animal_young_scaledata_df$oScarlet_status, levels = rev(sort(unique(norm_by_animal_young_scaledata_df$oScarlet_status))))
str(norm_by_animal_young_scaledata_df)

# oScarletLOW before oScarletHIGH in summary stats df (for plots)
summary_stats_norm_by_animal_young_scaledata_df$oScarlet_status <- factor(summary_stats_norm_by_animal_young_scaledata_df$oScarlet_status, levels = rev(sort(unique(summary_stats_norm_by_animal_young_scaledata_df$oScarlet_status))))
str(summary_stats_norm_by_animal_young_scaledata_df)

# Plot for normalized mean OVA of oScarletLOW vs. oScarletHIGH (males)
ggplot(data = subset(norm_by_animal_young_scaledata_df, norm_by_animal_young_scaledata_df$sex == "M"), aes(x = oScarlet_status, y = norm_mean_OVA)) +
  ylim(c(0, 20)) +
  geom_point(shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                             summary_stats_norm_by_animal_young_scaledata_df$sex == "M"),
               aes(x = as.numeric(oScarlet_status) - 0.1, xend = as.numeric(oScarlet_status) + 0.1,
                   y = norm_mean_OVA_mean, yend = norm_mean_OVA_mean), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                              summary_stats_norm_by_animal_young_scaledata_df$sex == "M"),
                aes(x = oScarlet_status, y = norm_mean_OVA_mean,
                    ymin = norm_mean_OVA_mean - norm_mean_OVA_se,
                    ymax = norm_mean_OVA_mean + norm_mean_OVA_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean OVA of oScarletLOW vs. oScarletHIGH (males)") +
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
norm_mean_OVA_oScarletLOWvsoScarletHIGH_wilcox_result_M <- wilcox.test(norm_by_animal_young_scaledata_df$norm_mean_OVA[norm_by_animal_young_scaledata_df$sex == "M" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletLOW"],
                                                                       norm_by_animal_young_scaledata_df$norm_mean_OVA[norm_by_animal_young_scaledata_df$sex == "M" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletHIGH"],
                                                                       paired = TRUE)
norm_mean_OVA_oScarletLOWvsoScarletHIGH_pval_M <- norm_mean_OVA_oScarletLOWvsoScarletHIGH_wilcox_result_M$p.value

# Plot for normalized mean OVA of oScarletLOW vs. oScarletHIGH (females)
ggplot(data = subset(norm_by_animal_young_scaledata_df, norm_by_animal_young_scaledata_df$sex == "F"), aes(x = oScarlet_status, y = norm_mean_OVA)) +
  ylim(c(0, 20)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                             summary_stats_norm_by_animal_young_scaledata_df$sex == "F"),
               aes(x = as.numeric(oScarlet_status) - 0.1, xend = as.numeric(oScarlet_status) + 0.1,
                   y = norm_mean_OVA_mean, yend = norm_mean_OVA_mean), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                              summary_stats_norm_by_animal_young_scaledata_df$sex == "F"),
                aes(x = oScarlet_status, y = norm_mean_OVA_mean,
                    ymin = norm_mean_OVA_mean - norm_mean_OVA_se,
                    ymax = norm_mean_OVA_mean + norm_mean_OVA_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean OVA of oScarletLOW vs. oScarletHIGH (females)") +
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
norm_mean_OVA_oScarletLOWvsoScarletHIGH_wilcox_result_F <- wilcox.test(norm_by_animal_young_scaledata_df$norm_mean_OVA[norm_by_animal_young_scaledata_df$sex == "F" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletLOW"],
                                                                       norm_by_animal_young_scaledata_df$norm_mean_OVA[norm_by_animal_young_scaledata_df$sex == "F" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletHIGH"],
                                                                       paired = TRUE)
norm_mean_OVA_oScarletLOWvsoScarletHIGH_pval_F <- norm_mean_OVA_oScarletLOWvsoScarletHIGH_wilcox_result_F$p.value

# Plot for normalized percent OVA+ of oScarletLOW vs. oScarletHIGH (males)
ggplot(data = subset(norm_by_animal_young_scaledata_df, norm_by_animal_young_scaledata_df$sex == "M"), aes(x = oScarlet_status, y = `norm_percent_OVA+`)) +
  ylim(c(0, 70)) +
  geom_point(shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                             summary_stats_norm_by_animal_young_scaledata_df$sex == "M"),
               aes(x = as.numeric(oScarlet_status) - 0.1, xend = as.numeric(oScarlet_status) + 0.1,
                   y = `norm_percent_OVA+_mean`, yend = `norm_percent_OVA+_mean`), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                              summary_stats_norm_by_animal_young_scaledata_df$sex == "M"),
                aes(x = oScarlet_status, y = `norm_percent_OVA+_mean`,
                    ymin = `norm_percent_OVA+_mean` - `norm_percent_OVA+_se`,
                    ymax = `norm_percent_OVA+_mean` + `norm_percent_OVA+_se`),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm percent OVA+ of oScarletLOW vs. oScarletHIGH (males)") +
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
`norm_percent_OVA+_oScarletLOWvsoScarletHIGH_wilcox_result_M` <- wilcox.test(norm_by_animal_young_scaledata_df$`norm_percent_OVA+`[norm_by_animal_young_scaledata_df$sex == "M" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletLOW"],
                                                                             norm_by_animal_young_scaledata_df$`norm_percent_OVA+`[norm_by_animal_young_scaledata_df$sex == "M" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletHIGH"],
                                                                             paired = TRUE)
`norm_percent_OVA+_oScarletLOWvsoScarletHIGH_pval_M` <- `norm_percent_OVA+_oScarletLOWvsoScarletHIGH_wilcox_result_M`$p.value

# Plot for normalized percent OVA+ of oScarletLOW vs. oScarletHIGH (females)
ggplot(data = subset(norm_by_animal_young_scaledata_df, norm_by_animal_young_scaledata_df$sex == "F"), aes(x = oScarlet_status, y = `norm_percent_OVA+`)) +
  ylim(c(0, 70)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                             summary_stats_norm_by_animal_young_scaledata_df$sex == "F"),
               aes(x = as.numeric(oScarlet_status) - 0.1, xend = as.numeric(oScarlet_status) + 0.1,
                   y = `norm_percent_OVA+_mean`, yend = `norm_percent_OVA+_mean`), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                              summary_stats_norm_by_animal_young_scaledata_df$sex == "F"),
                aes(x = oScarlet_status, y = `norm_percent_OVA+_mean`,
                    ymin = `norm_percent_OVA+_mean` - `norm_percent_OVA+_se`,
                    ymax = `norm_percent_OVA+_mean` + `norm_percent_OVA+_se`),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm percent OVA+ of oScarletLOW vs. oScarletHIGH (females)") +
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
`norm_percent_OVA+_oScarletLOWvsoScarletHIGH_wilcox_result_F` <- wilcox.test(norm_by_animal_young_scaledata_df$`norm_percent_OVA+`[norm_by_animal_young_scaledata_df$sex == "F" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletLOW"],
                                                                             norm_by_animal_young_scaledata_df$`norm_percent_OVA+`[norm_by_animal_young_scaledata_df$sex == "F" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletHIGH"],
                                                                             paired = TRUE)
`norm_percent_OVA+_oScarletLOWvsoScarletHIGH_pval_F` <- `norm_percent_OVA+_oScarletLOWvsoScarletHIGH_wilcox_result_F`$p.value

# Plot for normalized mean OVA of OVA+ of oScarletLOW vs. oScarletHIGH (males)
ggplot(data = subset(norm_by_animal_young_scaledata_df, norm_by_animal_young_scaledata_df$sex == "M"), aes(x = oScarlet_status, y = `norm_mean_OVA_of_OVA+`)) +
  ylim(c(0, 2)) +
  geom_point(shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                             summary_stats_norm_by_animal_young_scaledata_df$sex == "M"),
               aes(x = as.numeric(oScarlet_status) - 0.1, xend = as.numeric(oScarlet_status) + 0.1,
                   y = `norm_mean_OVA_of_OVA+_mean`, yend = `norm_mean_OVA_of_OVA+_mean`), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                              summary_stats_norm_by_animal_young_scaledata_df$sex == "M"),
                aes(x = oScarlet_status, y = `norm_mean_OVA_of_OVA+_mean`,
                    ymin = `norm_mean_OVA_of_OVA+_mean` - `norm_mean_OVA_of_OVA+_se`,
                    ymax = `norm_mean_OVA_of_OVA+_mean` + `norm_mean_OVA_of_OVA+_se`),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean OVA of OVA+ of oScarletLOW vs. oScarletHIGH (males)") +
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
`norm_mean_OVA_of_OVA+_oScarletLOWvsoScarletHIGH_wilcox_result_M` <- wilcox.test(norm_by_animal_young_scaledata_df$`norm_mean_OVA_of_OVA+`[norm_by_animal_young_scaledata_df$sex == "M" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletLOW"],
                                                                                 norm_by_animal_young_scaledata_df$`norm_mean_OVA_of_OVA+`[norm_by_animal_young_scaledata_df$sex == "M" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletHIGH"],
                                                                                 paired = TRUE)
`norm_mean_OVA_of_OVA+_oScarletLOWvsoScarletHIGH_pval_M` <- `norm_mean_OVA_of_OVA+_oScarletLOWvsoScarletHIGH_wilcox_result_M`$p.value

# Plot for normalized mean OVA of OVA+ of oScarletLOW vs. oScarletHIGH (females)
ggplot(data = subset(norm_by_animal_young_scaledata_df, norm_by_animal_young_scaledata_df$sex == "F"), aes(x = oScarlet_status, y = `norm_mean_OVA_of_OVA+`)) +
  ylim(c(0, 2)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), color = "#8b008b", alpha = 0.5) +
  geom_line(linetype = "blank") +
  geom_segment(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                             summary_stats_norm_by_animal_young_scaledata_df$sex == "F"),
               aes(x = as.numeric(oScarlet_status) - 0.1, xend = as.numeric(oScarlet_status) + 0.1,
                   y = `norm_mean_OVA_of_OVA+_mean`, yend = `norm_mean_OVA_of_OVA+_mean`), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_young_scaledata_df,
                              summary_stats_norm_by_animal_young_scaledata_df$sex == "F"),
                aes(x = oScarlet_status, y = `norm_mean_OVA_of_OVA+_mean`,
                    ymin = `norm_mean_OVA_of_OVA+_mean` - `norm_mean_OVA_of_OVA+_se`,
                    ymax = `norm_mean_OVA_of_OVA+_mean` + `norm_mean_OVA_of_OVA+_se`),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean OVA of OVA+ of oScarletLOW vs. oScarletHIGH (females)") +
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

# Wilcoxon test
`norm_mean_OVA_of_OVA+_oScarletLOWvsoScarletHIGH_wilcox_result_F` <- wilcox.test(norm_by_animal_young_scaledata_df$`norm_mean_OVA_of_OVA+`[norm_by_animal_young_scaledata_df$sex == "F" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletLOW"],
                                                                                 norm_by_animal_young_scaledata_df$`norm_mean_OVA_of_OVA+`[norm_by_animal_young_scaledata_df$sex == "F" & norm_by_animal_young_scaledata_df$oScarlet_status == "oScarletHIGH"],
                                                                                 paired = TRUE)
`norm_mean_OVA_of_OVA+_oScarletLOWvsoScarletHIGH_pval_F` <- `norm_mean_OVA_of_OVA+_oScarletLOWvsoScarletHIGH_wilcox_result_F`$p.value

# ------------------------------------------------------------------
# Young vs. old oScarletHIGH comparison
# ------------------------------------------------------------------

# For same-sex same-day normalization, first subset to oScarletHIGH only
by_animal_oScarletHIGH_scaledata_df <- subset(by_animal_and_oScarlet_status_scaledata_df, by_animal_and_oScarlet_status_scaledata_df$oScarlet_status == "oScarletHIGH")

# Same-sex same-day young (ssd_young) averages
ssd_young_avg_by_animal_oScarletHIGH_scaledata_df <- by_animal_oScarletHIGH_scaledata_df %>%
  group_by(sex, date) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(.x[age == "young"], na.rm = TRUE), 
                   .names = "ssd_young_avg_{.col}"))

# Create a norm df with normalized values 
norm_by_animal_oScarletHIGH_scaledata_df <- by_animal_oScarletHIGH_scaledata_df %>%
  left_join(ssd_young_avg_by_animal_oScarletHIGH_scaledata_df, by = c("sex", "date")) %>%
  mutate(across(starts_with("percent") | starts_with("mean"), 
                ~ . / get(paste0("ssd_young_avg_", cur_column())), 
                .names = "norm_{.col}")) %>%
  dplyr::select(age, sex, date, starts_with("norm_"))

# Summary statistics for norm df by age and sex
# Standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}

summary_stats_norm_by_animal_oScarletHIGH_scaledata_df <- norm_by_animal_oScarletHIGH_scaledata_df %>%
  group_by(age, sex) %>%
  summarise(across(starts_with("norm"), 
                   list(mean = mean, se = se), 
                   .names = "{.col}_{.fn}"))

# Young before old in norm df (for plots)
#norm_by_animal_oScarletHIGH_scaledata_df$age <- factor(norm_by_animal_oScarletHIGH_scaledata_df$age, levels = rev(sort(unique(norm_by_animal_oScarletHIGH_scaledata_df$age))))
str(norm_by_animal_oScarletHIGH_scaledata_df)

# Young before old in summary stats df (for plots)
#summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$age <- factor(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$age, levels = rev(sort(unique(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$age))))
str(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df)

# Plot for normalized mean oScarlet of oScarletHIGH (males)
ggplot(data = subset(norm_by_animal_oScarletHIGH_scaledata_df, norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"), aes(x = age, y = norm_mean_oScarlet, color = age)) +
  ylim(c(0, 1.3)) +
  geom_point(shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#ee2400", "#ffb09c")) +
  geom_segment(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                             summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"),
               aes(x = as.numeric(age) - 0.1, xend = as.numeric(age) + 0.1,
                   y = norm_mean_oScarlet_mean, yend = norm_mean_oScarlet_mean), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                              summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"),
                aes(x = age, y = norm_mean_oScarlet_mean,
                    ymin = norm_mean_oScarlet_mean - norm_mean_oScarlet_se,
                    ymax = norm_mean_oScarlet_mean + norm_mean_oScarlet_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) + 
  #ggtitle("norm mean oScarlet of oScarletHIGH (males)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

# Wilcoxon test (unpaired)
norm_mean_oScarlet_wilcox_result_M <- wilcox.test(norm_by_animal_oScarletHIGH_scaledata_df$norm_mean_oScarlet[norm_by_animal_oScarletHIGH_scaledata_df$sex == "M" & norm_by_animal_oScarletHIGH_scaledata_df$age == "young"],
                                                  norm_by_animal_oScarletHIGH_scaledata_df$norm_mean_oScarlet[norm_by_animal_oScarletHIGH_scaledata_df$sex == "M" & norm_by_animal_oScarletHIGH_scaledata_df$age == "old"])
norm_mean_oScarlet_pval_M <- norm_mean_oScarlet_wilcox_result_M$p.value

# Plot for normalized mean oScarlet of oScarletHIGH (females)
ggplot(data = subset(norm_by_animal_oScarletHIGH_scaledata_df, norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"), aes(x = age, y = norm_mean_oScarlet, color = age)) +
  ylim(c(0, 1.3)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#ee2400", "#ffb09c")) +
  geom_segment(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                             summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"),
               aes(x = as.numeric(age) - 0.1, xend = as.numeric(age) + 0.1,
                   y = norm_mean_oScarlet_mean, yend = norm_mean_oScarlet_mean),
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                              summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"),
                aes(x = age, y = norm_mean_oScarlet_mean,
                    ymin = norm_mean_oScarlet_mean - norm_mean_oScarlet_se,
                    ymax = norm_mean_oScarlet_mean + norm_mean_oScarlet_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean oScarlet of oScarletHIGH (females)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

# Wilcoxon test (unpaired)
norm_mean_oScarlet_wilcox_result_F <- wilcox.test(norm_by_animal_oScarletHIGH_scaledata_df$norm_mean_oScarlet[norm_by_animal_oScarletHIGH_scaledata_df$sex == "F" & norm_by_animal_oScarletHIGH_scaledata_df$age == "young"],
                                                  norm_by_animal_oScarletHIGH_scaledata_df$norm_mean_oScarlet[norm_by_animal_oScarletHIGH_scaledata_df$sex == "F" & norm_by_animal_oScarletHIGH_scaledata_df$age == "old"])
norm_mean_oScarlet_pval_F <- norm_mean_oScarlet_wilcox_result_F$p.value

# Fig 5S2b
# Plot for normalized mean OVA of oScarletHIGH (males)
# Save at 200W x 300H
ggplot(data = subset(norm_by_animal_oScarletHIGH_scaledata_df, norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"), aes(x = age, y = norm_mean_OVA, color = age)) +
  ylim(c(0, 2)) +
  geom_point(shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#8b008b", "#ff80ff")) +
  geom_segment(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                             summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"),
               aes(x = as.numeric(age) - 0.1, xend = as.numeric(age) + 0.1,
                   y = norm_mean_OVA_mean, yend = norm_mean_OVA_mean), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                              summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"),
                aes(x = age, y = norm_mean_OVA_mean,
                    ymin = norm_mean_OVA_mean - norm_mean_OVA_se,
                    ymax = norm_mean_OVA_mean + norm_mean_OVA_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean OVA of oScarletHIGH (males)") +
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
norm_mean_OVA_wilcox_result_M <- wilcox.test(norm_by_animal_oScarletHIGH_scaledata_df$norm_mean_OVA[norm_by_animal_oScarletHIGH_scaledata_df$sex == "M" & norm_by_animal_oScarletHIGH_scaledata_df$age == "young"],
                                             norm_by_animal_oScarletHIGH_scaledata_df$norm_mean_OVA[norm_by_animal_oScarletHIGH_scaledata_df$sex == "M" & norm_by_animal_oScarletHIGH_scaledata_df$age == "old"])
norm_mean_OVA_pval_M <- norm_mean_OVA_wilcox_result_M$p.value

# Fig 5S2b
# Plot for normalized mean OVA of oScarletHIGH (females)
# Save at 200W x 300H
ggplot(data = subset(norm_by_animal_oScarletHIGH_scaledata_df, norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"), aes(x = age, y = norm_mean_OVA, color = age)) +
  ylim(c(0, 2)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#8b008b", "#ff80ff")) +
  geom_segment(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                             summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"),
               aes(x = as.numeric(age) - 0.1, xend = as.numeric(age) + 0.1,
                   y = norm_mean_OVA_mean, yend = norm_mean_OVA_mean), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                              summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"),
                aes(x = age, y = norm_mean_OVA_mean,
                    ymin = norm_mean_OVA_mean - norm_mean_OVA_se,
                    ymax = norm_mean_OVA_mean + norm_mean_OVA_se),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean OVA of oScarletHIGH (females)") +
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
norm_mean_OVA_wilcox_result_F <- wilcox.test(norm_by_animal_oScarletHIGH_scaledata_df$norm_mean_OVA[norm_by_animal_oScarletHIGH_scaledata_df$sex == "F" & norm_by_animal_oScarletHIGH_scaledata_df$age == "young"],
                                             norm_by_animal_oScarletHIGH_scaledata_df$norm_mean_OVA[norm_by_animal_oScarletHIGH_scaledata_df$sex == "F" & norm_by_animal_oScarletHIGH_scaledata_df$age == "old"])
norm_mean_OVA_pval_F <- norm_mean_OVA_wilcox_result_F$p.value

# Fig 5S2b
# Plot for normalized percent OVA+ of oScarletHIGH (males)
# Save at 200W x 300H
ggplot(data = subset(norm_by_animal_oScarletHIGH_scaledata_df, norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"), aes(x = age, y = `norm_percent_OVA+`, color = age)) +
  ylim(c(0, 2)) +
  geom_point(shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#8b008b", "#ff80ff")) +
  geom_segment(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                             summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"),
               aes(x = as.numeric(age) - 0.1, xend = as.numeric(age) + 0.1,
                   y = `norm_percent_OVA+_mean`, yend = `norm_percent_OVA+_mean`), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                              summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"),
                aes(x = age, y = `norm_percent_OVA+_mean`,
                    ymin = `norm_percent_OVA+_mean` - `norm_percent_OVA+_se`,
                    ymax = `norm_percent_OVA+_mean` + `norm_percent_OVA+_se`),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm percent OVA+ of oScarletHIGH (males)") +
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
`norm_percent_OVA+_wilcox_result_M` <- wilcox.test(norm_by_animal_oScarletHIGH_scaledata_df$`norm_percent_OVA+`[norm_by_animal_oScarletHIGH_scaledata_df$sex == "M" & norm_by_animal_oScarletHIGH_scaledata_df$age == "young"],
                                                   norm_by_animal_oScarletHIGH_scaledata_df$`norm_percent_OVA+`[norm_by_animal_oScarletHIGH_scaledata_df$sex == "M" & norm_by_animal_oScarletHIGH_scaledata_df$age == "old"])
`norm_percent_OVA+_pval_M` <- `norm_percent_OVA+_wilcox_result_M`$p.value

# Fig 5S2b
# Plot for normalized percent OVA+ of oScarletHIGH (females)
# Save at 200W x 300H
ggplot(data = subset(norm_by_animal_oScarletHIGH_scaledata_df, norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"), aes(x = age, y = `norm_percent_OVA+`, color = age)) +
  ylim(c(0, 2)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#8b008b", "#ff80ff")) +
  geom_segment(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                             summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"),
               aes(x = as.numeric(age) - 0.1, xend = as.numeric(age) + 0.1,
                   y = `norm_percent_OVA+_mean`, yend = `norm_percent_OVA+_mean`), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                              summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"),
                aes(x = age, y = `norm_percent_OVA+_mean`,
                    ymin = `norm_percent_OVA+_mean` - `norm_percent_OVA+_se`,
                    ymax = `norm_percent_OVA+_mean` + `norm_percent_OVA+_se`),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm percent OVA+ of oScarletHIGH (females)") +
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
`norm_percent_OVA+_wilcox_result_F` <- wilcox.test(norm_by_animal_oScarletHIGH_scaledata_df$`norm_percent_OVA+`[norm_by_animal_oScarletHIGH_scaledata_df$sex == "F" & norm_by_animal_oScarletHIGH_scaledata_df$age == "young"],
                                                   norm_by_animal_oScarletHIGH_scaledata_df$`norm_percent_OVA+`[norm_by_animal_oScarletHIGH_scaledata_df$sex == "F" & norm_by_animal_oScarletHIGH_scaledata_df$age == "old"])
`norm_percent_OVA+_pval_F` <- `norm_percent_OVA+_wilcox_result_F`$p.value

# Fig 5S2b
# Plot for normalized mean OVA of OVA+ of oScarletHIGH (males)
# Save at 200W x 300H
ggplot(data = subset(norm_by_animal_oScarletHIGH_scaledata_df, norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"), aes(x = age, y = `norm_mean_OVA_of_OVA+`, color = age)) +
  ylim(c(0, 2)) +
  geom_point(shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#8b008b", "#ff80ff")) +
  geom_segment(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                             summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"),
               aes(x = as.numeric(age) - 0.1, xend = as.numeric(age) + 0.1,
                   y = `norm_mean_OVA_of_OVA+_mean`, yend = `norm_mean_OVA_of_OVA+_mean`), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                              summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "M"),
                aes(x = age, y = `norm_mean_OVA_of_OVA+_mean`,
                    ymin = `norm_mean_OVA_of_OVA+_mean` - `norm_mean_OVA_of_OVA+_se`,
                    ymax = `norm_mean_OVA_of_OVA+_mean` + `norm_mean_OVA_of_OVA+_se`),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean OVA of OVA+ of oScarletHIGH (males)") +
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
`norm_mean_OVA_of_OVA+_wilcox_result_M` <- wilcox.test(norm_by_animal_oScarletHIGH_scaledata_df$`norm_mean_OVA_of_OVA+`[norm_by_animal_oScarletHIGH_scaledata_df$sex == "M" & norm_by_animal_oScarletHIGH_scaledata_df$age == "young"],
                                                       norm_by_animal_oScarletHIGH_scaledata_df$`norm_mean_OVA_of_OVA+`[norm_by_animal_oScarletHIGH_scaledata_df$sex == "M" & norm_by_animal_oScarletHIGH_scaledata_df$age == "old"])
`norm_mean_OVA_of_OVA+_pval_M` <- `norm_mean_OVA_of_OVA+_wilcox_result_M`$p.value

# Fig 5S2b
# Plot for normalized mean OVA of OVA+ of oScarletHIGH (females)
# Save at 200W x 300H
ggplot(data = subset(norm_by_animal_oScarletHIGH_scaledata_df, norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"), aes(x = age, y = `norm_mean_OVA_of_OVA+`, color = age)) +
  ylim(c(0, 2)) +
  geom_point(shape = 17, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(linetype = "blank") +
  scale_color_manual(values = c("#8b008b", "#ff80ff")) +
  geom_segment(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                             summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"),
               aes(x = as.numeric(age) - 0.1, xend = as.numeric(age) + 0.1,
                   y = `norm_mean_OVA_of_OVA+_mean`, yend = `norm_mean_OVA_of_OVA+_mean`), 
               color = "black", size = 1, alpha = 0.5) +
  geom_errorbar(data = subset(summary_stats_norm_by_animal_oScarletHIGH_scaledata_df,
                              summary_stats_norm_by_animal_oScarletHIGH_scaledata_df$sex == "F"),
                aes(x = age, y = `norm_mean_OVA_of_OVA+_mean`,
                    ymin = `norm_mean_OVA_of_OVA+_mean` - `norm_mean_OVA_of_OVA+_se`,
                    ymax = `norm_mean_OVA_of_OVA+_mean` + `norm_mean_OVA_of_OVA+_se`),
                width = 0.1, color = "black", size = 1, alpha = 0.5) +
  #ggtitle("norm mean OVA of OVA+ of oScarletHIGH (females)") +
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
`norm_mean_OVA_of_OVA+_wilcox_result_F` <- wilcox.test(norm_by_animal_oScarletHIGH_scaledata_df$`norm_mean_OVA_of_OVA+`[norm_by_animal_oScarletHIGH_scaledata_df$sex == "F" & norm_by_animal_oScarletHIGH_scaledata_df$age == "young"],
                                                       norm_by_animal_oScarletHIGH_scaledata_df$`norm_mean_OVA_of_OVA+`[norm_by_animal_oScarletHIGH_scaledata_df$sex == "F" & norm_by_animal_oScarletHIGH_scaledata_df$age == "old"])
`norm_mean_OVA_of_OVA+_pval_F` <- `norm_mean_OVA_of_OVA+_wilcox_result_F`$p.value

# Save source data table
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/5S2b")
write.csv(by_animal_and_oScarlet_status_scaledata_df, "SourceDataTable8.csv")