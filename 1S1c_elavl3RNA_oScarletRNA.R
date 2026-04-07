rm(list=ls())
library(ggplot2)
library(tidyverse)
library(dplyr)

# Set working directory
setwd("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/1S1c/1S1c_QuPath_proj/")

# Read .csv
measurements_elavl3RNA_oScarletRNA_df <- read.csv("1S1c_elavl3RNA_oScarletRNA_measurements.csv")

# Add columns with image data
measurements_elavl3RNA_oScarletRNA_df <- measurements_elavl3RNA_oScarletRNA_df %>%
  separate(col = Image, into = c("Date", "nImage", "Slide", "Animal", "Section", "Region", "Age", "Sex", "Channel405", "Channel488", "Channel546", "Channel647", "Magnification", "Notes"), sep = "_", remove = FALSE) %>%
  mutate(Notes = sub("\\.[^.]+$", "", Notes))

# Fig 1S1c
# Scatterplot
# elavl3 RNA transcripts on x-axis
# oScarlet RNA transcripts on y-axis
# Save at 400W x 400H
ggplot(measurements_elavl3RNA_oScarletRNA_df,
       aes(x = Subcellular..Channel.2..Num.spots.estimated,
           y = Subcellular..Channel.4..Num.spots.estimated)) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  coord_fixed(ratio = 1) +
  geom_point(size = 2, color = "grey",
             position = position_jitter(width = 0, height = 0)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.line.x = element_line(color = "black", size = 0.25),
        axis.line.y = element_line(color = "black", size = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none")