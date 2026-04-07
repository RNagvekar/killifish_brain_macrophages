rm(list = ls())
library(ggplot2)

# Read in normcounts_df
AtlasBrain_normcounts_df <- read.csv("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/5S1b/CountsNormDESeq2_Brain_240708.csv")

# Subset to elavl3
elavl3Brain_normcounts_df <- AtlasBrain_normcounts_df[AtlasBrain_normcounts_df$X == "elavl3", ]
rownames(elavl3Brain_normcounts_df) <- elavl3Brain_normcounts_df$X
elavl3Brain_normcounts_df <- elavl3Brain_normcounts_df[, -1]

# Read in design_df
Atlas_design_df <- read.csv("C:/Users/Rahul/Dropbox/Rahul shared folder/Paper_F31_Rahul/Paper/V2/sourcedata/5S1b/ExperimentDesign_allbatches_combined_v7.csv")

# Subset to brain only
AtlasBrain_design_df <- Atlas_design_df[Atlas_design_df$tissue == "Brain", ]

# For plot
AtlasBrain_forplot_df <- AtlasBrain_design_df[, colnames(AtlasBrain_design_df) %in% c("X", "sex", "age_days", "tissue")]
rownames(AtlasBrain_forplot_df) <- AtlasBrain_forplot_df$X
AtlasBrain_forplot_df <- AtlasBrain_forplot_df[, -1]
AtlasBrain_forplot_df$elavl3normcounts <- as.numeric(elavl3Brain_normcounts_df[1, rownames(AtlasBrain_forplot_df)])

# Regression
model <- lm(elavl3normcounts ~ age_days, data = AtlasBrain_forplot_df)
coef(model)
summary(model)

# Pearson correlation
cor(AtlasBrain_forplot_df$age_days, AtlasBrain_forplot_df$elavl3normcounts, method = "pearson")

# Fig 5S1b
# Save at 500W x 700H
ggplot(data = AtlasBrain_forplot_df,
       aes(x = age_days, y = elavl3normcounts)) +
  scale_x_continuous(limits = c(0, 170)) +
  scale_y_continuous(limits = c(0, 45000), labels = scales::comma) +
  geom_point(aes(shape = sex), size = 5, position = position_jitter(width = 0, height = 0), color = "grey", alpha = 0.5) +
  scale_shape_manual(values = c("M" = 15, "F" = 17)) +
  geom_line(linetype = "blank") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5, linetype = "dashed") +
  #ggtitle("elavl3 norm counts by age (brain)") +
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