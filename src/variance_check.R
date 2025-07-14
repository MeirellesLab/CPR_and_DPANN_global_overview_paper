library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)

load("data/rdata/homogeneity_of_variances.RData")

# Extract distances to centroids
distances <- homogeneity$distances
sample_ids <- rownames(homogeneity$vectors) 
groups <- homogeneity$group[names(distances)]  # Align names properly


distances_df <- data.frame(
  Group = homogeneity$group,
  Distance = homogeneity$distances
)

anova <- anova(homogeneity)


# Create a boxplot to visualize the distribution of variances
plot <- ggplot(distances_df, aes(x = Group, y = Distance,)) +
  geom_boxplot() +
  labs(x = "Ecosystem", y = "Distance to centroid", title = "Boxplot of Variances by Group") +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_blank())
#Write the anova results in the plot 
plot <- plot + 
  annotate("text", x = 1.5, y = max(distances_df$Distance) * 0.9, 
           label = paste("p =", format.pval(anova$`Pr(>F)`[1], digits = 2)), 
           size = 5, hjust = 0.5)
# Save the plot
ggsave("results/figures/variance_check.svg", plot, width = 8, height = 6, dpi = 300)

# save the anova results in txt file
writeLines(capture.output(print(anova)), "results/variance_check.txt")
