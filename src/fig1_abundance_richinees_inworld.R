#' @title Panel 1
#' @description This script generates the panel 1 of the paper. It depends on
#'  the outputs from "R/1_general_prospect.R".
#' @author Camilo Ferreira, Leonardo Brait, Felipe Alexandre, Pablo Viana.
#' @date 2022

################################# Environment ##################################
set.seed(201094)

source("src/util/install_and_load.R")
install_and_load(libs = c(
  "tidyverse" = "2.0.0",
  "ggpubr" = "0.6.0",
  "rnaturalearth" = "0.1.0",
  "cowplot" = "1.0.0",
  "extrafont" = "0.19",
  "svglite" = "2.1.1"
))
font_import(pattern = "Arial", prompt = FALSE)
extrafont::loadfonts(device = c("all"))
source("src/util/visualization_settings.R")
source("src/util/visualization_treatment.R")

################################ Load data #####################################
phyla_abundances_long <- read_csv("data/treated/phyla_abundances_long.csv") %>%
  treatment()
phyla_abundances_wide <- read_csv("data/treated/phyla_abundances_wide.csv") %>%
  treatment()
microgroups_prevalence_persite <- read_csv(
  "data/summaries/microgroups_prevalence_persite.csv"
) %>%
  treatment()
lifestyle_microgroups_prevalence <- read_csv(
  "data/summaries/lifestyle_microgroups_prevalence.csv"
) %>%
  treatment()
ecosystem_microgroups_prevalence <- read_csv(
  "data/summaries/ecosystem_microgroups_prevalence.csv"
) %>%
  treatment()
load("data/rdata/nmds.RData")
load("data/rdata/permanova_lifestyle.RData")
load("data/rdata/permanova_ecosystem.RData")

############################ World Map of samples ##############################
grouped_samples <- phyla_abundances_long %>%
  group_by(latitude, longitude, ecosystem) %>%
  summarise(nsamples = n()) %>%
  ungroup() %>%
  mutate(
    nsamples = as.numeric(nsamples),
    gradient = case_when(
      latitude < -66.5 ~ "South Pole",
      latitude > 66.5 ~ "North Pole",
      latitude >= -66.5 & latitude <= -35.5 ~ "South Temperate",
      latitude >= -35.5 & latitude <= -23.5 ~ "South Subtropical",
      latitude >= -23.5 & latitude <= 0 ~ "South Tropical",
      latitude >= 0 & latitude <= 23.5 ~ "North Tropical",
      latitude >= 23.5 & latitude <= 35.5 ~ "North Subtropical",
      latitude >= 35.5 & latitude <= 66.5 ~ "NorthTemperate"
    )
  )
worldmap <- ggplot(ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(fill = "#f9f7f0") +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  geom_point(
    data = grouped_samples,
    aes(x = longitude, y = latitude, size = nsamples, color = ecosystem),
    alpha = 0.6,
    shape = 19
  ) +
  scale_color_manual(values = ecosystem_colors)

############################# Richness vs Latitude #############################
source("src/util/draw_latitude_gam.R")
bonafide_latitude_plot <- draw_latitude_gam(
  data = microgroups_prevalence_persite %>% filter(microgroup == "Bonafide"),
  breaks = c(41.4, 42.42, 43.5),
  title_y = "Richness",
  main_title = "Bonafide",
  title_x = ""
)
cpr_latitude_plot <- draw_latitude_gam(
  data = microgroups_prevalence_persite %>% filter(microgroup == "CPR"),
  breaks = c(89.5, 99.5, 110.4),
  title_y = "",
  main_title = "CPR",
  title_x = "Latitude"
)
dpann_latitude_plot <- draw_latitude_gam(
  data = microgroups_prevalence_persite %>% filter(microgroup == "DPANN"),
  breaks = c(3.4, 5.5, 7.75),
  title_y = "",
  main_title = "DPANN",
  title_x = ""
)
plot_latitude <- plot_grid(
  bonafide_latitude_plot,
  cpr_latitude_plot,
  dpann_latitude_plot,
  ncol = 3,
  rel_widths = c(1, 1, 1)
)

##################### barplots for richness and abundance ######################
source("src/util/draw_barplot_simple.R")
# Life Style -------------------------------------------------------------------
## Richness ----------------------------
bonafide_barplot_richness_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "Bonafide"),
  title = "Bonafide",
  x_var = "life_style",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "Richness",
  title_x = "",
  legend_title = "Life Style",
  legend_position = "none",
  breaks = c(0, 22, 44),
  break_labels = c("    0", "22", "44"),
  colors = life_style_colors
)
cpr_barplot_richness_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "CPR"),
  title = "CPR",
  x_var = "life_style",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "",
  legend_title = "Life Style",
  legend_position = "none",
  breaks = c(0, 54, 105),
  break_labels = c("    0", "54", "105"),
  colors = life_style_colors
)
dpann_barplot_richness_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "DPANN"),
  title = "DPANN",
  x_var = "life_style",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "",
  legend_title = "Life Style",
  legend_position = "none",
  breaks = c(0, 3.1, 6.4),
  break_labels = c("    0", "3.1", "6.4"),
  colors = life_style_colors
)

# #Abundance ---------------------------
bonafide_abundance_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "Bonafide"),
  title = "",
  x_var = "life_style",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "Relative abundance (%)",
  title_x = "",
  legend_title = "Life style",
  legend_position = "none",
  breaks = c(0, 0.50, 1),
  break_labels = c("    0", "50", "100"),
  colors = life_style_colors
)
cpr_abundance_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "CPR"),
  title = "",
  x_var = "life_style",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "Life style",
  legend_title = "Life style",
  legend_position = "none",
  breaks = c(0, 0.007, 0.0133),
  break_labels = c("    0", "0.7", "1.3"),
  colors = life_style_colors
)
dpann_abundance_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "DPANN"),
  title = "",
  x_var = "life_style",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "",
  legend_title = "Life style",
  legend_position = "top",
  breaks = c(0, 0.000104, 0.00021),
  break_labels = c("    0", "0.01", "0.02"),
  colors = life_style_colors
)

# Ecossystem ------------------------------------------------------------------
# Richness -----------------------------
bonafide_richness_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "Bonafide"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "Richness",
  title_x = "",
  legend_title = "Ecosystem",
  legend_position = "top",
  breaks = c(0, 22, 44),
  break_labels = c("    0", "22", "44"),
  colors = ecosystem_colors
)
cpr_richness_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "CPR"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "",
  legend_title = "ecosystem",
  legend_position = "none",
  breaks = c(0, 55, 110),
  break_labels = c("    0", "55", "110"),
  colors = ecosystem_colors
)
dpann_richness_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "DPANN"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "",
  legend_title = "ecosystem",
  breaks = c(0, 4, 8),
  break_labels = c("    0", "4", "8"),
  colors = ecosystem_colors
)

# Abundance ----------------------------
bonafide_abundance_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "Bonafide"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "Relative abundance (%)",
  title_x = "",
  legend_title = "",
  legend_position = "none",
  breaks = c(0, 0.5, 1),
  break_labels = c("    0", "50", "100"),
  colors = ecosystem_colors
)
cpr_abundance_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "CPR"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "Ecosystem",
  legend_title = "ecosystem",
  legend_position = "none",
  breaks = c(0, 0.035, 0.07),
  break_labels = c("    0", "3.5", "7"),
  colors = ecosystem_colors
)
dpann_abundance_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "DPANN"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "",
  legend_title = "ecosystem",
  legend_position = "none",
  breaks = c(0, 0.00060, 0.00120),
  break_labels = c("    0", "0.06", "0.12"),
  colors = ecosystem_colors
)

# Get legends ------------------------------------------------------------------
ecosystem_legend <- cowplot::get_legend(bonafide_richness_ecosystem)
bonafide_richness_ecosystem <- bonafide_richness_ecosystem +
  theme(legend.position = "none")
life_style_legend <- cowplot::get_legend(dpann_abundance_lifestyle)
dpann_abundance_lifestyle <- dpann_abundance_lifestyle +
  theme(legend.position = "none")

# Merge plots ------------------------------------------------------------------
barplot_abundance_lifestyle <- ggarrange(
  bonafide_abundance_lifestyle,
  cpr_abundance_lifestyle,
  dpann_abundance_lifestyle,
  ncol = 3
)
barplot_richness_lifestyle <- ggarrange(
  bonafide_barplot_richness_lifestyle,
  cpr_barplot_richness_lifestyle,
  dpann_barplot_richness_lifestyle,
  ncol = 3
)
barplot_richness_ecosystem <- ggarrange(
  bonafide_richness_ecosystem,
  cpr_richness_ecosystem,
  dpann_richness_ecosystem,
  ncol = 3
)
barplot_abundance_ecosystem <- ggarrange(
  bonafide_abundance_ecosystem,
  cpr_abundance_ecosystem,
  dpann_abundance_ecosystem,
  ncol = 3
)

############# nmds plot #######################################################

# Create dataframe ---------------------
nmds_df <-
  cbind(
    phyla_abundances_wide$ecosystem,
    as.data.frame(nmds$points),
    phyla_abundances_wide$samples
  ) %>%
  mutate(stress = nmds$stress)
colnames(nmds_df) <- c("ecosystem", "MDS1", "MDS2", "samples", "stress")

# Plot ---------------------------------
nmds_allsamples <-
  ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = ecosystem)) +
  theme_pubr() +
  theme(
    text = element_text(size = unit(8, "points"), family = "Arial"),
    plot.title = element_text(
      hjust = 0.5, family = "Arial", size = unit(8, "points")
    ),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      face = "bold", family = "Arial", size = unit(8, "points")
    ),
    legend.spacing.x = unit(0.1, "points"),
    legend.spacing.y = unit(0.1, "points"),
    axis.title.x = element_text(
      size = unit(8, "points"), face = "bold", family = "Arial"
    ),
    axis.title.y = element_text(
      size = unit(8, "points"), face = "bold", family = "Arial"
    ),
    legend.text = element_text(size = unit(8, "points"), family = "Arial")
  ) +
  geom_point(size = 1, shape = 20) +
  scale_color_manual(values = ecosystem_colors, name = "Ecosystem") +
  ggtitle("") +
  labs(
    x = paste0(
      "MDS1 (",
      round(attr(nmds$species, "shrinkage")[1] * 100, digits = 2),
      "%)"
    ),
    y = paste0(
      "MDS2 (",
      round(attr(nmds$species, "shrinkage")[2] * 100, digits = 2),
      "%)"
    )
  ) +
  annotate(
    "text",
    x = min(nmds_df$MDS1) + 2,
    y = max(nmds_df$MDS2) - 0.5,
    label = paste(
      "Stress =",
      round(unique(nmds_df$stress)[1], digits = 3),
      "\n",
      "Life style R² =", round(unique(permanova_lifestyle$R2)[1], digits = 3),
      "\n",
      "p-value < ", round(unique(permanova_lifestyle$"Pr(>F)")[1], digits = 4),
      "\n",
      "Ecossystem R² =", round(unique(permanova_ecosystem$R2)[1], digits = 3),
      "\n",
      "p-value <", round(unique(permanova_ecosystem$"Pr(>F)")[1], digits = 4)
    ),
    size = unit(3, "points"),
    family = "Arial"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4, shape = 16)
  )) +
  annotation_custom(
    life_style_legend, xmin = 1, xmax = 3.1, ymin = 1, ymax = 2
  ) +
  scale_y_continuous(breaks = c(-1, -0.2, 0.6, 1.4)) +
  scale_x_continuous(breaks = c(-1, 0, 1, 2))

############################### Merge plots ####################################
result_dir <- "results/figures/"
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}

## Panel 1 (Latitude GAM plots and NMDS)

panel_1 <- plot_grid(
  nmds_allsamples, plot_latitude,
  ncol = 1,
  nrow = 2,
  rel_widths = c(1, 1),
  rel_heights = c(3, 1)
)

## Sabe panel 1
ggsave(
  paste0(result_dir, "panel1.svg"),
  plot = panel_1,
  width = 16,
  height = 18,
  units = "cm"
)



#top_right <- plot_grid(
#  barplot_richness_lifestyle, barplot_abundance_lifestyle,
#  ncol = 1, align = "hv",
#  rel_widths = c(1, 1),
#  labels = c("d", "e"),
#  label_size = 13,
#  label_fontfamily = "Arial"
#)
#bottom_right <- plot_grid(
#  barplot_richness_ecosystem, barplot_abundance_ecosystem,
#  ncol = 1, align = "hv",
#  rel_heights = c(1, 1),
#  labels = c("f", "g"),
#  label_fontfamily = "Arial",
#  label_size = 13
#)
#right_panel <- plot_grid(
#  top_right,
#  bottom_right,
#  nrow = 2,
#  rel_widths = c(1, 1)
#)
#left_panel <- plot_grid(worldmap, nmds_allsamples, plot_latitude,
#  ncol = 1,
#  rel_heights = c(1.5, 2, 1),
#  labels = c("a", "b", "c"),
#  label_size = 13,
#  label_fontfamily = "Arial"
#)

#panel1 <- plot_grid(left_panel, right_panel, ncol = 2, rel_widths = c(1, 1))
#ggsave(
#  paste0(result_dir, "panel1.pdf"),
#  plot = panel1,
#  width = 11,
#  height = 8,
#)
#ggsave(
#  paste0(result_dir, "panel1.svg"),
#  plot = panel1,
#  width = 11,
#  height = 8,
#)
#ggsave(
#  paste0(result_dir, "panel1.png"),
#  plot = panel1,
#  width = 11,
#  height = 8,
#)
