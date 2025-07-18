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
  "svglite" = "2.1.1",
  "grid" = "4.1.0"
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
phyla_abundances_wide_bonafide <- read_csv("data/treated/phyla_abundances_wide_bonafide.csv") %>%
  treatment()
phyla_abundances_wide_cpr <- read_csv("data/treated/phyla_abundances_wide_cpr.csv") %>%
  treatment()
phyla_abundances_wide_dpann <- read_csv("data/treated/phyla_abundances_wide_dpann.csv") %>%
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
load("data/rdata/nmds_bonafide.RData")
load("data/rdata/nmds_cpr.RData")
load("data/rdata/nmds_dpann.RData")
load("data/rdata/permanova_lifestyle.RData")
load("data/rdata/permanova_ecosystem.RData")
load("data/rdata/permanova_lifestyle_bonafide.RData")
load("data/rdata/permanova_ecosystem_bonafide.RData")
load("data/rdata/permanova_lifestyle_cpr.RData")
load("data/rdata/permanova_ecosystem_cpr.RData")
load("data/rdata/permanova_lifestyle_dpann.RData")
load("data/rdata/permanova_ecosystem_dpann.RData")

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

##################### barplots for richness and abundance ######################
source("src/util/draw_barplot_simple.R")
# Life Style -------------------------------------------------------------------
## Richness ----------------------------
bonafide_barplot_richness_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "Bonafide"),
  title = "",
  x_var = "life_style",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "Richness",
  title_x = "",
  legend_title = "Life Style",
  legend_position = "none",
  breaks = c(20, 30, 40),
  break_labels = c("20", "30", "40"),
  limits = c(20, 45),
  colors = life_style_colors
)
cpr_barplot_richness_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "CPR"),
  title = "",
  x_var = "life_style",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "Life Style",
  legend_title = "Life Style",
  legend_position = "none",
  breaks = c(80, 95, 105),
  break_labels = c("80", "95", "105"),
  limits = c(80, 110),
  colors = life_style_colors
)
dpann_barplot_richness_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "DPANN"),
  title = "",
  x_var = "life_style",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "",
  legend_title = "Life Style",
  legend_position = "none",
  breaks = c(0, 3.1, 6.4),
  break_labels = c("0", "3.1", "6.4"),
  colors = life_style_colors
)

# #Abundance ---------------------------
bonafide_barplot_abundance_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "Bonafide"),
  title = "Bonafide",
  x_var = "life_style",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "Relative abundance",
  title_x = "",
  legend_title = "Life style",
  legend_position = "none",
  breaks = c(0.90, 0.95, 1),
  break_labels = c("90", "95", "100"),
  limits = c(0.90, 1),
  colors = life_style_colors
)
cpr_barplot_abundance_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "CPR"),
  title = "CPR",
  x_var = "life_style",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "",
  legend_title = "Life style",
  legend_position = "none",
  breaks = c(0, 0.007, 0.0133),
  break_labels = c("    0", "0.7", "1.3"),
  colors = life_style_colors
)
dpann_barplot_abundance_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microgroup == "DPANN"),
  title = "DPANN",
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
bonafide_barplot_richness_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "Bonafide"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "Richness",
  title_x = "",
  legend_title = "Ecosystem",
  legend_position = "top",
  breaks = c(35, 40, 45),
  break_labels = c("35", "40", "45"),
  limits = c(35, 45),
  colors = ecosystem_colors
)
cpr_barplot_richness_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "CPR"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "Ecosystem",
  legend_title = "ecosystem",
  legend_position = "none",
  breaks = c(80, 95, 110),
  break_labels = c("80", "95", "110"),
  limits = c(80, 110),
  colors = ecosystem_colors
)
dpann_barplot_richness_ecosystem <- draw_barplot_simple(
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
bonafide_barplot_abundance_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "Bonafide"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "Relative abundance",
  title_x = "",
  legend_title = "",
  legend_position = "none",
  breaks = c(0.9, 0.95, 1),
  break_labels = c("90", "95", "100"),
  limits = c(0.9, 1),
  colors = ecosystem_colors
)
cpr_barplot_abundance_ecosystem <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microgroup == "CPR"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "",
  legend_title = "ecosystem",
  legend_position = "none",
  breaks = c(0, 0.035, 0.07),
  break_labels = c("    0", "3.5", "7"),
  colors = ecosystem_colors
)
dpann_barplot_abundance_ecosystem <- draw_barplot_simple(
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


# # ##################### violinplots for richness and abundance ######################
# source("src/util/draw_violinplot.R")
# # Life Style -------------------------------------------------------------------
# ## Richness ----------------------------

# # Create a column with the richness by sample in each microgroup dataframe.
# # Bonafide
# bonafide_taxa_cols <- 2:(ncol(phyla_abundances_wide_bonafide)-5)
# phyla_abundances_wide_bonafide$richness <- rowSums(phyla_abundances_wide_bonafide[, bonafide_taxa_cols] > 0)


# #CPR
# cpr_taxa_cols <- 7:(ncol(phyla_abundances_wide_cpr))
# phyla_abundances_wide_cpr$richness <- rowSums(phyla_abundances_wide_cpr[, cpr_taxa_cols] > 0)

# #DPANN
# dpann_taxa_cols <- 7:(ncol(phyla_abundances_wide_dpann))
# phyla_abundances_wide_dpann$richness <- rowSums(phyla_abundances_wide_dpann[, dpann_taxa_cols] > 0)


# bonafide_violinplot_richness_lifestyle <- draw_violinplot(
#   data = phyla_abundances_wide_bonafide,
#   title = "Culturable",
#   x_var = "life_style",
#   y_var = "richness",
#   title_y = "Richness",
#   title_x = "",
#   legend_title = "Life Style",
#   legend_position = "none",
#   breaks = c(0, 22, 44),
#   break_labels = c("0", "22", "44"),
#   colors = life_style_colors,
#   add_jitter = FALSE
# ) + ylim(30, 50)

# cpr_violinplot_richness_lifestyle <- draw_violinplot(
#   data = phyla_abundances_wide_cpr,
#   title = "CPR",
#   x_var = "life_style",
#   y_var = "richness",
#   title_y = "Richness",
#   title_x = "",
#   legend_title = "Life Style",
#   legend_position = "none",
#   breaks = c(0, 22, 44),
#   break_labels = c("0", "22", "44"),
#   colors = life_style_colors,
#   add_jitter = FALSE
# ) + ylim(0, 125)

# dpann_violinplot_richness_lifestyle <- draw_violinplot(
#   data = phyla_abundances_wide_dpann,
#   title = "DPANN",
#   x_var = "life_style",
#   y_var = "richness",
#   title_y = "Richness",
#   title_x = "",
#   legend_title = "Life Style",
#   legend_position = "none",
#   breaks = c(0, 22, 44),
#   break_labels = c("0", "22", "44"),
#   colors = life_style_colors,
#   add_jitter = FALSE
# ) + ylim(0, 12)

# ##Abundance ---------------------------

# # Create a column with the total relative abundance of the microgroup in each sample

# # Bonafide
# phyla_abundances_wide_bonafide$total_abu <- rowSums(phyla_abundances_wide_bonafide[, bonafide_taxa_cols]) * 100
# # CPR
# phyla_abundances_wide_cpr$total_abu <- rowSums(phyla_abundances_wide_cpr[, cpr_taxa_cols]) * 100
# # DPANN
# phyla_abundances_wide_dpann$total_abu <- rowSums(phyla_abundances_wide_dpann[, dpann_taxa_cols]) * 100


# bonafide_violinplot_abundance_lifestyle <- draw_violinplot(
#   data = phyla_abundances_wide_bonafide,
#   title = "Culturable",
#   x_var = "life_style",
#   y_var = "total_abu",
#   title_y = "Relative Abundance (%)",
#   title_x = "",
#   legend_title = "Life Style",
#   legend_position = "none",
#   breaks = c(0, 50, 100),
#   break_labels = c("0", "50", "100"),
#   colors = life_style_colors,
#   add_jitter = FALSE
# ) + ylim(90, 100)

# cpr_violinplot_abundance_lifestyle <- draw_violinplot(
#   data = phyla_abundances_wide_cpr,
#   title = "CPR",
#   x_var = "life_style",
#   y_var = "total_abu",
#   title_y = "Relative Abundance (%)",
#   title_x = "",
#   legend_title = "Life Style",
#   legend_position = "none",
#   breaks = c(0, 2, 10),
#   break_labels = c("0", "2", "10"),
#   colors = life_style_colors,
#   add_jitter = FALSE
# ) + ylim(0, 5)

# dpann_violinplot_abundance_lifestyle <- draw_violinplot(
#   data = phyla_abundances_wide_dpann,
#   title = "DPANN",
#   x_var = "life_style",
#   y_var = "total_abu",
#   title_y = "Relative Abundance (%)",
#   title_x = "",
#   legend_title = "Life Style",
#   legend_position = "none",
#   breaks = c(0, 0.05, 0.075),
#   break_labels = c("0", "0.05", "0.075"),
#   colors = life_style_colors,
#   add_jitter = FALSE
# ) + ylim(0, 0.075)

# # Ecossystem ------------------------------------------------------------------
# # Richness -----------------------------
# bonafide_violinplot_richness_ecosystem <- draw_violinplot(
#   data = phyla_abundances_wide_bonafide,
#   title = "Culturable",
#   x_var = "ecosystem",
#   y_var = "richness",
#   title_y = "Richness",
#   title_x = "",
#   legend_title = "Ecosystem",
#   legend_position = "none",
#   breaks = c(0, 50, 100),
#   break_labels = c("0", "50", "100"),
#   colors = ecosystem_colors,
#   add_jitter = FALSE
# ) + ylim(30, 50)

# cpr_violinplot_richness_ecosystem <- draw_violinplot(
#   data = phyla_abundances_wide_cpr,
#   title = "CPR",
#   x_var = "ecosystem",
#   y_var = "richness",
#   title_y = "Richness",
#   title_x = "",
#   legend_title = "Ecosystem",
#   legend_position = "none",
#   breaks = c(0, 50, 100),
#   break_labels = c("0", "50", "100"),
#   colors = ecosystem_colors,
#   add_jitter = FALSE
# ) + ylim(0, 110)

# dpann_violinplot_richness_ecosystem <- draw_violinplot(
#   data = phyla_abundances_wide_dpann,
#   title = "DPANN",
#   x_var = "ecosystem",
#   y_var = "richness",
#   title_y = "Richness",
#   title_x = "",
#   legend_title = "Ecosystem",
#   legend_position = "none",
#   breaks = c(0, 5, 10),
#   break_labels = c("0", "5", "10"),
#   colors = ecosystem_colors,
#   add_jitter = FALSE
# ) + ylim(0, 10)

# # Abundance ----------------------------
# bonafide_violinplot_abundance_ecosystem <- draw_violinplot(
#   data = phyla_abundances_wide_bonafide,
#   title = "Culturable",
#   x_var = "ecosystem",
#   y_var = "total_abu",
#   title_y = "Relative Abundance (%)",
#   title_x = "",
#   legend_title = "Ecosystem",
#   legend_position = "none",
#   breaks = c(0, 50, 100),
#   break_labels = c("0", "50", "100"),
#   colors = ecosystem_colors,
#   add_jitter = FALSE
# ) + ylim(90, 100)

# cpr_violinplot_abundance_ecosystem <- draw_violinplot(
#   data = phyla_abundances_wide_cpr,
#   title = "CPR",
#   x_var = "ecosystem",
#   y_var = "total_abu",
#   title_y = "Relative Abundance (%)",
#   title_x = "",
#   legend_title = "Ecosystem",
#   legend_position = "none",
#   breaks = c(0, 2, 7),
#   break_labels = c("0", "2", "7"),
#   colors = ecosystem_colors,
#   add_jitter = FALSE
# ) + ylim(0, 7)

# dpann_violinplot_abundance_ecosystem <- draw_violinplot(
#   data = phyla_abundances_wide_dpann,
#   title = "DPANN",
#   x_var = "ecosystem",
#   y_var = "total_abu",
#   title_y = "Relative Abundance (%)",
#   title_x = "",
#   legend_title = "Ecosystem",
#   legend_position = "none",
#   breaks = c(0, 0.05, 0.1),
#   break_labels = c("0", "0.05", "0.1"),
#   colors = ecosystem_colors,
#   add_jitter = FALSE
# ) + ylim(0, 0.075)

# Get legends ------------------------------------------------------------------

# Ecosystem legend -------------------------------------------------------------
# Convert the plot to a grob
ecosystem_grob_plot <- ggplotGrob(bonafide_barplot_richness_ecosystem)

# Find which part of the grob layout contains the legend
ecosystem_legend_index <- which(ecosystem_grob_plot$layout$name == "guide-box-top")

# Extract the legend
ecosystem_legend_grob <- ecosystem_grob_plot$grobs[[ecosystem_legend_index]]


# Lifestyle legend -------------------------------------------------------------
# Convert the plot to a grob
lifestyle_grob_plot <- ggplotGrob(dpann_barplot_abundance_lifestyle)

# Find which part of the grob layout contains the legend
lifestyle_legend_index <- which(lifestyle_grob_plot$layout$name == "guide-box-top")

# Extract the legend
lifestyle_legend_grob <- lifestyle_grob_plot$grobs[[lifestyle_legend_index]]

dpann_barplot_abundance_lifestyle <- dpann_barplot_abundance_lifestyle +
  theme(legend.position = "none")
bonafide_barplot_richness_ecosystem <- bonafide_barplot_richness_ecosystem +
  theme(legend.position = "none")

# Merge plots ------------------------------------------------------------------

# Barplots ---------------------------------------------------------------------
bonafide_barplot_abundance_ecosystem <- bonafide_barplot_abundance_ecosystem +
  theme(axis.text.x = element_blank())
  
bonafide_barplot_richness_ecosystem <- bonafide_barplot_richness_ecosystem +
  theme(axis.text.x = element_blank(),
        plot.title = element_blank())

bonafide_barplot_abundance_lifestyle <- bonafide_barplot_abundance_lifestyle +
  theme(axis.text.x = element_blank())

bonafide_barplot_richness_lifestyle <- bonafide_barplot_richness_lifestyle +
  theme(axis.text.x = element_blank())

cpr_barplot_abundance_lifestyle <- cpr_barplot_abundance_lifestyle +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank())
cpr_barplot_richness_lifestyle <- cpr_barplot_richness_lifestyle +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_blank())
cpr_barplot_abundance_ecosystem <- cpr_barplot_abundance_ecosystem +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank())
cpr_barplot_richness_ecosystem <- cpr_barplot_richness_ecosystem +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_blank())

dpann_barplot_abundance_lifestyle <- dpann_barplot_abundance_lifestyle +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank())
dpann_barplot_richness_lifestyle <- dpann_barplot_richness_lifestyle +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_blank())
dpann_barplot_abundance_ecosystem <- dpann_barplot_abundance_ecosystem +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank())
dpann_barplot_richness_ecosystem <- dpann_barplot_richness_ecosystem +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_blank())


barplot_abundance_lifestyle <- ggarrange(
  bonafide_barplot_abundance_lifestyle,
  cpr_barplot_abundance_lifestyle,
  dpann_barplot_abundance_lifestyle,
  ncol = 3,
  labels = c("A", "B", "C")
)
barplot_richness_lifestyle <- ggarrange(
  bonafide_barplot_richness_lifestyle,
  cpr_barplot_richness_lifestyle,
  dpann_barplot_richness_lifestyle,
  ncol = 3,
  labels = c("D", "E", "F")
)
barplot_richness_ecosystem <- ggarrange(
  bonafide_barplot_richness_ecosystem,
  cpr_barplot_richness_ecosystem,
  dpann_barplot_richness_ecosystem,
  ncol = 3,
  labels = c("J", "K", "L")
)
barplot_abundance_ecosystem <- ggarrange(
  bonafide_barplot_abundance_ecosystem,
  cpr_barplot_abundance_ecosystem,
  dpann_barplot_abundance_ecosystem,
  ncol = 3,
  labels = c("G", "H", "I")
)

# # Violinplots ------------------------------------------------------------------
# bonafide_violinplot_abundance_ecosystem <- bonafide_violinplot_abundance_ecosystem +
#   theme(axis.text.x = element_blank())
  
# bonafide_violinplot_richness_ecosystem <- bonafide_violinplot_richness_ecosystem +
#   theme(axis.text.x = element_blank(),
#         plot.title = element_blank())

# bonafide_violinplot_abundance_lifestyle <- bonafide_violinplot_abundance_lifestyle +
#   theme(axis.text.x = element_blank())

# bonafide_violinplot_richness_lifestyle <- bonafide_violinplot_richness_lifestyle +
#   theme(axis.text.x = element_blank(),
#         plot.title = element_blank())

# cpr_violinplot_abundance_lifestyle <- cpr_violinplot_abundance_lifestyle +
#   theme(axis.title.y = element_blank(),
#         axis.text.x = element_blank())
# cpr_violinplot_richness_lifestyle <- cpr_violinplot_richness_lifestyle +
#   theme(axis.title.y = element_blank(),
#         axis.text.x = element_blank(),
#         plot.title = element_blank())
# cpr_violinplot_abundance_ecosystem <- cpr_violinplot_abundance_ecosystem +
#   theme(axis.title.y = element_blank(),
#         axis.text.x = element_blank())
# cpr_violinplot_richness_ecosystem <- cpr_violinplot_richness_ecosystem +
#   theme(axis.title.y = element_blank(),
#         axis.text.x = element_blank(),
#         plot.title = element_blank())

# dpann_violinplot_abundance_lifestyle <- dpann_violinplot_abundance_lifestyle +
#   theme(axis.title.y = element_blank(),
#         axis.text.x = element_blank())
# dpann_violinplot_richness_lifestyle <- dpann_violinplot_richness_lifestyle +
#   theme(axis.title.y = element_blank(),
#         axis.text.x = element_blank(),
#         plot.title = element_blank())
# dpann_violinplot_abundance_ecosystem <- dpann_violinplot_abundance_ecosystem +
#   theme(axis.title.y = element_blank(),
#         axis.text.x = element_blank())
# dpann_violinplot_richness_ecosystem <- dpann_violinplot_richness_ecosystem +
#   theme(axis.title.y = element_blank(),
#         axis.text.x = element_blank(),
#         plot.title = element_blank())

# violinplot_abundance_lifestyle <- ggarrange(
#   bonafide_violinplot_abundance_lifestyle,
#   cpr_violinplot_abundance_lifestyle,
#   dpann_violinplot_abundance_lifestyle,
#   ncol = 3,
#   labels = c("A", "B", "C")
# )
# violinplot_richness_lifestyle <- ggarrange(
#   bonafide_violinplot_richness_lifestyle,
#   cpr_violinplot_richness_lifestyle,
#   dpann_violinplot_richness_lifestyle,
#   ncol = 3,
#   labels = c("D", "E", "F")
# )
# violinplot_richness_ecosystem <- ggarrange(
#   bonafide_violinplot_richness_ecosystem,
#   cpr_violinplot_richness_ecosystem,
#   dpann_violinplot_richness_ecosystem,
#   ncol = 3,
#   labels = c("J", "K", "L")
# )
# violinplot_abundance_ecosystem <- ggarrange(
#   bonafide_violinplot_abundance_ecosystem,
#   cpr_violinplot_abundance_ecosystem,
#   dpann_violinplot_abundance_ecosystem,
#   ncol = 3,
#   labels = c("G", "H", "I")
# )

############# nmds plot #######################################################

# NMDS general plot -----------------------------------------------------------
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
    text = element_text(size = unit(12, "points"), family = "Arial"),
    plot.title = element_text(
      hjust = 0.5, family = "Arial", size = unit(18, "points"), face = "bold"
    ),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      face = "bold", family = "Arial", size = unit(15, "points")
    ),
    legend.spacing.x = unit(0.1, "points"),
    legend.spacing.y = unit(0.1, "points"),
    axis.title.x = element_text(
      size = unit(15, "points"), face = "bold", family = "Arial"
    ),
    axis.title.y = element_text(
      size = unit(15, "points"), face = "bold", family = "Arial"
    ),
    legend.text = element_text(size = unit(12, "points"), family = "Arial")
  ) +
  geom_point(size = 2, shape = 20) +
  scale_color_manual(values = ecosystem_colors, name = "Ecosystem") +
  ggtitle("Whole community") +
  labs(
    x = "MDS1",
    y = "MDS2"
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
    size = unit(4, "points"),
    family = "Arial"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4, shape = 16)
  )) +
  scale_y_continuous(breaks = c(-1, -2, 0, 1, 2)) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))

# NMDS by microgroup ---------------------------------------------------------

microgroups <- phyla_abundances_long %>%
  dplyr::select(taxon, microgroup) %>%
  distinct()

# NMDS bonafide plot ----------------------------------------------------------
# Create dataframe ---------------------
bonafide_groups <- microgroups %>%
  filter(microgroup == "Bonafide")

# Create phyla_abundances_wide_bonafide df 
phyla_abundances_wide_bonafide <- phyla_abundances_wide %>%
  dplyr::select(samples, life_style, ecosystem, habitat, latitude, longitude, bonafide_groups$taxon)

nmds_bonafide_df <-
  cbind(
    phyla_abundances_wide_bonafide$ecosystem,
    as.data.frame(nmds_bonafide$points),
    phyla_abundances_wide_bonafide$samples
  ) %>%
  mutate(stress = nmds_bonafide$stress)
colnames(nmds_bonafide_df) <- c("ecosystem", "MDS1", "MDS2", "samples", "stress")

# Plot ---------------------------------
nmds_bonafide <-
  ggplot(nmds_bonafide_df, aes(x = MDS1, y = MDS2, color = ecosystem)) +
  theme_pubr() +
  theme(
    text = element_text(size = unit(12, "points"), family = "Arial"),
    plot.title = element_text(
      hjust = 0.5, family = "Arial", size = unit(18, "points"), face = "bold"
    ),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      face = "bold", family = "Arial", size = unit(15, "points")
    ),
    legend.spacing.x = unit(0.1, "points"),
    legend.spacing.y = unit(0.1, "points"),
    axis.title.x = element_text(
      size = unit(15, "points"), face = "bold", family = "Arial"
    ),
    axis.title.y = element_text(
      size = unit(15, "points"), face = "bold", family = "Arial"
    ),
    legend.text = element_text(size = unit(12, "points"), family = "Arial")
  ) +
  geom_point(size = 2, shape = 20) +
  scale_color_manual(values = ecosystem_colors, name = "Ecosystem") +
  ggtitle("Bonafide") +
  labs(
    x = "MDS1",
    y = "MDS2"
  ) +
  annotate(
    "text",
    x = min(nmds_bonafide_df$MDS1) + 2,
    y = max(nmds_bonafide_df$MDS2) - 0.5,
    label = paste(
      "Stress =",
      round(unique(nmds_bonafide_df$stress)[1], digits = 3),
      "\n",
      "Life style R² =",
      round(unique(permanova_lifestyle_bonafide$R2)[1], digits = 3),
      "\n",
      "p-value < ", round(unique(permanova_lifestyle_bonafide$"Pr(>F)")[1], digits = 4),
      "\n",
      "Ecossystem R² =",
      round(unique(permanova_ecosystem_bonafide$R2)[1], digits = 3),
      "\n",
      "p-value < ", round(unique(permanova_ecosystem_bonafide$"Pr(>F)")[1], digits = 4)
    ),
    size = unit(4, "points"),
    family = "Arial"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4, shape = 16)
  )) +
  scale_y_continuous(breaks = c(-1, -2, 0, 1, 2)) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))


# NMDS CPR plot --------------------------------------------------------------
# Create dataframe ---------------------
cpr_groups <- microgroups %>%
  filter(microgroup == "CPR")

# Create phyla_abundances_wide_bonafide df 
phyla_abundances_wide_cpr <- phyla_abundances_wide %>%
  dplyr::select(samples, life_style, ecosystem, habitat, latitude, longitude, cpr_groups$taxon)

nmds_cpr_df <-
  cbind(
    phyla_abundances_wide_cpr$ecosystem,
    as.data.frame(nmds_cpr$points),
    phyla_abundances_wide_cpr$samples
  ) %>%
  mutate(stress = nmds_cpr$stress)
colnames(nmds_cpr_df) <- c("ecosystem", "MDS1", "MDS2", "samples", "stress")

# Plot ---------------------------------
nmds_cpr <-
  ggplot(nmds_cpr_df, aes(x = MDS1, y = MDS2, color = ecosystem)) +
  theme_pubr() +
  theme(
    text = element_text(size = unit(12, "points"), family = "Arial"),
    plot.title = element_text(
      hjust = 0.5, family = "Arial", size = unit(18, "points"), face = "bold"
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
      size = unit(15, "points"), face = "bold", family = "Arial"
    ),
    axis.title.y = element_text(
      size = unit(15, "points"), face = "bold", family = "Arial"
    ),
    legend.text = element_text(size = unit(12, "points"), family = "Arial")
  ) +
  geom_point(size = 2, shape = 20) +
  scale_color_manual(values = ecosystem_colors, name = "Ecosystem") +
  ggtitle("CPR") +
  labs(
    x = "MDS1",
    y = "MDS2"
  ) +
  annotate(
    "text",
    x = min(nmds_cpr_df$MDS1) + 2,
    y = max(nmds_cpr_df$MDS2) - 0.5,
    label = paste(
      "Stress =",
      round(unique(nmds_cpr_df$stress)[1], digits = 3),
      "\n",
      "Life style R² =",
      round(unique(permanova_lifestyle_cpr$R2)[1], digits = 3),
      "\n",
      "p-value < ", round(unique(permanova_lifestyle_cpr$"Pr(>F)")[1], digits = 4),
      "\n",
      "Ecossystem R² =",
      round(unique(permanova_ecosystem_cpr$R2)[1], digits = 3),
      "\n",
      "p-value < ", round(unique(permanova_ecosystem_cpr$"Pr(>F)")[1], digits = 4)
    ),
    size = unit(4, "points"),
    family = "Arial"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4, shape = 16)
  )) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))

# NMDS DPANN plot ------------------------------------------------------------
# Create dataframe ---------------------
dpann_groups <- microgroups %>%
  filter(microgroup == "DPANN")

# Create phyla_abundances_wide_bonafide df 
phyla_abundances_wide_dpann <- phyla_abundances_wide %>%
  dplyr::select(samples, life_style, ecosystem, habitat, latitude, longitude, dpann_groups$taxon)

#Chek if there is any row with all zeros in phyla_abundances_persite_long_dpann
print("Checking if there is any row with all zeros in phyla_abundances_wide_dpann...")
print(sum(apply(phyla_abundances_wide_dpann[, -c(1, 2, 3, 4, 5, 6)], 1, sum) == 0))

print("Had 39 zeros rows, so we will remove them")
phyla_abundances_wide_dpann <- phyla_abundances_wide_dpann[apply(phyla_abundances_wide_dpann[, -c(1, 2, 3, 4, 5, 6)], 1, sum) != 0, ]


nmds_dpann_df <-
  cbind(
    phyla_abundances_wide_dpann$ecosystem,
    as.data.frame(nmds_dpann$points),
    phyla_abundances_wide_dpann$samples
  ) %>%
  mutate(stress = nmds_dpann$stress)

colnames(nmds_dpann_df) <- c("ecosystem", "MDS1", "MDS2", "samples", "stress")

# Plot ---------------------------------
nmds_dpann <- 
  ggplot(nmds_dpann_df, aes(x = MDS1, y = MDS2, color = ecosystem)) +
  theme_pubr() +
  theme(
    text = element_text(size = unit(12, "points"), family = "Arial"),
    plot.title = element_text(
      hjust = 0.5, family = "Arial", size = unit(18, "points"), face = "bold"
    ),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      face = "bold", family = "Arial", size = unit(15, "points")
    ),
    legend.spacing.x = unit(0.1, "points"),
    legend.spacing.y = unit(0.1, "points"),
    axis.title.x = element_text(
      size = unit(15, "points"), face = "bold", family = "Arial"
    ),
    axis.title.y = element_text(
      size = unit(15, "points"), face = "bold", family = "Arial"
    ),
    legend.text = element_text(size = unit(12, "points"), family = "Arial")
  ) +
  geom_point(size = 2, shape = 20) +
  scale_color_manual(values = ecosystem_colors, name = "Ecosystem") +
  ggtitle("DPANN") +
  labs(
    x = "MDS1",
    y = "MDS2"
  ) +
  annotate(
    "text",
    x = min(nmds_dpann_df$MDS1) + 2,
    y = max(nmds_dpann_df$MDS2) - 0.5,
    label = paste(
      "Stress =",
      round(unique(nmds_dpann_df$stress)[1], digits = 3),
      "\n",
      "Life style R² =",
      round(unique(permanova_lifestyle_dpann$R2)[1], digits = 3),
      "\n",
      "p-value < ", round(unique(permanova_lifestyle_dpann$"Pr(>F)")[1], digits = 4),
      "\n",
      "Ecossystem R² =",
      round(unique(permanova_ecosystem_dpann$R2)[1], digits = 3),
      "\n",
      "p-value < ", round(unique(permanova_ecosystem_dpann$"Pr(>F)")[1], digits = 4)
    ),
    size = unit(4, "points"),
    family = "Arial"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4, shape = 16)
  )) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))

############################### Merge plots ####################################
result_dir <- "results/figures/"
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}

## Panel 1 (all NMDS)

## extract nmds allsamples legend
nmds_allsamples_grob_plot <- ggplotGrob(nmds_allsamples)

# Find which part of the grob layout contains the legend
nmds_allsamples_grob_index <- which(nmds_allsamples_grob_plot$layout$name == "guide-box-right")

# Extract the legend
nmds_allsamples_legend_grob <- nmds_allsamples_grob_plot$grobs[[nmds_allsamples_grob_index]]


nmds_bonafide <- nmds_bonafide + theme(legend.position = "none")
nmds_cpr <- nmds_cpr + theme(legend.position = "none")
nmds_dpann <- nmds_dpann + theme(legend.position = "none")
nmds_allsamples <- nmds_allsamples + theme(legend.position = "none")

nmds_panel <- plot_grid(
  nmds_allsamples, nmds_bonafide, nmds_cpr, nmds_dpann,
  ncol = 2,
  nrow = 2,
  rel_widths = c(1, 1),
  rel_heights = c(1, 1),
  labels = c("B", "C", "D", "E"),
  label_size = 13,
  label_fontfamily = "Arial"
)

nmds_panel_legend <- plot_grid(
  nmds_panel,
  nmds_allsamples_legend_grob,
  ncol = 2,
  rel_widths = c(1, 0.2)
# rel_widths = c(1, 0.1)
)


panel_1 <- plot_grid(worldmap, nmds_panel_legend,
 ncol = 2,
 rel_widths = c(1, 0.85),
 labels = c("A", ""),
 label_size = 13,
 label_fontfamily = "Arial"
)

# Saving panel 1
ggsave(
 paste0(result_dir, "panel_1.svg"),
 plot = panel_1,
 width = 60,
 height = 20,
 units = "cm"
)
ggsave(
 paste0(result_dir, "panel_1.png"),
 plot = panel_1,
 width = 60,
 height = 20,
 units = "cm"
)

ggsave(
 paste0(result_dir, "panel_1.pdf"),
 plot = panel_1,
 width = 60,
 height = 20,
 units = "cm"
)

# Panel 2 (Richnes and abundance graphs)

# With barplot
ls_barplot <- plot_grid(
 barplot_abundance_lifestyle, barplot_richness_lifestyle, 
 ncol = 1,
 rel_widths = c(1, 1),
 label_size = 13,
 label_fontfamily = "Arial"
)

ls_barplot_legend <- cowplot::plot_grid(
  ls_barplot,
  lifestyle_legend_grob,
  ncol = 1,
  rel_heights = c(1, 0.3)
)

eco_barplot <- plot_grid(
 barplot_abundance_ecosystem, barplot_richness_ecosystem, 
 ncol = 1, 
 rel_heights = c(1, 1),
 label_fontfamily = "Arial",
 label_size = 13
)

eco_barplot_legend <- cowplot::plot_grid(
  eco_barplot,
  ecosystem_legend_grob,
  ncol = 1,
  rel_heights = c(1, 0.3)
)


panel_3_barplot <- plot_grid(
 ls_barplot_legend,
 eco_barplot_legend,
 nrow = 2,
 rel_widths = c(1, 1),
 rel_heights = c(1, 1)
) 



ggsave(
  paste0(result_dir, "panel_3_barplot.svg"),
  plot = panel_3_barplot,
  width = 27.6,
  height = 23,
  units = "cm"
)

ggsave(
  paste0(result_dir, "panel_3_barplot.png"),
  plot = panel_3_barplot,
  width = 27.6,
  height = 23,
  units = "cm"
)

ggsave(
  paste0(result_dir, "panel_3_barplot.pdf"),
  plot = panel_3_barplot,
  width = 27.6,
  height = 23,
  units = "cm"
)

# # With violinplot
# top_right_violinplot <- plot_grid(
#  violinplot_abundance_lifestyle, violinplot_richness_lifestyle, 
#  ncol = 1, align = "hv",
#  rel_widths = c(1, 1),
#  label_size = 13,
#  label_fontfamily = "Arial"
# )
# bottom_right_violinplot <- plot_grid(
#  violinplot_abundance_ecosystem, violinplot_richness_ecosystem, 
#  ncol = 1, align = "hv",
#  rel_heights = c(1, 1),
#  label_fontfamily = "Arial",
#  label_size = 13
# )

# bottom_right_violinplot_legend <- cowplot::plot_grid(
#   legends,
#   bottom_right_violinplot,
#   ncol = 1,
#   rel_heights = c(0.1, 1),
#   rel_widths = c(0.1, 1),
#   align = "h"
# )

# panel_2_violinplot <- plot_grid(
#  top_right_violinplot,
#  bottom_right_violinplot_legend,
#  nrow = 2,
#  rel_widths = c(1, 1),
#  rel_heights = c(1, 1)
# ) 


# ggsave(
#   paste0(result_dir, "panel_2_violinplot.svg"),
#   plot = panel_2_violinplot,
#   width = 20,
#   height = 26,
#   units = "cm"
# )

# ggsave(
#   paste0(result_dir, "panel_2_violinplot.png"),
#   plot = panel_2_violinplot,
#   width = 20,
#   height = 26,
#   units = "cm"
# )

# ggsave(
#   paste0(result_dir, "panel_2_violinplot.pdf"),
#   plot = panel_2_violinplot,
#   width = 20,
#   height = 26,
#   units = "cm"
# )

# ggsave(
#   paste0(result_dir, "panel_2_violinplot_wide.pdf"),
#   plot = panel_2_violinplot,
#   width = 30,
#   height = 26,
#   units = "cm"
# )

# Panel 4 (GAM plots)
panel_2 <- plot_grid(
  cpr_latitude_plot, dpann_latitude_plot, bonafide_latitude_plot,
  ncol = 3,
  rel_widths = c(1, 1, 1),
  labels = c("A", "B", "C"),
  label_size = 13,
  label_fontfamily = "Arial"
)

ggsave(
  paste0(result_dir, "panel_2.svg"),
  plot = panel_2,
  width = 18,
  height = 6,
  units = "cm"
)


