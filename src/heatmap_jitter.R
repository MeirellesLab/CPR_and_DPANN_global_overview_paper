#' @title Panel 3
#' @description Generate heatmap and jitter plots for panel 3.
#' @author Felipe Alexandre and Leonardo Brait
#' @date 2023
################################## Environment #################################
set.seed(201094)
source("src/util/install_and_load.R")
source("src/util/visualization_settings.R")
install_and_load(libs = c(
  "tidyverse" = "2.0.0",
  "ggpubr" = "0.2.3",
  "extrafont" = "0.19",
  "scales" = "1.0.0",
  "cowplot" = "1.0.0",
  "svglite" = "2.1.1"
))
font_import(pattern = "Arial", prompt = FALSE)
extrafont::loadfonts(device = c("all"))

################################## Load data ###################################
phyla_abundances <-
  read_csv("data/treated/phyla_abundances_long.csv")

simper_result <-
  read_csv("data/statistics/simper_ecosystem.csv") %>%
  rename(taxon = OTU) %>%
  separate(comparison, into = c("comparison_1", "comparison_2"), sep = "_")

radiation <- read_csv("data/taxa_compositions/radiation_phyla.csv")

###################### Order and treat names of phyla ##########################
source("src/util/visualization_treatment.R")
phyla_order <- phyla_abundances %>%
  group_by(taxon) %>%
  summarise(mean_abundance = mean(abundance)) %>%
  arrange(mean_abundance) %>%
  mutate(taxon = str_replace(.$taxon, "Candidatus", "Ca.")) %>%
  mutate(taxon = str_replace(.$taxon, "candidate", "Ca.")) %>%
  pull(taxon) %>%
  as.character()
phyla_abundances <- phyla_abundances %>%
  mutate(taxon = str_replace(.$taxon, "Candidatus", "Ca.")) %>%
  mutate(taxon = str_replace(.$taxon, "candidate", "Ca.")) %>%
  treatment() %>%
  mutate(taxon = factor(taxon, levels = phyla_order))

cpr_list <- phyla_abundances %>%
  filter(microgroup == "CPR") %>%
  pull(taxon) %>%
  as.character()
dpann_list <- phyla_abundances %>%
  filter(microgroup == "DPANN") %>%
  pull(taxon) %>%
  as.character()
simper_result <- simper_result %>%
  mutate(taxon = str_replace(.$taxon, "Candidatus", "Ca.")) %>%
  mutate(taxon = str_replace(.$taxon, "candidate", "Ca.")) %>%
  mutate(microgroup = factor(case_when(
    taxon %in% dpann_list ~ "DPANN",
    taxon %in% cpr_list ~ "CPR",
    TRUE ~ "Bonafide"
  ))) %>%
  mutate(taxon = factor(taxon, levels = phyla_order)) %>%
  treatment()

################################# Jitter plots #################################

jitter_candidate <-
  ggplot(
    data = subset(phyla_abundances, microgroup != "Bonafide"),
    aes(
      y = abundance * 100,   x = taxon,
      fill = ecosystem,      color = ecosystem_colors
    )
  ) +
  labs(x = "", y = "Relative abundance (%)") +
  geom_jitter(
    alpha = 0.25, aes(group = ecosystem, color = ecosystem), size = .8
  ) +
  geom_hline(
    yintercept = 0.01 * 100,   linetype = "dashed",
    colour = "#737373",      size = 1
  ) +
  annotate(
    "text", y = 2,    x = 12,   label = "Rarity Threshold",
    angle = 90,    size = 10,   colour = "#737373", family = "Arial"
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  guides(colour = guide_legend(
    override.aes = list(shape = 16, nrow = 1, alpha = 1, size = 6),
    nrow = 3,
    title.position = "left",
    title.hjust = 1,
  )) +
  theme(
    text = element_text(size = unit(20, "points"), family = "Arial"),
    axis.text.x = element_text(
      size = unit(20, "points"),  hjust = 1, family = "Arial"
    ),
    legend.key = element_rect(fill = "transparent"),
    axis.text.y = element_text(
      size = unit(20, "points"),  hjust = 1, family = "Arial"
    ),
    legend.key.size = unit(1.3, "points"),
    legend.direction = "horizontal",
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = unit(20, "points"), face = "bold", family = "Arial"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title = element_text(
      face = "bold", family = "Arial", size = unit(32, "points")
    ),
    panel.spacing = unit(0.1, "lines"),
    legend.title = element_text(
      face = "bold", family = "Arial", size = unit(20, "points")
    ),
    plot.margin = margin(t = 6, r = 0, b = 2, l = 0)
  ) +
  scale_color_manual(name = "Ecosystem", values = ecosystem_colors) +
  scale_fill_manual(name = "Ecosystem", values = ecosystem_colors) +
  coord_flip()

# Bonafide -----------------------------
jitter_bonafide <-
  ggplot(
    data = subset(phyla_abundances, microgroup == "Bonafide"),
    aes(
      y = abundance * 100,   x = taxon,
      fill = ecosystem,      color = ecosystem_colors
    )
  ) +
  labs(x = "", y = "Relative abundance (%)") +
  geom_jitter(
    alpha = 0.25, aes(group = ecosystem, color = ecosystem), size = .8
  ) +
  geom_hline(
    yintercept = 0.01 * 100,   linetype = "dashed",
    colour = "#737373",      size = 1
  ) +
  annotate(
    "text", y = 2,    x = 12,   label = "Rarity Threshold",
    angle = 90,    size = 10,   colour = "#737373", family = "Arial"
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  theme(
    text = element_text(size = unit(20, "points"), family = "Arial"),
    axis.text.x = element_text(
      size = unit(20, "points"),  hjust = 1, family = "Arial"
    ),
    axis.text.y = element_text(
      size = unit(20, "points"),  hjust = 1, family = "Arial"
    ),
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = unit(20, "points"), face = "bold", family = "Arial"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title = element_text(
      size = unit(20, "points"), face = "bold", family = "Arial"
    ),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "none"
  ) +
  scale_color_manual(name = "Ecosystem", values = ecosystem_colors) +
  scale_fill_manual(name = "Ecosystem", values = ecosystem_colors) +
  coord_flip()

################################### Heatmaps ###################################
# Create summarized df for heatmap that group by taxon and comparison_1 and create a column mean_contribution that is the mean of contribution~
simper_result_sum <- simper_result %>%
  group_by(comparison_1, taxon, microgroup) %>%
  summarise(mean_contribution = mean(contribution)) %>%
  ungroup()

## Candidate --------------------------

heatmap_contribution_candidate <- ggplot(data = subset(simper_result_sum, microgroup != "Bonafide")) +
  labs(x = "Ecosystem", y = "", title = "") +
  geom_tile(
    aes(
      x = comparison_1, y = taxon, fill = mean_contribution
    ),
    colour = "gray80"
  ) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", na.value = "gray", name = "Contribution") +  # Change gradient here
  #scale_x_discrete(expand = c(0, 0)) +  # Decrease spacing between x-axis categories
  theme_pubr() +
  theme(
    axis.text.x = element_text(
      hjust = 1, size = unit(24, "points"), angle = 90,
      vjust = 0.5, family = "Arial"
    ),
    axis.title.x = element_text(
      size = unit(25, "points"), face = "bold", family = "Arial"
    ),
    axis.line = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = unit(25, "points"), family = "Arial"),
    # legend.key.size = unit(1, "points"),
    legend.title = element_text(
      size = unit(30, "points"), face = "bold", family = "Arial"
    ),
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    plot.margin = margin(t = 6, r = 50, b = 2, l = 0),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    # axis.text.y = element_text(
    #   size = unit(10, "points"), family = "Arial"
    # ),
    axis.ticks.y = element_blank()
  )

# ggsave(
#   heatmap_contribution_candidate,
#   filename = "results/figures/heatmap_candidate.png",
#   width = 190, height = 200,
#   units = "mm", scale = 4, dpi = 500
# )

# heatmap_key_candidate <-
#   ggplot(data = subset(simper_result, microgroup != "Bonafide")) +
#   labs(x = "comparison_1", y = "comparison_2") +
#   #facet_grid(~ecosystem, scales = "free_x", space = "free_x", drop = TRUE) +
#   geom_tile(
#     aes(
#       x = habitat, y = taxon, alpha = factor(LIASP_isKeystone), fill = ecosystem
#     ),
#     colour = "gray80"
#   ) +
#   scale_fill_manual("Ecosystem", values = ecosystem_colors, guide = FALSE) +
#   scale_alpha_discrete(
#     "Keystones",
#     labels = c(
#       "1" = "It is Keystone", "0" = "It isn't Keystone", "NA" = "Absence"
#     ),
#     guide = FALSE
#   ) +
#   labs(x = "Habitats") +
#   theme(
#     axis.text.x = element_text(
#       hjust = 1,   size = unit(20, "points"), angle = 90,
#       vjust = 0.5, family = "Arial"
#     ),
#     axis.title.x = element_text(
#       size = unit(20, "points"), face = "bold", family = "Arial"
#     ),
#     axis.line = element_blank(),
#     legend.position = NULL,
#     legend.key = element_rect(fill = NA, colour = "black"),
#     legend.background = element_rect(
#       color = "black", fill = "white", size = 0.3, linetype = "solid"
#     ),
#     strip.placement = "outside",
#     panel.grid.major = element_blank(),
#     plot.margin = margin(t = 6, r = 1, b = 0.6, l = 3.8),
#     strip.background = element_rect(colour = "White", fill = "white"),
#     strip.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   )

## Bonafide ----------------------------
heatmap_contribution_bonafide <- ggplot(data = subset(simper_result_sum, microgroup == "Bonafide")) +
  labs(x = "Ecosystem", y = "Taxon", title = "") +
  geom_tile(
    aes(
      x = comparison_1, y = taxon, fill = mean_contribution
    ),
    colour = "gray80"
  ) +
  scale_fill_viridis_c(option = "cividis", name = "Contribution") +
  #scale_x_discrete(expand = c(0, 0)) +  # Decrease spacing between x-axis categories
  theme_pubr() +
  theme(
    axis.text.x = element_text(
      hjust = 1, size = unit(10, "points"), angle = 90,
      vjust = 0.5, family = "Arial"
    ),
    axis.title.x = element_text(
      size = unit(20, "points"), face = "bold", family = "Arial"
    ),
    axis.line = element_blank(),
    legend.position = NULL,
    legend.key = element_rect(fill = NA, colour = "black"),
    legend.background = element_rect(
      color = "black", fill = "white", size = 0.3, linetype = "solid"
    ),
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    plot.margin = margin(t = 6, r = 1, b = 0.6, l = 3.8),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    # axis.text.y = element_text(
    #   size = unit(10, "points"), family = "Arial"
    # ),
    axis.ticks.y = element_blank()
  )


########################## Holobiont VS free-living ############################
holo_vs_free_candidate <-
  subset(phyla_abundances, microgroup != "Bonafide") %>%
  group_by(taxon, life_style) %>%
  summarise(testando = mean(abundance)) %>%
  ggplot(aes(x = taxon, y = testando, fill = life_style)) +
  geom_col(position = "fill") +
  geom_hline(yintercept =  0.5) +
  coord_flip() +
  scale_y_continuous(
    breaks = c(0.1, 0.5, 0.9), labels = c("0%", "50%", "100%")
  ) +
  labs(y = NULL, x = NULL, fill = NULL) +
  scale_fill_manual(values = life_style_colors) +
  theme_void() +
  theme(
    axis.title = element_text(
      size = unit(20, "points"), face = "bold", family = "Arial"
    ),
    axis.text.x = element_text(size = unit(25, "points"), family = "Arial"),
    legend.position = "bottom",
    legend.title = element_text(
      size = unit(20, "points"), face = "bold", family = "Arial"
    ),
    legend.key.size = unit(20, "points"),
    legend.text = element_text(size = unit(20, "points"), family = "Arial"),
    plot.margin = margin(t = 6, r = 0, b = 2, l = 6)
  ) +
  labs(fill = "Life Style")

# Bonafide -----------------------------
holo_vs_free_bonafide <-
  subset(phyla_abundances, microgroup == "Bonafide") %>%
  group_by(taxon, life_style) %>%
  summarise(testando = mean(abundance)) %>%
  ggplot(aes(x = taxon, y = testando, fill = life_style)) +
  geom_col(position = "fill") +
  geom_hline(yintercept =  0.5) +
  coord_flip() +
  scale_y_continuous(
    breaks = c(0.1, 0.5, 0.9),
    labels = c("0%", "50%", "100%")
  ) +
  labs(y = NULL, x = NULL, fill = NULL) +
  scale_fill_manual(values = life_style_colors) +
  theme_void() +
  theme(
    axis.title = element_text(
      size = unit(20, "points"), face = "bold", family = "Arial"
    ),
    axis.text.x = element_text(size = unit(19, "points"), family = "Arial"),
    legend.position = "none",
    plot.margin = margin(5, 3.5, 1, 1)
  )

################################# Final plots ##################################
supplementary_dir <- "results/figures/"
if (!dir.exists(supplementary_dir)) {
  dir.create(supplementary_dir, recursive = TRUE)
}

# Get legends --------------------------

lifestyle_grob_plot <- ggplotGrob(holo_vs_free_candidate)
# Find which part of the grob layout contains the legend
lifestyle_legend_index <- which(lifestyle_grob_plot$layout$name == "guide-box-bottom")
# Extract the legend
lifestyle_legend <- lifestyle_grob_plot$grobs[[lifestyle_legend_index]]
holo_vs_free_candidate <- holo_vs_free_candidate +
  theme(legend.position = "none")

ecosystem_legend <- get_legend(jitter_candidate)
jitter_candidate <- jitter_candidate +
  theme(legend.position = "none")

# Panel 3 ------------------------------
panel_3 <- plot_grid(
  jitter_candidate,
  holo_vs_free_candidate,
  heatmap_contribution_candidate,
  rel_widths = c(0.7, 0.1, 0.3),
  ncol = 3, labels = c("A", "B", "C"),
  label_fontfamily = "Arial",
  label_x = c(0, -0.12, -0.017),
  label_size = 26,
  align = "h",
  axis = "l"
) + theme(plot.background = element_rect(fill = "white", colour = NA))

panel_3 <- ggdraw() +
  draw_plot(panel_3) +
  draw_plot(ecosystem_legend, x = .24, y = - .46, width = .004) +
  draw_plot(lifestyle_legend, x = .50, y = - .46, width = .004)

formats <- c("pdf", "svg")
for(format in formats){
  ggsave(
    panel_3, filename = paste0(supplementary_dir, "panel_3.", format),
    width = 190, height = 200,
    units = "mm", scale = 4, dpi = 500
  )
}

# Panel supplementary ------------------
supplementary <- plot_grid(
  jitter_bonafide,
  holo_vs_free_bonafide,
  heatmap_key_bonafide,
  rel_widths = c(1, 0.15, 1),
  ncol = 3, labels = c("a", "b", "c"),
  label_size = 26,
  label_fontfamily = "Arial",
  align = "h",
  axis = "bt",
  label_x = c(0, -.12, -.017)
) + theme(plot.background = element_rect(fill = "white", colour = NA))

supplementary <- ggdraw() +
  draw_plot(supplementary) +
  draw_plot(ecosystem_legend, x = .24, y = - .46, width = .004) +
  draw_plot(life_style_legend, x = .50, y = - .46, width = .004)

for (format in formats){
  ggsave(
    supplementary,
    filename = paste0(supplementary_dir, "bonafides.", format),
    width = 190, height = 110,
    units = "mm", scale = 4
  )
}