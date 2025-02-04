#' @title General prospect
#' @description Count samples, Analyses phyla structure through NMDS, Permanova,
#'   peranova, simper and gam. These analyses are mostly used in Panel 1.
#' @author Leonardo Brait and Camilo Ferreira
#' @date 2023-09-20

################################# Environment ##################################
set.seed(201094)
source("R/src/install_and_load.R")
install_and_load(libs = c(
  "lmPerm" = "2.1.0",
  "multcomp" = "1.4-25",
  "mgcv" = "1.9.0",
  "vegan" = "2.6-4",
  "tidyverse" = "2.0.0",
  "funrar" = "1.5.0"
))

################################# Load data ####################################
phyla_abundances <-
  read_csv("data/taxa_compositions/kraken_biomedb_relative_phyla.csv")
metadata <- read_csv("data/metadata/biome_classification.csv") %>%
  select(samples, life_style, ecosystem, habitat, latitude, longitude)
radiation <- read_csv("data/taxa_compositions/radiation_phyla.csv")

################################# Merge and Treat ##############################
treated_dir <- "data/treated/"
if (!dir.exists(treated_dir)) {
  dir.create(treated_dir)
}

# Merge --------------------------------
merged_df <- phyla_abundances %>%
  inner_join(metadata, by = "samples") %>%
  group_by(ecosystem, habitat) %>%
  mutate(total_samples = n()) %>%
  filter(total_samples >= 5) %>%
  select(-total_samples) %>%
  ungroup()
source("R/src/visualization_treatment.R")
merged_df <- treatment(merged_df)
write_csv(merged_df, paste0(treated_dir, "phyla_abundances_wide.csv"))

# Wide to long -------------------------
phyla_abundances_long <- merged_df %>%
  gather(
    taxon, abundance,
    -c(samples, ecosystem, habitat, life_style, latitude, longitude)
  )
################################# Countings ####################################
summary_dir <- "data/summaries/"
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir)
}

# Count samples ------------------------
total_samples <- merged_df %>%
  mutate(category = "total") %>%
  group_by(category) %>%
  summarise(total_samples = n())
lifestyle_samples <- merged_df %>%
  mutate(category = life_style) %>%
  group_by(life_style) %>%
  summarise(total_samples = n()) %>%
  rename(category = life_style)
ecosystem_samples <- merged_df %>%
  mutate(category = ecosystem) %>%
  group_by(ecosystem) %>%
  summarise(total_samples = n()) %>%
  rename(category = ecosystem)
habitats_samples <- merged_df %>%
  mutate(category = habitat) %>%
  group_by(habitat) %>%
  summarise(total_samples = n()) %>%
  rename(category = habitat)

# Merge and save -----------------------
samples_summary <- bind_rows(
  total_samples,      lifestyle_samples,
  ecosystem_samples,   habitats_samples,
)
write_csv(samples_summary, paste0(summary_dir, "samples_summary.csv"))

rm(
  total_samples,     lifestyle_samples,
  ecosystem_samples,  habitats_samples,
  samples_summary
)

############################ Communities Structures ############################
rdata_dir <- "data/rdata/"
if (!dir.exists(rdata_dir)) {
  dir.create(rdata_dir)
}
statistics_dir <- "data/statistics/"
if (!dir.exists(statistics_dir)) {
  dir.create(statistics_dir)
}

# NMDS ---------------------------------
standarized_abundances <- merged_df %>%
  select(-c(samples, latitude, longitude, ecosystem, life_style, habitat)) %>%
  decostand(method = "hellinger") %>%
  vegdist(method = "jaccard")
if (!file.exists(paste0(rdata_dir, "nmds.RData"))) {
  print("NMDS not found! Running...")
  nmds <- metaMDS(
    standarized_abundances,
    distance = "jaccard", try = 1, trymax = 4999, stress = 1, parallel = 60
  )
  save(nmds, file = paste0(rdata_dir, "nmds.RData"))
  print(paste0("NMDS saved in: ", rdata_dir))
} else {
  print("NMDS found! Jumping!")
}

# Permanova ----------------------------
## Process data
phyla_abundances_persite <- phyla_abundances_long %>%
  group_by(latitude, longitude, ecosystem, life_style, habitat, taxon) %>%
  summarise(abundance = mean(abundance)) %>%
  ungroup() %>%
  spread(taxon, abundance, fill = 0)
standarized_abundances_persite <- phyla_abundances_persite %>%
  select(-c(latitude, longitude, ecosystem, life_style, habitat)) %>%
  decostand(method = "hellinger")
standarized_abundances_persite_matrix <- standarized_abundances_persite %>%
  vegdist(method = "jaccard")

## Execute for lifestyle
if (!file.exists(paste0(rdata_dir, "permanova_lifestyle.RData"))) {
  print("Permanova for lifestyle not found! Running...")
  permanova_lifestyle <- adonis2(
    standarized_abundances_persite_matrix ~ phyla_abundances_persite$life_style,
    permutations = 4999, parallel = 30
  )
  save(
    permanova_lifestyle, file = paste0(rdata_dir, "permanova_lifestyle.RData")
  )
  capture.output(
    permanova_lifestyle,
    file = paste0(statistics_dir, "permanova_lifestyle.txt")
  )
  print(paste0("Permanova for lifestyle saved in: ", rdata_dir))
} else {
  print("Permanova for lifestyle found! Jumping...")
}

## Execute for ecosystem
if (!file.exists(paste0(rdata_dir, "permanova_ecosystem.RData"))) {
  print("Permanova for ecosystem not found! Running...")
  permanova_ecosystem <- adonis2(
    standarized_abundances_persite_matrix ~ phyla_abundances_persite$ecosystem,
    permutations = 4999, parallel = 30
  )
  save(
    permanova_ecosystem, file = paste0(rdata_dir, "permanova_ecosystem.RData")
  )
  capture.output(
    permanova_ecosystem,
    file = paste0(statistics_dir, "permanova_ecosystem.txt")
  )
  print(paste0("Permanova for ecosystem saved in: ", rdata_dir))
} else {
  print("Permanova for ecosystem found! Jumping!")
  load(paste0(rdata_dir, "permanova_ecosystem.RData"))
}

# Simper -------------------------------
if (!file.exists(paste0(rdata_dir, "simper.RData"))) {
  print("Simper not found!, running...")
  category <- phyla_abundances_persite[["ecosystem"]]
  raw_simper <-  simper(
    standarized_abundances_persite,
    group = category,
    parallel = 1,
    permutations = 4999
  )
  save(raw_simper, file = paste0(rdata_dir, "simper.RData"))
  print("Done!")
} else {
  print("Simper output already exists. Generating tables!")
  load(paste0(rdata_dir, "simper.RData"))
  source("R/simper_ranking.R")
}

############################## Attach Phyla Groups #############################
dpann_groups <- radiation %>%
  filter(microgroup == "DPANN") %>%
  pull(taxon)
cpr_groups <- radiation %>%
  filter(microgroup == "CPR") %>%
  pull(taxon)
phyla_abundances_long <- phyla_abundances_long %>%
  mutate(microgroup = case_when(
    taxon %in% dpann_groups ~ "DPANN",
    taxon %in% cpr_groups ~ "CPR",
    TRUE ~ "Bonafide"
  ))
write_csv(
  phyla_abundances_long,
  paste0(treated_dir, "phyla_abundances_long.csv")
)


# Process ------------------------------
phyla_abundances_persite_long <-
  phyla_abundances_long %>%
  group_by(
    taxon,       life_style,  ecosystem,
    habitat,       latitude,  longitude,
    microgroup
  ) %>%
  summarise(abundance = mean(abundance)) %>%
  ungroup()
write_csv(
  phyla_abundances_persite_long,
  paste0(treated_dir, "phyla_abundances_persite_long.csv")
)

###################### Summaries of phyla groups prevalences ###################
microgroups_totalabundance <- phyla_abundances_persite_long %>%
  group_by(microgroup) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  spread(microgroup, abundance) %>%
  as.matrix() %>%
  funrar:::make_relative() %>%
  as.data.frame() %>%
  gather(microgroup, abundance)
microgroups_counts <- phyla_abundances_persite_long %>%
  group_by(microgroup) %>%
  summarise(taxons = n_distinct(taxon)) %>%
  ungroup() %>%
  mutate(abundance = microgroups_totalabundance$abundance)

write_csv(
  microgroups_counts,
  paste0(summary_dir, "microgroups_general.csv")
)

microgroups_abundance <- phyla_abundances_persite_long %>%
  group_by(
    life_style, ecosystem,     habitat,
    latitude,   longitude,   microgroup
  ) %>%
  summarise(abundance = sum(abundance)) %>%
  complete(microgroup) %>%
  replace_na(list(abundance = 0)) %>%
  ungroup()
microgroups_richness <- phyla_abundances_persite_long %>%
  group_by(
    life_style, ecosystem, habitat,
    latitude,   longitude, microgroup
  ) %>%
  summarise(richness = sum(abundance > 0)) %>%
  complete(microgroup) %>%
  replace_na(list(richness = 0)) %>%
  ungroup()

# Merge and save -----------------------
microgroups_prevalence_persite <- microgroups_abundance %>%
  inner_join(
    microgroups_richness,
    by = c(
      "life_style", "ecosystem", "habitat",
      "latitude",   "longitude", "microgroup"
    )
  )
write_csv(
  microgroups_prevalence_persite,
  paste0(summary_dir, "microgroups_prevalence_persite.csv")
)

# Summarizations -----------------------
lifestyle_microgroups_prevalence <- microgroups_prevalence_persite %>%
  group_by(life_style, microgroup) %>%
  summarise(
    mean_abu_life = mean(abundance),
    sd_abu_life = sd(abundance),
    se_abu_life = sd_abu_life / sqrt(n()),
    mean_rich_life = mean(richness),
    sd_rich_life = sd(richness),
    se_rich_life = sd_rich_life / sqrt(n())
  ) %>%
  ungroup()
ecosystem_microgroups_prevalence <- microgroups_prevalence_persite %>%
  group_by(ecosystem, microgroup) %>%
  summarise(
    mean_abu_life = mean(abundance),
    sd_abu_life = sd(abundance),
    se_abu_life = sd_abu_life / sqrt(n()),
    mean_rich_life = mean(richness),
    sd_rich_life = sd(richness),
    se_rich_life = sd_rich_life / sqrt(n())
  ) %>%
  ungroup()

# save ---------------------------------
write_csv(
  lifestyle_microgroups_prevalence,
  paste0(summary_dir, "lifestyle_microgroups_prevalence.csv")
)
write_csv(
  ecosystem_microgroups_prevalence,
  paste0(summary_dir, "ecosystem_microgroups_prevalence.csv")
)

################################# GAM Analysis #################################

# Richness -----------------------------
bonafide_richness_gam <- gam(
  richness ~ s(latitude, bs = "tp", k = 20),
  data = microgroups_prevalence_persite %>% filter(microgroup == "Bonafide"),
  method = "REML"
)
capture.output(
  summary(bonafide_richness_gam),
  file = paste0(statistics_dir, "gam_general_bonafide_richness.txt")
)
cpr_richness_gam <- gam(
  richness ~ s(latitude, bs = "tp", k = 20),
  data = microgroups_prevalence_persite %>% filter(microgroup == "CPR"),
  method = "REML"
)
capture.output(
  summary(cpr_richness_gam),
  file = paste0(statistics_dir, "gam_general_cpr_richness.txt")
)
dpann_richness_gam <- gam(
  richness ~ s(latitude, bs = "tp", k = 20),
  data = microgroups_prevalence_persite %>% filter(microgroup == "DPANN"),
  method = "REML"
)
capture.output(
  summary(dpann_richness_gam),
  file = paste0(statistics_dir, "gam_general_dpann_richness.txt")
)

# Abundance ----------------------------
bonafide_abundance_gam <- gam(
  abundance ~ s(latitude, bs = "tp", k = 20),
  data = microgroups_prevalence_persite %>% filter(microgroup == "Bonafide"),
  method = "REML"
)
capture.output(
  summary(bonafide_abundance_gam),
  file = paste0(statistics_dir, "gam_general_bonafide_abundance.txt")
)
cpr_abundance_gam <- gam(
  abundance ~ s(latitude, bs = "tp", k = 20),
  data = microgroups_prevalence_persite %>% filter(microgroup == "CPR"),
  method = "REML"
)
capture.output(
  summary(cpr_abundance_gam),
  file = paste0(statistics_dir, "gam_general_cpr_abundance.txt")
)
dpann_abundance_gam <- gam(
  abundance ~ s(latitude, bs = "tp", k = 20),
  data = microgroups_prevalence_persite %>% filter(microgroup == "DPANN"),
  method = "REML"
)
capture.output(
  summary(dpann_abundance_gam),
  file = paste0(statistics_dir, "gam_general_dpann_abundance.txt")
)

############################# Peranova Analysis ################################
source("R/src/do_peranova.R")

for (i in c("Bonafide", "CPR", "DPANN")) {
  print(paste("Doing:", i, "general richness"))
  do_peranova(
    data = subset(microgroups_richness, microgroup == i),
    response = "richness",
    full_model = "life_style + ecosystem / habitat",
    predictor_letters = "ecosystem",
    file_name = paste0("data/statistics/peranova_general_richness_", i),
    method = "euclidian",
    parallel = 4,
    permutations = 4999
  )
}
for (i in c("Bonafide", "CPR", "DPANN")) {
  print(paste("Doing:", i, "general abundance"))
  do_peranova(
    data = subset(microgroups_abundance, microgroup == i),
    response = "abundance",
    full_model = "life_style + ecosystem / habitat",
    predictor_letters = "ecosystem",
    file_name = paste0("data/statistics/peranova_general_abundance_", i),
    method = "euclidian",
    parallel = 4,
    permutations = 4999
  )
}
