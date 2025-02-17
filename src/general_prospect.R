#' @title General prospect
#' @description Count samples, Analyses phyla structure through NMDS, Permanova,
#'   peranova, simper and gam. These analyses are mostly used in Panel 1.
#' @author Leonardo Brait and Camilo Ferreira
#' @date 2023-09-20

################################# Environment ##################################
set.seed(201094)
print("Loading libraries...")
source("src/util/install_and_load.R")
install_and_load(libs = c(
  "lmPerm" = "2.1.0",
  "multcomp" = "1.4-25",
  "mgcv" = "1.9.0",
  "vegan" = "2.6-4",
  "tidyverse" = "2.0.0",
  "funrar" = "1.5.0"
))

################################# Load data ####################################
print("Loading data...")
phyla_abundances <-
  read_csv("data/taxa_compositions/kraken_biomedb_relative_phyla.csv")

metadata <- read_csv("data/metadata/biome_classification.csv") %>%
  dplyr::select(samples, life_style, ecosystem, habitat, latitude, longitude)
  
radiation <- read_csv("data/taxa_compositions/radiation_phyla.csv")

################################# Merge and Treat ##############################
treated_dir <- "data/treated/"
if (!dir.exists(treated_dir)) {
  dir.create(treated_dir)
}

# Merge --------------------------------
print("Merging data...")
merged_df <- phyla_abundances %>%
  inner_join(metadata, by = "samples") %>%
  group_by(ecosystem, habitat) %>%
  mutate(total_samples = n()) %>%
  filter(total_samples >= 5) %>%
  select(-total_samples) %>%
  ungroup()
source("src/util/visualization_treatment.R")
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
print("Counting samples...")
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

# Main MDS ---------------------------------
print("Running general NMDS...")
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
print("Running general Permanova...")
## Process data
phyla_abundances_persite <- phyla_abundances_long %>%
  group_by(latitude, longitude, ecosystem, life_style, habitat, taxon) %>%
  summarise(abundance = mean(abundance)) %>%
  ungroup() %>%
  spread(taxon, abundance, fill = 0)
standarized_abundances_persite <- phyla_abundances_persite %>%
  dplyr::select(-c(latitude, longitude, ecosystem, life_style, habitat)) %>%
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
print("Running Simper ecosystems...")
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
 source("src/util/simper_ranking.R")
}

print("Running Simper habitats...")
if (!file.exists(paste0(rdata_dir, "simper_habitats.RData"))) {
 print("Simper habitats not found!, running...")
 category <- phyla_abundances_persite[["habitat"]]
 raw_simper_habitats <- simper(
   standarized_abundances_persite,
   group = category,
   parallel = 1,
   permutations = 4999
 )
 save(raw_simper_habitats, file = paste0(rdata_dir, "simper_habitats.RData"))
 print("Done!")
} else {
 print("Simper habitats output already exists. Generating tables!")
 load(paste0(rdata_dir, "simper_habitats.RData"))
 source("src/util/simper_ranking_habitats.R")
}

############################## Attach Phyla Groups #############################
print("Attaching phyla groups...")
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
print("Summarizing phyla groups prevalences...")
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
################################# NMDS by microgroup ###########################
print("Running NMDS by microgroup...")

# Creating abundance table by microgroup (the taxa are in the columns, filter them by the variables cpr_groups and dpann_groups, and the rest id bonafide)
print("Creating abundance table by microgroup...")
# Remember that the taxa are the columns 2 to 160
merged_df_cpr <- merged_df %>%
  dplyr::select(one_of(c("samples", "ecosystem", "habitat", "life_style", "latitude", "longitude"), cpr_groups))

merged_df_dpann <- merged_df %>%
  dplyr::select(one_of(c("samples", "ecosystem", "habitat", "life_style", "latitude", "longitude"), dpann_groups))

all_unculturable <- c(cpr_groups, dpann_groups)
# Selecione as colunas desejadas, excluindo as colunas em exclude_groups
merged_df_bonafide <- merged_df %>%
  dplyr::select(setdiff(colnames(merged_df), all_unculturable))

#Chek if there is any row with all zeros in each df 
print("Checking if there is any row with all zeros in each df...")
print("CPR")
print(sum(apply(merged_df_cpr[, -c(1, 2, 3, 4, 5, 6)], 1, sum) == 0))

print("DPANN")
print(sum(apply(merged_df_dpann[, -c(1, 2, 3, 4, 5, 6)], 1, sum) == 0))

print("bonafide")
print(sum(apply(merged_df_bonafide[, -c(1, 45, 46, 47, 48, 49)], 1, sum) == 0))

print("Only DPANN df has rows with all zeros, so we will remove them")
merged_df_dpann <- merged_df_dpann[apply(merged_df_dpann[, -c(1, 2, 3, 4, 5, 6)], 1, sum) != 0, ]

# Bonafide NMDS -----------------------
print("Running NMDS Bonafide...")

standarized_abundances_bonafide <- merged_df_bonafide %>%
  dplyr::select(-c(samples, latitude, longitude, ecosystem, life_style, habitat)) %>%
  decostand(method = "hellinger") %>%
  vegdist(method = "jaccard")
if (!file.exists(paste0(rdata_dir, "nmds_bonafide.RData"))) {
  print("NMDS Bonafide not found! Running...")
  nmds_bonafide <- metaMDS(
    standarized_abundances_bonafide,
    distance = "jaccard", try = 1, trymax = 4999, stress = 1, parallel = 60
  )
  save(nmds_bonafide, file = paste0(rdata_dir, "nmds_bonafide.RData"))
  print(paste0("NMDS Bonafide saved in: ", rdata_dir))
} else {
  print("NMDS Bonafide found! Jumping!")
}

# CPR NMDS ---------------------------
print("Running NMDS CPR...")
standarized_abundances_cpr <- merged_df_cpr %>%
  dplyr::select(-c(samples, latitude, longitude, ecosystem, life_style, habitat)) %>%
  decostand(method = "hellinger") %>%
  vegdist(method = "jaccard")
if (!file.exists(paste0(rdata_dir, "nmds_cpr.RData"))) {
  print("NMDS CPR not found! Running...")
  nmds_cpr <- metaMDS(
    standarized_abundances_cpr,
    distance = "jaccard", try = 1, trymax = 4999, stress = 1, parallel = 60
  )
  save(nmds_cpr, file = paste0(rdata_dir, "nmds_cpr.RData"))
  print(paste0("NMDS CPR saved in: ", rdata_dir))
} else {
  print("NMDS CPR found! Jumping!")
}

# DPANN NMDS -------------------------
print("Running NMDS DPANN...")
standarized_abundances_dpann <- merged_df_dpann %>%
  dplyr::select(-c(samples, latitude, longitude, ecosystem, life_style, habitat)) %>%
  decostand(method = "hellinger") %>%
  vegdist(method = "jaccard")
if (!file.exists(paste0(rdata_dir, "nmds_dpann.RData"))) {
  print("NMDS DPANN not found! Running...")
  nmds_dpann <- metaMDS(
    standarized_abundances_dpann,
    distance = "jaccard", try = 1, trymax = 4999, stress = 1, parallel = 60
  )
  save(nmds_dpann, file = paste0(rdata_dir, "nmds_dpann.RData"))
  print(paste0("NMDS DPANN saved in: ", rdata_dir))
} else {
  print("NMDS DPANN found! Jumping!")
}

################################# Permanova by microgroup #######################
print("Running Permanova by microgroup...")
# Bonafide Permanova -----------------
print("Running Permanova Bonafide...")
# Process ------------------------------
phyla_abundances_persite_long_bonafide <- phyla_abundances_persite_long %>%
  filter(microgroup == "Bonafide") %>%
  dplyr::select(-c(microgroup)) %>%
  group_by(latitude, longitude, ecosystem, life_style, habitat, taxon) %>%
  dplyr::summarise(abundance = mean(abundance)) %>%
  ungroup() %>%
  spread(taxon, abundance, fill = 0)
standarized_abundances_persite_matrix_bonafide <- phyla_abundances_persite_long_bonafide %>%
  dplyr::select(-c(latitude, longitude, ecosystem, life_style, habitat)) %>%
  decostand(method = "hellinger") %>%
  vegdist(method = "jaccard")

# Execute for lifestyle -----------------------------
if (!file.exists(paste0(rdata_dir, "permanova_lifestyle_bonafide.RData"))) {
  print("Permanova for lifestyle Bonafide not found! Running...")
  permanova_lifestyle_bonafide <- adonis2(
    standarized_abundances_persite_matrix_bonafide ~ phyla_abundances_persite_long_bonafide$life_style,
    permutations = 4999, parallel = 30
  )
  save(
    permanova_lifestyle_bonafide,
    file = paste0(rdata_dir, "permanova_lifestyle_bonafide.RData")
  )
  capture.output(
    permanova_lifestyle_bonafide,
    file = paste0(statistics_dir, "permanova_lifestyle_bonafide.txt")
  )
  print(paste0("Permanova for lifestyle Bonafide saved in: ", rdata_dir))
} else {
  print("Permanova for lifestyle Bonafide found! Jumping...")
}

# Execute for ecosystem -----------------------------
if (!file.exists(paste0(rdata_dir, "permanova_ecosystem_bonafide.RData"))) {
  print("Permanova for ecosystem Bonafide not found! Running...")
  permanova_ecosystem_bonafide <- adonis2(
    standarized_abundances_persite_matrix_bonafide ~ phyla_abundances_persite_long_bonafide$ecosystem,
    permutations = 4999, parallel = 30
  )
  save(
    permanova_ecosystem_bonafide,
    file = paste0(rdata_dir, "permanova_ecosystem_bonafide.RData")
  )
  capture.output(
    permanova_ecosystem_bonafide,
    file = paste0(statistics_dir, "permanova_ecosystem_bonafide.txt")
  )
  print(paste0("Permanova for ecosystem Bonafide saved in: ", rdata_dir))
} else {
  print("Permanova for ecosystem Bonafide found! Jumping...")
}

# CPR Permanova ------------------------
print("Running Permanova CPR...")
# Process ------------------------------
phyla_abundances_persite_long_cpr <-
  phyla_abundances_persite_long %>%
  filter(microgroup == "CPR")%>%
  dplyr::select(-c(microgroup)) %>%
  group_by(latitude, longitude, ecosystem, life_style, habitat, taxon) %>%
  dplyr::summarise(abundance = mean(abundance)) %>%
  ungroup() %>%
  spread(taxon, abundance, fill = 0)
standarized_abundances_persite_matrix_cpr <-
  phyla_abundances_persite_long_cpr %>%
  dplyr::select(-c(latitude, longitude, ecosystem, life_style, habitat)) %>%
  decostand(method = "hellinger") %>%
  vegdist(method = "jaccard")

# Execute for lifestyle -----------------------------
if (!file.exists(paste0(rdata_dir, "permanova_lifestyle_cpr.RData"))) {
  print("Permanova for lifestyle CPR not found! Running...")
  permanova_lifestyle_cpr <- adonis2(
    standarized_abundances_persite_matrix_cpr ~ phyla_abundances_persite_long_cpr$life_style,
    permutations = 4999, parallel = 30
  )
  save(
    permanova_lifestyle_cpr,
    file = paste0(rdata_dir, "permanova_lifestyle_cpr.RData")
  )
  capture.output(
    permanova_lifestyle_cpr,
    file = paste0(statistics_dir, "permanova_lifestyle_cpr.txt")
  )
  print(paste0("Permanova for lifestyle CPR saved in: ", rdata_dir))
} else {
  print("Permanova for lifestyle CPR found! Jumping...")
}

# Execute for ecosystem -----------------------------
if (!file.exists(paste0(rdata_dir, "permanova_ecosystem_cpr.RData"))) {
  print("Permanova for ecosystem CPR not found! Running...")
  permanova_ecosystem_cpr <- adonis2(
    standarized_abundances_persite_matrix_cpr ~ phyla_abundances_persite_long_cpr$ecosystem,
    permutations = 4999, parallel = 30
  )
  save(
    permanova_ecosystem_cpr,
    file = paste0(rdata_dir, "permanova_ecosystem_cpr.RData")
  )
  capture.output(
    permanova_ecosystem_cpr,
    file = paste0(statistics_dir, "permanova_ecosystem_cpr.txt")
  )
  print(paste0("Permanova for ecosystem CPR saved in: ", rdata_dir))
} else {
  print("Permanova for ecosystem CPR found! Jumping...")
}

# DPANN Permanova -----------------------
print("Running Permanova DPANN...")
# Process ------------------------------
phyla_abundances_persite_long_dpann <-
  phyla_abundances_persite_long %>%
  filter(microgroup == "DPANN")%>%
  dplyr::select(-c(microgroup)) %>%
  group_by(latitude, longitude, ecosystem, life_style, habitat, taxon) %>%
  dplyr::summarise(abundance = mean(abundance)) %>%
  ungroup() %>%
  spread(taxon, abundance, fill = 0)

#Chek if there is any row with all zeros in phyla_abundances_persite_long_dpann
print("Checking if there is any row with all zeros in phyla_abundances_persite_long_dpann...")
print(sum(apply(phyla_abundances_persite_long_dpann[, -c(1, 2, 3, 4, 5, 6)], 1, sum) == 0))

print("Had 3 zeros rows, so we will remove them")
phyla_abundances_persite_long_dpann <- phyla_abundances_persite_long_dpann[apply(phyla_abundances_persite_long_dpann[, -c(1, 2, 3, 4, 5, 6)], 1, sum) != 0, ]

standarized_abundances_persite_matrix_dpann <-
  phyla_abundances_persite_long_dpann %>%
  dplyr::select(-c(latitude, longitude, ecosystem, life_style, habitat)) %>%
  decostand(method = "hellinger") %>%
  vegdist(method = "jaccard")

# Execute for lifestyle -----------------------------
if (!file.exists(paste0(rdata_dir, "permanova_lifestyle_dpann.RData"))) {
  print("Permanova for lifestyle DPANN not found! Running...")
  permanova_lifestyle_dpann <- adonis2(
    standarized_abundances_persite_matrix_dpann ~ phyla_abundances_persite_long_dpann$life_style,
    permutations = 4999, parallel = 30
  )
  save(
    permanova_lifestyle_dpann,
    file = paste0(rdata_dir, "permanova_lifestyle_dpann.RData")
  )
  capture.output(
    permanova_lifestyle_dpann,
    file = paste0(statistics_dir, "permanova_lifestyle_dpann.txt")
  )
  print(paste0("Permanova for lifestyle DPANN saved in: ", rdata_dir))
} else {
  print("Permanova for lifestyle DPANN found! Jumping...")
}

# Execute for ecosystem -----------------------------
if (!file.exists(paste0(rdata_dir, "permanova_ecosystem_dpann.RData"))) {
  print("Permanova for ecosystem DPANN not found! Running...")
  permanova_ecosystem_dpann <- adonis2(
    standarized_abundances_persite_matrix_dpann ~ phyla_abundances_persite_long_dpann$ecosystem,
    permutations = 4999, parallel = 30
  )
  save(
    permanova_ecosystem_dpann,
    file = paste0(rdata_dir, "permanova_ecosystem_dpann.RData")
  )
  capture.output(
    permanova_ecosystem_dpann,
    file = paste0(statistics_dir, "permanova_ecosystem_dpann.txt")
  )
  print(paste0("Permanova for ecosystem DPANN saved in: ", rdata_dir))
} else {
  print("Permanova for ecosystem DPANN found! Jumping...")
}

################################# GAM Analysis #################################
print("Running GAM analysis...")
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
print("Running Peranova analysis...")
source("src/util/do_peranova.R")

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

print("General prospect done!")
