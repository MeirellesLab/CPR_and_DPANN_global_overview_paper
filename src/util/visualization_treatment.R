#' @title Visualization Treatment
#' @description This script contains the functions used to treat the data for
#' visualization. Change names of habitats, ecosystems and life styles, then
#' reorder the factors according to holobiont first.
#' @param data Dataframe to be treated.
library(tidyverse)
library(dplyr)
library(forcats)
treatment <- function(data) {

  if ("life_style" %in% colnames(data)) {
    data <- data %>%
      mutate(
        life_style = case_when(
          life_style == "free-living" ~ "Free-living",
          life_style == "host-associated" ~ "Holobiont",
          TRUE ~ life_style
        ),
        life_style = fct_relevel(
          life_style,
          "Holobiont",
          "Free-living"
        )
      )
    data <- data
  }
  if ("ecosystem" %in% colnames(data)) {
    data <- data %>%
      mutate(
        ecosystem = case_when(
          ecosystem == "human_host-associated" ~ "Human Host",
          ecosystem == "animal_host-associated" ~ "Animal Host",
          ecosystem == "plant_associated" ~ "Plant Host",
          ecosystem == "groundwater" ~ "Groundwater",
          ecosystem == "freshwater" ~ "Freshwater",
          ecosystem == "wastewater" ~ "Wastewater",
          ecosystem == "saline_water" ~ "Saline Water",
          ecosystem == "sediment" ~ "Sediment",
          ecosystem == "soil" ~ "Soil",
          TRUE ~ ecosystem
        ),
        ecosystem = fct_relevel(
          ecosystem,
          "Human Host",
          "Animal Host",
          "Plant Host",
          "Groundwater",
          "Freshwater",
          "Wastewater",
          "Saline Water",
          "Sediment",
          "Soil"
        )
      )
    data <- data
  }
  if ("habitat" %in% colnames(data)) {
    data <- data %>%
      mutate(habitat = case_when(
        habitat == "agricultural_soil" ~ "Agricultural Soil",
        habitat == "wood_decomposition" ~ "Wood Decomposition",
        habitat == "animal_skin" ~ "Animal Skin",
        habitat == "alkaline_environment" ~ "Alkaline Environment",
        habitat == "alkaline_salt_lake" ~ "Alkaline Salt Lake",
        habitat == "anaerobic_sediment" ~ "Anaerobic Sediment",
        habitat == "animal_feces" ~ "Animal Feces",
        habitat == "animal_gut" ~ "Animal Gut",
        habitat == "animal_manure" ~ "Animal Manure",
        habitat == "anoxic_seawater" ~ "Anoxic Seawater",
        habitat == "aqueous_humour" ~ "Aqueous Humour",
        habitat == "beef" ~ "Beef",
        habitat == "bodily_fluid" ~ "Bodily Fluid",
        habitat == "bone" ~ "Bone",
        habitat == "sulfur_spring" ~ "Sulfur Spring",
        habitat == "fish" ~ "Fish",
        habitat == "mummy" ~ "Mummy",
        habitat == "chyme" ~ "Chyme",
        habitat == "coastal_marine_sediment" ~ "Coastal Marine Sediment",
        habitat == "coastal_seawater" ~ "Coastal Seawater",
        habitat == "compost" ~ "Compost",
        habitat == "contaminated_soil" ~ "Contaminated Soil",
        habitat == "contaminated_water" ~ "Contaminated Water",
        habitat == "continental_slope_sediment" ~ "Continental Slope Sediment",
        habitat == "coral" ~ "Coral",
        habitat == "coral_reef_biofilm" ~ "Coral Reef Biofilm",
        habitat == "coral_reef_seawater" ~ "Coral Reef Seawater",
        habitat == "desert_sediment" ~ "Desert Sediment",
        habitat == "desert_soil" ~ "Desert Soil",
        habitat == "drinking_water_reservoir" ~ "Drinking Water Reservoir",
        habitat == "estuarine_seawater" ~ "Estuarine Seawater",
        habitat == "farm_soil" ~ "Farm Soil",
        habitat == "forest_soil" ~ "Forest Soil",
        habitat == "freshwater_biofilm" ~ "Freshwater Biofilm",
        habitat == "garden_soil" ~ "Garden Soil",
        habitat == "grassland" ~ "Grassland",
        habitat == "hot_spring" ~ "Hot Spring",
        habitat == "human-gut" ~ "Human Gut",
        habitat == "human-skin" ~ "Human Skin",
        habitat == "human_associated_biofilm" ~ "Human Biofilm",
        habitat == "human_feces" ~ "Human Feces",
        habitat == "hydrothermal_vent" ~ "Hydrothermal Vent",
        habitat == "hypersaline_water" ~ "Hypersaline Water",
        habitat == "karst_porous" ~ "Karst Porous",
        habitat == "lake" ~ "Lake",
        habitat == "lake_sediment" ~ "Lake Sediment",
        habitat == "lung" ~ "Lung",
        habitat == "mangrove_sediment" ~ "Mangrove Sediment",
        habitat == "marine_biofilm" ~ "Marine Biofilm",
        habitat == "meadow_soil" ~ "Meadow Soil",
        habitat == "oceanic_seawater" ~ "Oceanic Seawater",
        habitat == "oil_affected_sediment" ~ "Oil Affected Sediment",
        habitat == "park_soil" ~ "Park Soil",
        habitat == "pasture_soil" ~ "Pasture Soil",
        habitat == "peat_soil" ~ "Peat Soil",
        habitat == "permafrost" ~ "Permafrost",
        habitat == "plant-associated" ~ "Plant Associated",
        habitat == "plant-associated_hypersaline_water" ~ "Hypersaline Dew",
        habitat == "polar_seawater" ~ "Polar Seawater",
        habitat == "prairie" ~ "Prairie",
        habitat == "reactor_sludge" ~ "Reactor Sludge",
        habitat == "reclaimed_water" ~ "Reclaimed Water",
        habitat == "rhizosphere" ~ "Rhizosphere",
        habitat == "river" ~ "River",
        habitat == "river_saline" ~ "River Saline",
        habitat == "river_sediment" ~ "River Sediment",
        habitat == "rumen" ~ "Rumen",
        habitat == "saline_evaporation_pond" ~ "Saline Evaporation Pond",
        habitat == "saline_marsh" ~ "Saline Marsh",
        habitat == "saline_marsh_sediment" ~ "Saline Marsh Sediment",
        habitat == "saline_water_biofilm" ~ "Saline Water Biofilm",
        habitat == "saliva" ~ "Saliva",
        habitat == "salt_lake" ~ "Salt Lake",
        habitat == "savanna_soil" ~ "Savanna Soil",
        habitat == "sea_worm" ~ "Sea Worm",
        habitat == "sewage" ~ "Sewage",
        habitat == "shrubland_soil" ~ "Shrubland",
        habitat == "sludge" ~ "Sludge",
        habitat == "sponge" ~ "Sponge",
        habitat == "stream" ~ "Stream",
        habitat == "stream_biofilm" ~ "Stream Biofilm",
        habitat == "street_sweepings" ~ "Street Sweepings",
        habitat == "tundra_soil" ~ "Tundra Soil",
        habitat == "wastewater_treatment_plant" ~ "Water Treatment Plant",
        habitat == "wood_fall" ~ "Wood Fall",
        habitat == "geyser" ~ "Geyser",
        habitat == "mine" ~ "Mine",
        habitat == "porous" ~ "Porous",
        habitat == "porous_contaminated" ~ "Porous Contaminated",
        habitat == "subsurface_saline" ~ "Saline",
        TRUE ~ habitat
      ))
    data <- data
  }
  return(data)

}