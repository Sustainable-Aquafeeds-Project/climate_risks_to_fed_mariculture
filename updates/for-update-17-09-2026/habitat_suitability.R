library(tidyverse)
library(here)
library(qs2)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthhires)
library(metR)

source(here("src/dirs.R"))
source(here("src/other_functions.R"))
source(here("src/model_functions.R"))

production_files <- list.files(file.path(outs_path, "data/farm_production"), full.names = T, pattern = "MC5000.qs")
farms_done <- as.integer(str_extract(basename(production_files), "(?<=farmID_)\\d+"))

assigned_farms_ITL <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs_read() %>% 
  filter(country == "Italy") %>% 
  mutate(done = case_when(farm_ID %in% farms_done ~ T, T ~ F))

# Feed and species params ----------------------------------------------------------------
species_ITL <- unique(assigned_farms_ITL$model_name)

species_params <- file.path(prepdata_path, "species_parameters") %>% 
  list.files(full.names = T) %>% 
  str_subset(pattern = paste0(species_ITL, collapse = "|"))

species_params_nms <- species_params %>% basename() %>% str_remove_all(".qs|species_parameters_")

species_params <- map(species_params, qd_read) %>% 
  setNames(species_params_nms)

feed_params <- file.path(prepdata_path, "feed_parameters") %>% 
  list.files(full.names = T)

feed_names <- feed_params %>% 
  basename() %>% 
  str_remove_all(".qs")

feed_params <- map(feed_params, qs_read) %>% 
  setNames(feed_names)

feed_params <- list(
  "japanese_amberjack" = feed_params[["japanese_amberjack-japan"]],
  "milkfish"           = feed_params[["milkfish-philippines"]],
  "european_seabass"   = feed_params[["european_seabass-europe"]],
  "seabream"           = feed_params[["seabream-europe"]],
  "tiger_prawn"        = feed_params[["tiger_prawn-world"]],
  "whiteleg_shrimp"    = feed_params[["whiteleg_shrimp-china"]]
)

# Feeding response in the EEZ ------------------------------------------------------------
tos_tifs <- list.files(here("for-update-14-07-26/tos_data"), full.names = T) %>% sort()
years <- str_extract(tos_tifs, "(?<=_)\\d{4}(?=_\\d{3}\\.tif$)")
tos_tifs <- split(tos_tifs, years)

# yr_tifs <- tos_tifs[[1]]
feeding_response <- map(
  tos_tifs, 
  function(yr_tifs) {
    tos_rast <- rast(yr_tifs)
    resp_rast <- get_species_responses_r(
      temp_stack = tos_rast,
      species_params = species_params[["seabream"]],
      feed_params = feed_params[["seabream"]],
      ref_weight = unname(species_params[["seabream"]]["CS"])/2,
      out_cols = c("feeding_response")
    )[[1]]
    app(resp_rast, "sum")
  },
  .progress = T
) %>% setNames(names(tos_tifs))

decades <- list(2020:2029, 2030:2039, 2040:2049, 2050:2059, 2060:2069, 2070:2079, 2080:2089, 2090:2099)


land <- ne_countries(scale = "large", returnclass = "sf")
brks <- seq(0, 375, by = 75)

df <- as.data.frame(feeding_response[[1]], xy = TRUE, na.rm = FALSE) %>% 
  rename(lon = x, lat = y) 

ITA_bbox <- bbox("ITA")

df %>% 
  ggplot(aes(lon, lat, z = sum)) +
  geom_contour_fill(breaks = brks, na.fill = FALSE) +
  geom_contour2(breaks = brks, colour = "grey30", linewidth = 0.15) +
  geom_text_contour(breaks = brks, stroke = 0.15, size = 2.5, skip = 1) +
  geom_sf(data = land, fill = "grey80", inherit.aes = FALSE) +
  scale_fill_divergent(name = "Feeding response") +   # continuous, centred on 0
  # scale_fill_divergent_discretised() +
  coord_sf(xlim = c(ITA_bbox['xmin'], ITA_bbox['xmax']), ylim = c(ITA_bbox['ymin'], ITA_bbox['ymax']), expand = FALSE)
