# Description ----------------------------------------------------------------------------------------------
# This document prepares the aquaculture locations from:
# Clawson, G., Kuempel, C. D., Frazier, M., Blasco, G., Cottrell, R. S., Froehlich, H. E., Metian, M., Nash, K. L., Többen, J., Verstaen, J., Williams, D. R., & Halpern, B. S. (2022). Mapping the spatial distribution of global mariculture production. Aquaculture, 553, 738066. https://doi.org/10.1016/j.aquaculture.2022.738066 and combines it with FAO fishing regions, ready to be paired with FAO production data.

library(tidyverse)
library(terra)
library(here)
library(purrr)
library(future)
library(furrr)

here("src", "dirs.R") %>% source()
here("src", "functions.R") %>% source()

# Raw data -------------------------------------------------------------------------------------------------
mari_locs <- file.path(rawdata_path, "all_marine_aquaculture_farms_sources_final.csv") %>% 
  read.csv()

mari_locs_sf <- mari_locs %>% 
  select(-c(data_type, data_type_2, details, data_source, data_year)) %>% 
  st_as_sf(coords = c("X", "Y"), crs = "EPSG:4326") %>% 
  mutate(farm_ID = row_number()) %>% 
  group_by(farm_ID) %>% 
  group_split()

# Get FAO fishing/country data -----------------------------------------------------------------------------
FAO_shp <- file.path(bigdata_path, "fao/FAO_AREAS_CWP_NOCOASTLINE", "FAO_AREAS_CWP_NOCOASTLINE.shp") %>% 
  read_sf() %>% 
  filter(F_LEVEL == "MAJOR") %>% 
  select(-c(NAME_FR, NAME_ES, F_AREA, F_SUBAREA, F_DIVISION, F_SUBDIVIS, F_SUBUNIT, F_LEVEL, OCEAN, SUBOCEAN, ID, F_STATUS))

# Just the paired codes and countries
FAO_code_key <- FAO_shp %>% 
  distinct(F_CODE, NAME_EN)
write.csv(FAO_code_key, file.path(prepdata_path, "FAO_code_key.csv"))

# Intersect farm points with FAO fishing code polygons ----------------------------------------------------
plan(
  strategy = "multisession", 
  workers = parallel::detectCores()-2
)
# this_farm <- mari_locs_sf[[1]]

farm_fao_intersect <- furrr::future_map(
  mari_locs_sf, 
  function(this_farm){st_intersection(this_farm, FAO_shp)
  },
  .progress = T,
  .options = furrr_options(seed = TRUE)
  ) %>% 
  bind_rows() %>% 
  mutate(
    species_group = as.factor(species_group),
    F_CODE =        as.integer(F_CODE),
    NAME_EN =       as.factor(NAME_EN)
  )

plan(sequential)

# Check that all farms have a FAO fishing region code and name
farm_fao_intersect %>% 
  filter(is.na(F_CODE) | is.na(NAME_EN))

# Save
qd_save(farm_fao_intersect, file.path(prepdata_path, "mariculture_locations_FAO_codes.qs"))
