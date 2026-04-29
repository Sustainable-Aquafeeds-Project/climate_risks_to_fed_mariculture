# Description ----------------------------------------------------------------------------------------------
# This script downloads and prepares Exclusive Economic Zone (EEZ) data from the Marine Regions
# EEZ land union dataset (v4, October 2024). It:
#   1. Loads EEZ polygons, filtering to uncontested "Union EEZ and country" areas
#   2. Applies a country-code fix for Mayotte (miscoded as Comoros in the source data)
#   3. Subsets to only the countries/regions required for the analysis
#   4. Dissolves polygons into a single multipolygon per country
#   5. Creates a 50km coastal buffer and clips EEZ data accordingly
#   6. Intersects EEZ areas with FAO major fishing regions and removes land areas
#   7. Saves the final EEZ-by-FAO-region dataset for downstream use

library(tidyverse)
library(sf)

here("src", "dirs.R") %>% source()
here("src", "functions.R") %>% source()

# Get EEZ data for each country ----------------------------------------------------------------------------
EEZ_data <- file.path(bigdata_path, "marine-regions", "EEZ_land_union_v4_202410", "EEZ_land_union_v4_202410.shp") %>% 
  read_sf() %>% 
  filter(POL_TYPE %in% c("Union EEZ and country")) %>% # filter out areas with overlapping/conflicting claims
  select(geometry, ISO_TER1) %>% 
  filter(!is.na(ISO_TER1)) %>% 
  st_transform(crs = "+proj=moll")

# Mayotte is coded as Comoros - fix
EEZ_data <- rbind(
  EEZ_data,
  EEZ_data %>% filter(ISO_TER1 == "COM") %>% mutate(ISO_TER1 = "MYT")
) %>% 
  rename(ISO3_Code = ISO_TER1)

# # Narrow down to only the regions required
# EEZ_data <- left_join(
#   regions_required_3, 
#   EEZ_data, 
#   by = c("ISO3_Code" = "ISO_TER1"), 
#   relationship = "many-to-many"
#   ) %>% 
#   st_as_sf()

# Dissolve into a single multipolygon per country - the data says it splits by FAO fishing zone but it LIES!
EEZ_data <- EEZ_data %>%
  group_by(ISO3_Code) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# Load pre-saved land polygons -----------------------------------------------------------------------------
# The code below was used to generate the saved land object (takes a long time due to st_snap/st_make_valid)
# land <- ne_countries(scale = "large", returnclass = "sf") %>%
#   st_transform(crs = "+proj=moll") %>%
#   st_union() %>% 
#   st_snap(., ., tolerance = 1e-9) %>% 
#   st_make_valid()
# qd_save(land, file.path(data_path, "intermediates", "land.qs"))

land <- qd_read(file.path(data_path, "intermediates", "land.qs"))

EEZ_data <- EEZ_data %>% 
  st_make_valid()

# Check it looks ok
ggplot(EEZ_data) +
  geom_sf(aes(geometry = geometry, colour = ISO3_Code)) +
  coord_sf() + 
  theme_void() + theme(legend.position = "none")

# Restrict EEZ areas to within 30km of the coast -----------------------------------------------------------
coast <- land %>%
  st_buffer(dist = 30000)

EEZ_coast <- EEZ_data %>% 
  st_difference(coast)

# Check it looks ok
ggplot(EEZ_coast) +
  geom_sf(aes(geometry = geometry, fill = ISO3_Code)) +
  coord_sf() + 
  theme_void() + theme(legend.position = "none")

# Split EEZ data by FAO fishing region ---------------------------------------------------------------------
FAO_shp <- file.path(bigdata_path, "fao/FAO_AREAS_CWP_NOCOASTLINE", "FAO_AREAS_CWP_NOCOASTLINE.shp") %>% 
  read_sf() %>% 
  filter(F_LEVEL == "MAJOR") %>% 
  select(F_CODE, geometry) %>% 
  st_transform(crs = "+proj=moll")

EEZ_by_FAO <- st_intersection(EEZ_dissolved, FAO_shp) %>%
  st_make_valid() %>%
  st_difference(land)

# Check it looks ok
ggplot(filter(EEZ_by_FAO, ISO3_Code == "AUS")) +
  geom_sf(aes(geometry = geometry, fill = F_CODE)) +
  coord_sf() + 
  theme_void()

# Save
qd_save(EEZ_by_FAO, file.path(prepdata_path, "EEZ_data_with_FAO_regions.qs"))
