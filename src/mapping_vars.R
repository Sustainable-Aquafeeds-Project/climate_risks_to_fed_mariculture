# Setup -------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

source("src/dirs.R")

# Modified rnaturalearth data ---------------------------------------------------------------------
my_world_path <- file.path(prepdata_path, "my_world.qs")

if (!file.exists(my_world_path)) {
  required_isos <- file.path(prepdata_path, "assigned_farms.qs") %>% 
    qs_read() %>% pull(ISO3_Code) %>% unique()
  
  # ne_countries() merges the French overseas departments (GUF, GLP, MTQ, MYT, REU) into the
  # France polygon and codes Palestine as PSX. Swap the France row for the separate French
  # map units, then build a single ISO3 column (iso_a3, falling back to adm0_a3 where it's -99).
  french_units <- ne_countries(scale = "large", type = "map_units", returnclass = "sf") %>% 
    filter(adm0_a3 == "FRA")
  
  my_world <- ne_countries(scale = "large", returnclass = "sf") %>% 
    filter(adm0_a3 != "FRA") %>% 
    bind_rows(french_units) %>% 
    mutate(ISO3_Code = case_when(
      adm0_a3 == "PSX" ~ "PSE",
      iso_a3 != "-99"  ~ iso_a3,
      .default = adm0_a3
    )) %>% 
    select(ISO3_Code, name, name_long, admin, continent, region_un, subregion, geometry)
  
  stopifnot(!any(duplicated(my_world$ISO3_Code)))
  stopifnot(length(setdiff(required_isos, my_world$ISO3_Code)) == 0)
  qs_save(my_world, my_world_path)
}

# Drop-in replacement for ne_countries(scale = "large", returnclass = "sf"); match on ISO3_Code
my_world <- qs_read(my_world_path)


# The Mercator projection -------------------------------------------------------------------------
crs_mercat <- "+proj=merc +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
worldmap_mercator <- my_world
graticules_mercator <- st_graticule(worldmap_mercator, lon = seq(-180, 180, 30), lat = seq(-90, 90, 30)) 

# ggplot() +
#   geom_sf(data = graticules_mercator, color = "gray80", size = 0.3) +
#   geom_sf(data = worldmap_mercator, fill = "white", color = "dimgray") +
#   coord_sf()

# The Robinson projection -------------------------------------------------------------------------
crs_robin <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
worldmap_robinson <- worldmap_mercator %>% 
  st_transform(crs = crs_robin)
graticules_robinson <- st_graticule(worldmap_robinson, lon = seq(-180, 180, 30), lat = seq(-90, 90, 30))

p_bigmap_robinson <- ggplot() +
  geom_sf(data = graticules_robinson, color = "gray80", size = 0.3) +
  geom_sf(data = worldmap_robinson, fill = "white", color = "dimgray") +
  coord_sf() +
  theme_void()

robinson_coord_sf <- function(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat)) {
  bbox <- st_bbox(
    c(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]),
    crs = 4326
  ) %>%
    st_as_sfc() %>%
    st_segmentize(dfMaxLength = 1) %>% 
    st_transform(crs_robin) %>%
    st_bbox()
  
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           crs = crs_robin)
}

# The Mollweide projection ------------------------------------------------------------------------
crs_moll <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

worldmap_mollweide <- worldmap_mercator %>% 
   st_transform(crs = crs_moll)
graticules_mollweide <- st_graticule(worldmap_mollweide, lon = seq(-180, 180, 30), lat = seq(-90, 90, 30)) 

# ggplot() +
#   geom_sf(data = graticules_mollweide, color = "gray80", size = 0.3) +
#   geom_sf(data = worldmap_mollweide, fill = "white", color = "dimgray") +
#   coord_sf()

mollweide_coord_sf <- function(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat)) {
  bbox <- st_bbox(
    c(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]),
    crs = 4326
  ) %>%
    st_as_sfc() %>%
    st_segmentize(dfMaxLength = 1) %>% 
    st_transform(crs_moll) %>%
    st_bbox()
  
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           crs = crs_moll)
}
