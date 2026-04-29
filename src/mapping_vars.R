# Setup -------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

# The Mercator projection -------------------------------------------------------------------------
crs_mercat <- "+proj=merc +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
worldmap_mercator <- ne_countries(scale = "large", returnclass = "sf")
graticules_mercator <- st_graticule(worldmap_mercator, lon = seq(-180, 180, 30), lat = seq(-90, 90, 30)) 

p_bigmap_mercator <- ggplot() +
  # geom_sf(data = graticules_mercator, color = "gray80", size = 0.3) +
  geom_sf(data = worldmap_mercator, fill = "white", color = "dimgray") +
  coord_sf() +
  labs(y = "Latitude", x = "Longitude")

# The Robinson projection -------------------------------------------------------------------------

crs_robin <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
worldmap_robinson <- worldmap_mercator %>% 
  st_transform(crs = crs_robin)
graticules_robinson <- st_graticule(worldmap_robinson, lon = seq(-180, 180, 30), lat = seq(-90, 90, 30)) 

p_bigmap_robinson <- ggplot() +
  # geom_sf(data = graticules_robinson, color = "gray80", size = 0.3) +
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

p_bigmap_mollweide <- ggplot() +
  # geom_sf(data = graticules_mollweide, color = "gray80", size = 0.3) +
  geom_sf(data = worldmap_mollweide, fill = "white", color = "dimgray") +
  coord_sf() +
  labs(y = "Latitude", x = "Longitude")

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
