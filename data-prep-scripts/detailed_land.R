# The purpose of this script is to take detailed polygons of coastlines and save them into more manageable files for each country.
# Unlike most of the data prep scripts, this one was conceived of quite late in the process so requires that many other variables/scripts already exist.

library(sf)
library(here)
library(tidyverse)
library(qs2)
library(furrr)
library(future)

here("src/dirs.R") %>% source()
here("src/country_bboxes.R") %>% source()

land_path <- file.path(prepdata_path, "detailed_land_polys")
coast_path <- file.path(prepdata_path, "coastal_zone_polys")
dir.create(land_path, showWarnings = F)
dir.create(coast_path, showWarnings = F)

assigned_farms_final <- file.path(prepdata_path, "assigned_farms_final.qs") %>% qs_read()
isos <- assigned_farms_final %>% pull(ISO3_Code) %>% unique() %>% sort()

land_files <- file.path(land_path, paste0(isos, ".qs"))
coast_files <- file.path(coast_path, paste0(isos, ".qs"))

# Only create files that don't already exist (delete the files to regenerate)
land_isos <- isos[!file.exists(land_files)]
coast_isos <- isos[!file.exists(coast_files)]

# The OSM split land polygons are planar-valid but not always valid on the sphere,
# so s2 throws "Loop 0 is not valid: Edge crosses edge" for some countries.
sf_use_s2(FALSE)

# 1. Detailed land polygons ----
walk(
  land_isos,
  function(this_iso) {
    this_bbox <- bbox_lims(this_iso)
    this_bbox <- st_bbox(c(
      xmin = this_bbox[['xlim']][1], xmax = this_bbox[['xlim']][2], 
      ymin = this_bbox[['ylim']][1], ymax = this_bbox[['ylim']][2]
    ), crs = st_crs(4326))

    this_land <- st_read(
      file.path(rawdata_path, "land-polygons-split-4326", "land_polygons.shp"),
      wkt_filter = st_as_text(st_as_sfc(this_bbox)),
      quiet = T
    )
    this_land <- st_make_valid(this_land)
    # crop before union so the union only works on what we keep
    this_land <- st_crop(this_land, this_bbox)
    this_land <- st_union(this_land)          # merge the tiles the file is split into
    this_land <- st_make_valid(this_land)

    qs_save(this_land, file.path(land_path, paste0(this_iso, ".qs")))
  },
  .progress = T
)

# 2. Coastal zone polygons ----
# Buffer around the land polygons (in degrees, since s2 is off), read from the files saved above
plan(multisession, workers = 15)
future_walk(
  coast_isos,
  function(this_iso) {
    this_land <- file.path(land_path, paste0(this_iso, ".qs")) %>% qs_read()

    this_coast <- this_land %>% 
      st_simplify(dTolerance = 0.075, preserveTopology = TRUE) %>%   # ~10% of the buffer
      st_buffer(0.75, nQuadSegs = 5)  # This is 0.75 degrees
    this_coast <- st_make_valid(this_coast)

    qs_save(this_coast, file.path(coast_path, paste0(this_iso, ".qs")))
  },
  .progress = T
)
sf_use_s2(TRUE)
plan(sequential)
