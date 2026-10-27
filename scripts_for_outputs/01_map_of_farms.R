# Map of assigned farm locations for each focus country.
# Assumes scripts_for_outputs/00_run_all.R has already been run (provides
# `assigned_farms`, `focus_countries`, `reports_dir`, and paths from src/dirs.R).

library(sf)
library(rnaturalearth)

here("src", "species_colours.R") %>% source()    # species_pal
here("src", "other_functions.R") %>% source()    # prettyplot()
here("src", "country_bboxes.R")  %>% source()    # bbox_lims()

# Dissolving the FAO-region pieces into one geometry per country fails under
# spherical (s2) geometry, so do the union planar and restore s2 afterwards.
sf_use_s2(FALSE)
EEZ_by_FAO_full <- file.path(prepdata_path, "EEZ_data_with_FAO_regions.qs") %>%
  qd_read() %>%
  st_transform(crs = 4326) %>%
  st_make_valid() %>% 
  group_by(ISO3_Code, country) %>% 
  summarise(.groups = "drop")
sf_use_s2(TRUE)

farm_location_map <- function(this_iso, out_dir) {
  these_farms  <- assigned_farms %>% filter(ISO3_Code == this_iso)
  this_country <- unique(these_farms$country)
  this_eez <- if (this_iso == "HKG") {
    EEZ_by_FAO_full %>% filter(ISO3_Code == "CHN")
  } else {
    EEZ_by_FAO_full %>% filter(ISO3_Code == this_iso)
  }

  # Countries straddling 180° (e.g. FJI) need their data shifted to 0..360 to match the bbox
  these_farms <- shift_for_bbox(these_farms, this_iso)
  this_eez    <- shift_for_bbox(this_eez, this_iso, union = TRUE)
  this_land   <- shift_for_bbox(land, this_iso)

  # Models present in this country, with a colour for each. Models not in `models_pal` (e.g. tuna) get a grey fallback.
  models_here <- sort(unique(these_farms$model_name))
  country_pal <- models_pal[models_here]
  names(country_pal) <- models_here
  country_pal[is.na(country_pal)] <- "grey50"

  # `eez_layer` is drawn between the land and the farms so the border doesn't cover points
  build_map <- function(eez_layer = NULL) {
    ggplot() +
      geom_sf(data = this_land, fill = "grey85", colour = "grey55", linewidth = 0.2) +
      eez_layer +
      geom_sf(data = these_farms, aes(colour = model_name), size = 2) +
      scale_colour_manual(values = country_pal, labels = models_lab, name = "Model") +
      do.call(coord_sf, bbox_lims(this_iso)) +
      labs(
        title = paste0("Assigned farms: ", country_info(this_iso)$name, " (", this_iso, ")"),
        x = NULL, y = NULL
      ) +
      prettyplot() +
      theme(legend.position = "bottom", legend.title = element_blank())
  }

  # ggsave(
  #   filename = file.path(out_dir, paste0("farm_location_map_", this_iso, ".png")),
  #   plot = build_map(), width = 8.5, height = 6.5, dpi = 300, bg = "white"
  # )

  # Create the same map but with a border around the country's EEZ
  p_eez <- build_map(
    geom_sf(data = this_eez, fill = NA, colour = "grey30", linewidth = 0.4, linetype = "dashed")
  )

  ggsave(
    filename = file.path(out_dir, paste0("farm_location_map_EEZ_", this_iso, ".png")),
    plot = p_eez, width = 8.5, height = 6.5, dpi = 300, bg = "white"
  )
}

walk2(focus_countries$ISO3_Code, focus_countries$dir_name, farm_location_map)

