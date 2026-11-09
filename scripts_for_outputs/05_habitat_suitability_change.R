# Maps of relative change in habitat suitability (HS) at each global warming
# level (GLW) vs the 2015-2024 reference period, for each focus country:
#   - one png per species model farmed in the country (with farm markers)
#   - one production-weighted png combining all species (no farm markers)
# Assumes scripts_for_outputs/00_run_all.R has already been run (provides
# `assigned_farms`, `land`, `focus_countries`, `GLW_ls`, `GLW_labels`,
# `report_theme()`, and paths from src/dirs.R).

library(sf)
library(terra)
library(tidyterra)
library(cmocean)
library(here)

here("src", "species_colours.R") %>% source()    # models_lab
here("src", "country_bboxes.R")  %>% source()    # bbox_lims(), shift_for_bbox(), crosses_antimeridian()

HS_dir    <- file.path(outs_path, "data", "habitat_suitability_decadal")
lyr_names <- c("ref_2015_2024", names(GLW_ls))

# Gap-fill one cell around the coastline so there are no white gaps between the
# land polygon and the HS data. This extrapolates slightly, so it is for plotting only.
fill_coastal_gaps <- function(r, coast_zone, w = 3) {
  if (is.null(coast_zone)) return(r)
  filled <- focal(r, w, "mean", na.rm = TRUE, na.policy = "only")
  cover(r, mask(filled, coast_zone))
}

# Shift a -180..180 raster to 0..360 (terra::rotate() only goes the other way)
shift_raster_360 <- function(r) {
  e <- ext(r)
  west <- crop(r, ext(-180, 0, e$ymin, e$ymax))
  east <- crop(r, ext(0, 180, e$ymin, e$ymax))
  ext(west) <- c(180, 360, e$ymin, e$ymax)
  merge(east, west)
}

# Mean across species rasters, weighted by each model's current production
weighted_HS <- function(rasts, weights) {
  w <- weights / sum(weights)
  Reduce(`+`, map2(rasts, w, \(r, wi) r * wi))
}

HS_change_theme <- function() {
  report_theme() +
    theme(
      legend.position     = "right",
      axis.text           = element_text(size = 7),
      legend.key.height   = unit(4.5, "lines"),
      legend.key.width    = unit(1.2, "lines"),
      legend.frame        = element_rect(colour = "black"),
      legend.ticks        = element_line(colour = "black"),
      legend.ticks.length = unit(0.3, "lines")
    )
}

# this_iso <- focus_countries$ISO3_Code[1]
# out_dir <- focus_countries$dir_name[1]
HS_change_maps <- function(this_iso, out_dir) {
  these_farms <- assigned_farms %>% filter(ISO3_Code == this_iso)
  models_here <- sort(unique(these_farms$model_name))

  # Species farmed here that have a complete set of HS rasters (ref + every GLW)
  HS_path   <- file.path(HS_dir, this_iso)
  HS_models <- keep(models_here, \(m) {
    all(file.exists(file.path(HS_path, m, paste0("MR_", lyr_names, ".tif"))))
  })
  if (length(HS_models) == 0) {
    warning(this_iso, ": no habitat suitability rasters, skipping")
    return(invisible(NULL))
  }
  no_HS <- setdiff(models_here, HS_models)
  if (length(no_HS) > 0) {
    message(this_iso, ": no complete HS rasters for ", paste(no_HS, collapse = ", "))
  }

  land_file  <- file.path(prepdata_path, "detailed_land_polys", paste0(this_iso, ".qs"))
  coast_file <- file.path(prepdata_path, "coastal_zone_polys",  paste0(this_iso, ".qs"))
  if (!file.exists(land_file) || !file.exists(coast_file)) {
    warning(this_iso, ": missing detailed land or coastal zone polygons, skipping")
    return(invisible(NULL))
  }

  # Countries straddling 180° (e.g. FJI) need everything shifted to 0..360 to match the bbox
  shift_raster <- crosses_antimeridian(this_iso)
  these_farms  <- shift_for_bbox(these_farms, this_iso)
  this_land    <- qs_read(land_file)
  coast_zone   <- qs_read(coast_file)

  # One SpatRaster per species: layers = ref + each GLW, scaled to 0-1
  HS_rasts <- map(
    set_names(HS_models),
    function(m) {
      r <- rast(file.path(HS_path, m, paste0("MR_", lyr_names, ".tif"))) / 365
      names(r) <- lyr_names
      if (shift_raster) r <- shift_raster_360(r)
      names(r) <- lyr_names
      r
    }
  )

  # Relative change at each GLW vs the reference period
  HS_change <- map(HS_rasts, \(r) (r[[names(GLW_ls)]] - r[["ref_2015_2024"]]) / r[["ref_2015_2024"]])

  HS_weights <- these_farms %>%
    st_drop_geometry() %>%
    distinct(model_name, value) %>%
    filter(model_name %in% HS_models)
  HS_weights <- set_names(HS_weights$value, HS_weights$model_name)[HS_models]
  HS_change_weighted <- weighted_HS(HS_change, HS_weights)

  # `farm_layer` is NULL for the weighted map, which has no farm markers
  plot_HS_change <- function(r, farm_layer = NULL) {
    lim <- max(abs(minmax(r)), na.rm = T)
    ggplot() +
      geom_spatraster(data = fill_coastal_gaps(r, coast_zone)) +
      geom_sf(data = this_land, fill = "grey85", colour = "grey55", linewidth = 0.2) +
      farm_layer +
      facet_wrap(~lyr, nrow = 2, labeller = as_labeller(GLW_labels)) +
      # Reversed so blue = increase, red = decrease
      scale_fill_gradientn(
        limits = c(-lim, lim), colours = rev(cmocean("balance")(256)),
        labels = scales::label_percent(), na.value = "transparent",
        name = "Change in\nhabitat\nsuitability"
      ) +
      do.call(coord_sf, bbox_lims(this_iso)) +
      labs(x = NULL, y = NULL) +
      HS_change_theme()
  }

  for (m in HS_models) {
    ggsave(
      filename = file.path(out_dir, paste0("HS_change_", m, "_", this_iso, ".png")),
      plot = plot_HS_change(
        HS_change[[m]],
        geom_sf(data = filter(these_farms, model_name == m), shape = 4, size = 1)
      ),
      width = 8.5, height = 6.5, dpi = 300, bg = "white"
    )
  }

  ggsave(
    filename = file.path(out_dir, paste0("HS_change_weighted_", this_iso, ".png")),
    plot = plot_HS_change(HS_change_weighted),
    width = 8.5, height = 6.5, dpi = 300, bg = "white"
  )
}

walk2(focus_countries$ISO3_Code, focus_countries$dir_name, HS_change_maps)
