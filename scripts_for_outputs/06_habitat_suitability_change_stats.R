# Summary statistics of habitat suitability (HS) for each focus country: mean HS
# in the reference period (2015-2024) and at each global warming level (GLW), plus
# the % change from the reference period, per species model and production-weighted.
# Two zones per country, one csv per country (for copying into report tables):
#   - "EEZ":          every cell of the country's HS rasters
#   - "coastal_50km": sea cells within ~50 km of the country's coastline. Cells are kept
#                     if their centre falls inside a 62.5 km (50 + 25/2) buffer around the
#                     country's land polygon, so cells near other countries' coasts are excluded.
# Stats are calculated on the raw rasters (no coastal gap-filling, which extrapolates).
# Assumes scripts_for_outputs/00_run_all.R has already been run (provides
# `assigned_farms`, `focus_countries`, `GLW_ls`, and paths from src/dirs.R).

library(sf)
library(terra)
library(rnaturalearth)
library(here)

here("src", "country_bboxes.R") %>% source()    # country_info()

HS_dir    <- file.path(outs_path, "data", "habitat_suitability_decadal")
lyr_names <- c("ref_2015_2024", names(GLW_ls))
coast_buffer_m <- 62500

# Country polygons used to define the coastal zone (same source as the documentation qmds)
world <- ne_countries(scale = "large", returnclass = "sf")

# Mean across species rasters, weighted by each model's current production
weighted_HS <- function(rasts, weights) {
  w <- weights / sum(weights)
  Reduce(`+`, map2(rasts, w, \(r, wi) r * wi))
}

# Mean HS per period for each species raster plus the production-weighted stack,
# and the % change from the reference period. One row per model, wide by period.
summarise_HS <- function(rasts, weighted) {
  map2_dfr(
    c(rasts, list(production_weighted = weighted)),
    c(names(rasts), "production_weighted"),
    function(r, nm) {
      global(r, "mean", na.rm = TRUE) %>%
        rownames_to_column("period") %>%
        mutate(model_name = nm)
    }
  ) %>%
    rename(mean_HS = mean) %>%
    group_by(model_name) %>%
    mutate(
      change_pct = 100 * (mean_HS - mean_HS[period == "ref_2015_2024"]) / mean_HS[period == "ref_2015_2024"]
    ) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = model_name,
      names_from = period,
      values_from = c(mean_HS, change_pct),
      names_glue = "{.value}_{period}"
    ) %>%
    dplyr::select(-change_pct_ref_2015_2024)
}

# this_iso <- focus_countries$ISO3_Code[1]
# out_dir <- focus_countries$dir_name[1]
HS_change_stats <- function(this_iso, out_dir) {
  these_farms <- assigned_farms %>%
    st_drop_geometry() %>%
    filter(ISO3_Code == this_iso)
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

  # One SpatRaster per species: layers = ref + each GLW, scaled to 0-1
  HS_rasts <- map(
    set_names(HS_models),
    function(m) {
      r <- rast(file.path(HS_path, m, paste0("MR_", lyr_names, ".tif"))) / 365
      names(r) <- lyr_names
      r
    }
  )

  HS_weights <- these_farms %>%
    distinct(model_name, value) %>%
    filter(model_name %in% HS_models)
  HS_weights  <- set_names(HS_weights$value, HS_weights$model_name)[HS_models]
  HS_weighted <- weighted_HS(HS_rasts, HS_weights)

  # 1. Whole EEZ
  eez_stats <- summarise_HS(HS_rasts, HS_weighted) %>%
    mutate(zone = "EEZ")

  # 2. Within ~50 km of this country's coast
  this_land <- world %>% filter(iso_a3 == this_iso | adm0_a3 == this_iso)
  coastal_stats <- if (nrow(this_land) == 0) {
    warning(this_iso, ": no country polygon found, coastal stats not calculated")
    tibble()
  } else {
    coast_zone <- this_land %>%
      vect() %>%
      project(crs(HS_weighted)) %>%
      buffer(coast_buffer_m)

    summarise_HS(
      map(HS_rasts, \(r) mask(r, coast_zone)),
      mask(HS_weighted, coast_zone)
    ) %>%
      mutate(zone = "coastal_50km")
  }

  out <- bind_rows(eez_stats, coastal_stats) %>%
    mutate(
      ISO3_Code = this_iso,
      country = country_info(this_iso)$name,
      production_weight = unname(HS_weights[model_name] / sum(HS_weights)),
      across(c(starts_with("mean_HS"), starts_with("change_pct")), \(x) round(x, 3))
    ) %>%
    dplyr::select(
      zone, ISO3_Code, country, model_name, production_weight,
      starts_with("mean_HS"), starts_with("change_pct")
    )

  write_csv(out, file.path(out_dir, paste0("HS_change_stats_", this_iso, ".csv")), na = "")
}

walk2(focus_countries$ISO3_Code, focus_countries$dir_name, HS_change_stats)
