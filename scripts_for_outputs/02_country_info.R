# Collate country-level summary tables (farms per model, excluded farms, FAO
# production coverage) into one csv per focus country, for copying into report tables.
# Assumes scripts_for_outputs/00_run_all.R has already been run (provides
# `assigned_farms`, `focus_countries`, and paths from src/dirs.R).

library(sf)
library(readxl)
library(here)

here("src", "other_functions.R") %>% source()    # country_summary()
here("src", "country_bboxes.R")  %>% source()    # country_info()

# Farms excluded after the production run error check. These are no longer in
# `assigned_farms` (the final set), so map farm_ID -> ISO via the pre-exclusion file.
farm_iso_lookup <- file.path(prepdata_path, "assigned_farms.qs") %>%
  qs_read() %>%
  st_drop_geometry() %>%
  distinct(farm_ID, ISO3_Code)

excluded_farms <- file.path(outs_path, "production_check_high_error_files_decisions.xlsx") %>%
  read_excel() %>%
  filter(Conclusion == "Exclude") %>%
  mutate(farm_ID = as.numeric(farmID)) %>%
  left_join(farm_iso_lookup, by = join_by(farm_ID))

country_info_table <- function(this_iso, out_dir) {
  these_farms <- assigned_farms %>%
    st_drop_geometry() %>%
    filter(ISO3_Code == this_iso)
  models_here <- sort(unique(these_farms$model_name))

  # FAO reports "tuna", but farms are assigned "tuna_full" (JPN) or "tuna_ranch" (elsewhere)
  is_included <- function(m) {
    m %in% models_here | (m == "tuna" & any(str_starts(models_here, "tuna_")))
  }

  # 1. Assigned farms per model and the production each model represents
  farms_per_model <- these_farms %>%
    group_by(model_name) %>%
    reframe(
      n_farms = n(),
      production_t = unique(value),
      production_perfarm_t = unique(prod_perfarm)
    ) %>%
    arrange(desc(production_t)) %>%
    mutate(table = "farms_per_model")

  # 2. Farms excluded due to production run errors
  failed_farms <- excluded_farms %>%
    filter(ISO3_Code == this_iso) %>%
    group_by(model_name, reason) %>%
    reframe(n_farms = n()) %>%
    mutate(table = "excluded_farms")

  # 3 & 4. FAO-reported production (2014-2023 mean) by model and by species within model.
  # country_summary() prints to the console; keep only its return value.
  fao <- tryCatch(
    {
      invisible(capture.output(res <- country_summary(this_iso)))
      res
    },
    error = function(e) NULL
  )

  fao_models <- if (is.null(fao)) tibble() else fao$models %>%
    transmute(
      model_name,
      production_t = value,
      pct_of_country_total = pct,
      included = is_included(model_name),
      table = "fao_production_by_model"
    )

  fao_species <- if (is.null(fao)) tibble() else fao$species %>%
    transmute(
      model_name,
      species = Scientific_Name,
      production_t = value,
      pct_of_model = pct,
      included = is_included(model_name),
      table = "fao_production_by_species"
    )

  out <- bind_rows(farms_per_model, failed_farms, fao_models, fao_species) %>%
    mutate(
      ISO3_Code = this_iso,
      country = country_info(this_iso)$name,
      across(where(is.numeric), \(x) round(x, 1))
    ) %>%
    dplyr::select(
      table, ISO3_Code, country, model_name, any_of(c(
        "species", "n_farms", "production_t", "production_perfarm_t",
        "pct_of_country_total", "pct_of_model", "included", "reason"
      ))
    )

  write_csv(out, file.path(out_dir, paste0("country_info_", this_iso, ".csv")), na = "")
}

walk2(focus_countries$ISO3_Code, focus_countries$dir_name, country_info_table)
