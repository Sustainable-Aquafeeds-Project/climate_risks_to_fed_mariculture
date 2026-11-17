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

country_info_table <- function(this_iso = "JPN", out_dir) {
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
      # FAO_production_t = unique(value),
      production_perfarm_t = unique(prod_perfarm)
    ) %>%
    arrange(desc(n_farms))

  # 2. Farms excluded due to production run errors
  failed_farms <- excluded_farms %>%
    filter(ISO3_Code == this_iso) %>%
    group_by(model_name) %>%
    summarise(
      farms_excl = n(),
      reasons = paste(unique(reason), collapse = "; "),
      .groups = "drop"
    )
  
  # Join to the farms per model table
  farms_per_model <- farms_per_model %>% 
    left_join(failed_farms, by = join_by(model_name)) %>% 
    mutate(model_name = case_when(model_name %in% c("tuna_full", "tuna_ranch") ~ "tuna", T ~ model_name))

  # 3 & 4. FAO-reported production (2015-2023 mean): country total with the share covered by included models, and by species within model. country_summary() prints to the console; keep only its return value.
  fao <- tryCatch(
    {
      invisible(capture.output(res <- country_summary(this_iso)))
      res
    },
    error = function(e) NULL
  )

  fao_country <- if (is.null(fao)) tibble() else fao$models %>%
    reframe(
      production_t = sum(value),
      pct_included = 100 * sum(value[is_included(model_name)]) / sum(value)
    ) %>%
    mutate(model_name = "ALL", included = TRUE)

  fao_species <- if (is.null(fao)) tibble() else fao$species %>%
    transmute(
      model_name,
      species = Scientific_Name,
      production_t = value,
      pct_of_model = pct,
      included = is_included(model_name)
    )

  out <- bind_rows(fao_country, fao_species) %>%
    left_join(farms_per_model, by = join_by(model_name), all) %>% 
    mutate(
      ISO3_Code = this_iso,
      country = country_info(this_iso)$name,
      across(where(is.numeric), \(x) round(x, 1))
    )
  
  m_order <- out %>% arrange(-production_t) %>% pull(model_name)

  out %>% 
    arrange(
      -included, 
      factor(model_name, levels = m_order), 
      pct_of_model
    )

  write_csv(out, file.path(out_dir, paste0("country_info_", this_iso, ".csv")), na = "")
}

walk2(focus_countries$ISO3_Code, focus_countries$dir_name, country_info_table, .progress = T)
