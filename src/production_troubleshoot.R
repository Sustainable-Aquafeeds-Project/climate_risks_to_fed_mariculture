# production_troubleshoot.R
#
# Investigate a single farm flagged by `check_production_files.R`.
#
# Set `farm_ID_to_check` below and run the whole script. It will:
#   1. Print the farm's country and assigned model species
#   2. Plot the daily temperature time series used to drive the model
#   3. Re-run the growth model for every harvest year (2025-2099), producing
#      the FULL set of diagnostic output variables (not just the reduced set
#      saved in the production files), and save one `plot_sim_results()` plot
#      per year into outputs/production_troubleshooting/. This folder is
#      overwritten each run, since only one farm is investigated at a time.
#
# Re-running (rather than just reading the saved production file) is
# necessary because the saved files only keep a reduced set of ~15 response
# variables - useful diagnostics like SGR, O2, NH4, and the elemental budgets
# aren't stored, but they are exactly what `plot_sim_results()` is designed
# to show.

library(here)
library(tidyverse)
library(qs2)
library(readxl)
library(future)
library(furrr)
library(slider)

source(here("src", "dirs.R"))
source(here("src", "model_functions.R"))
source(here("src", "other_functions.R"))

# =========================================================================
# ---- SET THE FARM TO INVESTIGATE HERE ----
farm_ID_to_check <- 00008951

# production_file <- file.path(outs_path, "data", "farm_production") %>% 
#   list.files(recursive = T, full.names = T, pattern = "MC5000.qs") %>% 
#   str_subset(fix_int(farm_ID_to_check, 8))

production_file <- file.path(outs_path, "data", "farm_production", "new_local", "production_results_farmID_00008951_MC100.qs")


# =========================================================================
troubleshoot_path <- file.path(outs_path, "production_troubleshooting")
dir.create(troubleshoot_path, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. Farm info ----------
assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") |>
  qs_read() |>
  sf::st_drop_geometry()

farm_info <- assigned_farms |> filter(farm_ID == farm_ID_to_check)

if (nrow(farm_info) == 0) {
  stop("farm_ID ", farm_ID_to_check, " was not found in assigned_farms.qs")
}
if (nrow(farm_info) > 1) {
  warning(
    "Multiple rows found for farm_ID ", farm_ID_to_check,
    " in assigned_farms.qs; using the first."
  )
  farm_info <- farm_info[1, ]
}

# Attach the feed assignment (mirrors 05_01_run_production.qmd)
feed_assignments <- file.path(rawdata_path, "feed_profiles_and_composition.xlsx") |>
  read_excel(sheet = "feed_assignments") |>
  select(model_name, country, feed_class)

farm_info <- farm_info |>
  left_join(feed_assignments, by = join_by(model_name, country)) |>
  mutate(feed_name = paste(model_name, feed_class, sep = "-"))

if (is.na(farm_info$feed_class)) {
  stop(
    "No feed_class assignment found for model_name = '", farm_info$model_name,
    "', country = '", farm_info$country, "'. Check feed_assignments in ",
    "feed_profiles_and_composition.xlsx."
  )
}

cat(sprintf(
  "farm_ID:       %d\nCountry:       %s (%s)\nModel species: %s\nFeed:          %s\n\n",
  farm_info$farm_ID, farm_info$country, farm_info$ISO3_Code,
  farm_info$model_name, farm_info$feed_name
))

# ---------- 2. Species and feed parameters ----------
sp <- file.path(prepdata_path, "species_parameters", paste0("species_parameters_", farm_info$model_name, ".qs")) |>
  qd_read()

fp <- file.path(prepdata_path, "feed_parameters", paste0(farm_info$feed_name, ".qs")) |>
  qs_read()

gd <- sp[["gd"]] |> as.integer()

# ---------- 3. Feed components ----------
fp_df <- tibble(
  Macronutrient = c("Protein", "Lipid", "Carbohydrate"),
  Total = c(fp[["protein"]]$macro, fp[["lipid"]]$macro, fp[["carb"]]$macro),
  Digestible = c(fp[["protein"]]$digestibility, fp[["lipid"]]$digestibility, fp[["carb"]]$digestibility)
) %>% 
  mutate(Digestible = Digestible * Total) %>% 
  pivot_longer(cols = c(Total, Digestible), names_to = "component", values_to = "value")

ggsave(
  file.path(troubleshoot_path, "feed_components.png"),
  ggplot(fp_df, aes(x = Macronutrient, y = value, fill = component)) +
    geom_col(alpha = 0.5, colour = "black", position = "identity") +
    prettyplot() +
    theme(legend.position = "right") +
    scale_fill_brewer(palette = "Set1"),
  height = 6.5,
  width = 7.5
)

# ---------- 4. Daily temperatures ----------
temp_file <- file.path(prepdata_path, "assigned_farms_alltemps", paste0("farm_", fix_int(farm_info$farm_ID, 8), ".qs"))

if (!file.exists(temp_file)) {
  stop("Could not find a temperature file for farm_ID ", farm_ID_to_check, " at: ", temp_file)
}

temps_all <- qs_read(temp_file) |>
  mutate(
    doy  = as.integer(format(date, "%j")),
    year = year(date)
  )

ggsave(
  file.path(troubleshoot_path, "daily_SST.png"),
  ggplot(temps_all, aes(x = doy, y = sst, colour = as.factor(year))) +
    geom_line() +
    geom_hline(yintercept = sp["Toa"]) +
    geom_hline(yintercept = sp["Taa"], linetype = "dashed") +
    geom_hline(yintercept = sp["Tma"], linetype = "dashed") +
    theme_classic() +
    labs(
      title = sprintf(
        "Farm %d (%s, %s): daily sea surface temperature",
        farm_info$farm_ID, farm_info$country, farm_info$model_name
      ),
      x = "Date", y = "Temperature (°C)"
    ) +
    scale_y_continuous(breaks = seq(-10, 40, 2.5)) +
    scale_colour_viridis_d() +
    theme(legend.position = "none"),
  height = 7,
  width = 10
)


# ---------- 5. Existing plots ---------------
# Copy some existing plots into the troubleshooting folder
file.copy(
  from = file.path(outs_path, "plots", "farm_production_checks", "farm_biomass", paste0("harvested_biomass_farmID_", farm_ID_to_check, ".png")), 
  to = file.path(troubleshoot_path, paste0("harvested_biomass_farmID_", farm_ID_to_check, ".png")), 
  overwrite = T
)

file.copy(
  from = file.path(outs_path, "plots", "farm_production_checks", "feed_required", paste0("feed_required_farmID_", farm_ID_to_check, ".png")), 
  to = file.path(troubleshoot_path, paste0("feed_required_farmID_", farm_ID_to_check, ".png")), 
  overwrite = T
)

file.copy(
  from = file.path(outs_path, "plots", "farm_production_checks", "feed_required_per_biomass", paste0("feed_required_perbiomass_farmID_", farm_ID_to_check, ".png")), 
  to = file.path(troubleshoot_path, paste0("feed_required_perbiomass_farmID_", farm_ID_to_check, ".png")), 
  overwrite = T
)

# ---------- 6. Examine the saved production output for every year ----------
# The production files only retain a reduced set of ~15 response variables, so we use plot_sim_results_troubleshooting() (see other_functions.R) - a trimmed copy of plot_sim_results() with the panels for variables that aren't saved (SGR, E_assim, E_somat, water_temp, food_enc, ing_pot, O2, NH4) removed.
if (length(production_file) == 0) {
  stop("No MC5000 production file found for farm_ID ", farm_ID_to_check, " under ", file.path(outs_path, "data", "farm_production"))
}
if (length(production_file) > 1) {
  warning("Multiple production files matched farm_ID ", farm_ID_to_check, "; using the first: ", production_file[1])
  production_file <- production_file[1]
}

prod_data <- qs_read(production_file)

walk(
  names(prod_data),
  function(yr) {
    yr_data <- prod_data[[yr]]
    days <- yr_data$days

    plot_df <- imap_dfr(yr_data$stats, function(mat, nm) {
      tibble(days = days, output_var = nm, value = mat[, 1])
    })

    p <- plot_sim_results_troubleshooting(plot_df, CS = sp[["CS"]])

    ggsave(
      file.path(troubleshoot_path, paste0("year_", yr, ".png")),
      plot = p,
      width = 14,
      height = 10,
      dpi = 150
    )
  },
  .progress = TRUE
)

cat(sprintf("\nPlots saved to: %s\n", troubleshoot_path))
