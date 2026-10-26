# production_troubleshoot.R
#
# Investigate a single farm flagged by `check_production_files.R`.
#
# Set `farm_ID_to_check` below and run the whole script. It will:
#   1. Print the farm's country and assigned model species
#   2. Plot the daily temperature time series used to drive the model
#   3. Read the saved production file and save one `plot_sim_results_troubleshooting()` plot per selected harvest year (see `years_to_plot`; default is 2025, 2099, and one year per decade between) into outputs/production_troubleshooting/. This folder is overwritten each run, since only one farm is investigated at a time.

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
farm_ID_to_check <- 88018

# Harvest years to plot in step 6. Default: 2025, 2099, and one per decade between.
years_to_plot <- c(2025, seq(2030, 2090, by = 10), 2099)

# production_file <- file.path(outs_path, "data", "farm_production") %>% 
#   list.files(recursive = T, full.names = T, pattern = "MC5000.qs") %>% 
#   str_subset(fix_int(farm_ID_to_check, 8))

production_file <- file.path(
  outs_path, "data", "farm_production",  
  paste0("production_results_farmID_", fix_int(farm_ID_to_check, 8), "_MC5000.qs")
)

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

# ---------- 1b. Other farms in the same country with the same model species ----------
# The issues CSV only records the production file basename, so recover farm_ID
# from the zero-padded number in the filename. One row per farm-year, hence distinct().
issue_farm_IDs <- file.path(outs_path, "production_check_issues.csv") |>
  read_csv(show_col_types = FALSE) |>
  mutate(farm_ID = str_extract(file, "(?<=farmID_)\\d+") |> as.integer()) |>
  distinct(farm_ID) |>
  pull(farm_ID)

other_farms <- assigned_farms |>
  filter(
    country == farm_info$country,
    model_name == farm_info$model_name,
    farm_ID != farm_ID_to_check
  ) |>
  distinct(farm_ID) |>
  mutate(has_issue = farm_ID %in% issue_farm_IDs)

cat(sprintf(
  "Other %s farms in %s: %d\n  ...of which flagged in production_check_issues.csv: %d\nThis farm flagged: %s\n\n",
  farm_info$model_name, farm_info$country,
  nrow(other_farms), sum(other_farms$has_issue),
  if (farm_ID_to_check %in% issue_farm_IDs) "yes" else "no"
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

years_available <- intersect(as.character(years_to_plot), names(prod_data))
if (length(years_available) < length(years_to_plot)) {
  warning(
    "Requested years not found in production file: ",
    paste(setdiff(as.character(years_to_plot), names(prod_data)), collapse = ", ")
  )
}

walk(
  years_available,
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
