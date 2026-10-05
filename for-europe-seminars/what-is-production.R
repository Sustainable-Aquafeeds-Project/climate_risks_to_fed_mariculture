library(tidyverse)
library(here)
library(qs2)
library(units)

source("src/dirs.R")
source("src/model_functions.R")
source("src/other_functions.R")

prod_files <- file.path(outs_path, "data", "farm_production") %>% 
  list.files(recursive = T, full.names = T, pattern = "MC5000.qs")
farms_done <- as.integer(str_extract(basename(prod_files), "(?<=farmID_)\\d+"))

assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs_read() %>% 
  sf::st_drop_geometry() %>% 
  mutate(prod_done = farm_ID %in% farms_done)

# The purpose of this script is to understand what production MEANS in terms of the model names

FAO_with_models <- file.path(prepdata_path, "FAO_with_models.qs") %>% 
  qs_read() %>% 
  filter(period %in% 2014:2023) %>% 
  group_by(ISO3_Code, country, model_name, Scientific_Name) %>% 
  reframe(value = sum(value))

country_summary <- function(data, iso, digits = 1) {
  d <- data |>
    filter(ISO3_Code == iso) |>
    mutate(value = drop_units(value))  # drop [t] units

  if (nrow(d) == 0) stop("No rows found for ISO3_Code '", iso, "'")

  total <- sum(d$value, na.rm = TRUE)

  # % of country total per model
  models <- d |>
    summarise(value = sum(value, na.rm = TRUE), .by = model_name) |>
    mutate(pct = 100 * value / total) |>
    arrange(desc(value))

  # % within each model per species
  species <- d |>
    summarise(value = sum(value, na.rm = TRUE), .by = c(model_name, Scientific_Name)) |>
    mutate(pct = 100 * value / sum(value), .by = model_name) |>
    arrange(model_name, desc(value))

  fmt_t <- \(x) format(round(x), big.mark = ",", trim = TRUE)
  fmt_p <- \(x) sprintf("%.*f%%", digits, x)

  cat(sprintf("%s (%s)\n", d$country[1], iso))
  cat(sprintf("Total production: %s t\n", fmt_t(total)))
  cat(sprintf("Number of models: %d\n", nrow(models)))

  for (i in seq_len(nrow(models))) {
    m <- models$model_name[i]
    cat(sprintf("\n  %s: %s t (%s of country total)\n",
                m, fmt_t(models$value[i]), fmt_p(models$pct[i])))

    sp <- filter(species, model_name == m)
    cat(sprintf("    - %s: %s t (%s of model)\n",
                sp$Scientific_Name, fmt_t(sp$value), fmt_p(sp$pct)), sep = "")
  }

  invisible(list(total = total, models = models, species = species))
}

country_summary(FAO_with_models, "JPN")

