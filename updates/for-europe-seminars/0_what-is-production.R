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

disallowed_models <- c("tuna", "nile_tilapia", "olive_flounder", "drum_croaker", "mud_crab", "golden_pompano")
assigned_farms %>% 
  group_by(ISO3_Code, country) %>% 
  mutate(models = n_distinct(model_name)) %>% 
  ungroup() %>% 
  group_by(ISO3_Code, country, models) %>% 
  summarise(disallowed = any(model_name %in% disallowed_models, na.rm = TRUE), .groups = "drop") |>
  filter(!disallowed, models > 1) |>
  select(ISO3_Code, country, models) %>% 
  arrange(ISO3_Code)

# The purpose of this script is to understand what production MEANS in terms of the model names

