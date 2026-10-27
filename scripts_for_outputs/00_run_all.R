library(tidyverse)
library(qs2)
library(here)
library(janitor)

source(here("src/dirs.R"))

# General data shared across all scripts
reports_dir    <- file.path(outs_path, "for_report")
assigned_farms <- file.path(prepdata_path, "assigned_farms_final.qs") %>% qs_read()
land           <- file.path(prepdata_path, "my_world.qs") %>% qs_read()

GLW_ls <- list(
  "GLW_1.5" = 2018:2037,
  "GLW_2.0" = 2032:2051,
  "GLW_3.0" = 2055:2074,
  "GLW_4.0" = 2075:2094
)
GLW_labels <- c("GLW_1.5" = "1.5°C", "GLW_2.0" = "2°C", "GLW_3.0" = "3°C", "GLW_4.0" = "4°C")
year_breaks <- seq(2025, 2095, 10)

# Step 1: Choose which countries you want to run each script for
focus_countries <- file.path(prepdata_path, "assigned_farms_final.qs") %>% 
  qs_read() %>% 
  distinct(country, ISO3_Code) %>% 
  mutate(
    clean_country = make_clean_names(country, case = "title"),
    dir_name = file.path(reports_dir, paste0(ISO3_Code, " - ", clean_country))
  ) %>% 
  mutate(dir_name = file.path(reports_dir))        # temporary so I can see them all together
walk(focus_countries$dir_name, dir.create, showWarnings = F)

# Create maps of the farm locations
source(here("scripts_for_outputs/01_map_of_farms.R"))

