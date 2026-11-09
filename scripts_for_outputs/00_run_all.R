library(tidyverse)
library(qs2)
library(here)
library(janitor)

source(here("src/dirs.R"))

# General data shared across all scripts
reports_dir    <- file.path(outs_path, "for_report")
assigned_farms <- file.path(prepdata_path, "assigned_farms_final.qs") %>% 
  qs_read() %>% 
  arrange(ISO3_Code)
# assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") %>% qs_read() # TEMPORARY
# land           <- file.path(prepdata_path, "my_world.qs") %>% qs_read() # may not need this

GLW_ls <- list(
  "GLW_1.5" = 2018:2037,
  "GLW_2.0" = 2032:2051,
  "GLW_3.0" = 2055:2074,
  "GLW_4.0" = 2075:2094
)
GLW_labels <- c("GLW_1.5" = "1.5°C", "GLW_2.0" = "2°C", "GLW_3.0" = "3°C", "GLW_4.0" = "4°C")
year_breaks <- seq(2025, 2095, 10)

report_theme <- function() {
  theme_classic(base_size = 11) +
    theme(
      # plot.title            = element_text(face = "bold", size = 14),
      # plot.subtitle         = element_text(size = 10, margin = margin(b = 6)),
      # plot.title.position   = "plot",
      # plot.caption.position = "plot",
      # plot.caption          = element_text(hjust = 0, colour = "grey45", size = 9),
      legend.position       = "top",
      legend.justification  = "center",
      legend.text           = element_text(size = 10),
      legend.title          = element_text(face = "bold"),
      legend.margin         = margin(b = -4),
      panel.grid.major.y    = element_blank(),
      panel.grid.minor      = element_blank(),
      panel.grid.major.x    = element_line(colour = "grey92", linewidth = 0.4),
      axis.text.y           = element_text(colour = "grey25", size = 10),
      axis.title.x          = element_text(margin = margin(t = 6)),
      plot.margin           = margin(10, 15, 10, 10)
    )
  }

# Height (inches) for the dot-and-line GLW plots, scaled to the number of models on the y axis.
# `base` covers the title, legend and x axis; `min_height` stops 1-2 model plots looking squashed.
GLW_plot_height <- function(n_models, per_model = 0.65, base = 1.6, min_height = 3) {
  max(min_height, base + per_model * n_models)
}

# Step 1: Choose which countries you want to run each script for
focus_countries <- file.path(prepdata_path, "assigned_farms_final.qs") %>% 
  qs_read() %>% 
  distinct(country, ISO3_Code) %>% 
  mutate(
    clean_country = make_clean_names(country, case = "title"),
    dir_name = file.path(reports_dir, paste0(clean_country, " - ", ISO3_Code))
  ) # %>% mutate(dir_name = file.path(reports_dir))        # temporary so I can see them all together
walk(focus_countries$dir_name, dir.create, showWarnings = F)

# Step 2: Create maps of the farm locations
source(here("scripts_for_outputs/01_map_of_farms.R"))

# Step 3: Collate country summary tables (farms per model, excluded farms, FAO production) into one csv per country
source(here("scripts_for_outputs/02_country_info.R"))

  # Step 3a: consolidate all the stats together so it's easier to transfer
  file.path(outs_path, "for_report") %>% 
    list.files(recursive = T, pattern = "country_info_", full.name = T) %>% 
    str_subset("all.csv", negate = T) %>% 
    map_dfr(read_csv, show_col_types = FALSE) %>% 
    # Exclude these for now, once decisions.xlsx is up to date these will be correct
    filter(table != "excluded_farms") %>% 
    arrange(table, -production_t) %>% 
    dplyr::select(-any_of(c("n_farms", "production_perfarm_t", "reason"))) %>% 
    write.csv(file.path(outs_path, "for_report", "country_info_all.csv"), na = "", row.names = FALSE)

# Step 4: Dot-and-line plot of projected biomass change per model at each GLW
source(here("scripts_for_outputs/03_production_change.R"))

# Step 5: Dot-and-line plot of projected change in feed per unit biomass per model at each GLW
source(here("scripts_for_outputs/04_feed_per_biomass_change.R"))

# Step 6: Maps of change in habitat suitability at each GLW, per species and production-weighted
source(here("scripts_for_outputs/05_habitat_suitability_change.R"))

