library(tidyverse)
library(here)
library(qs2)
library(ggh4x)
library(tidyverse)

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

GLW_ls <- list(
  "GLW_1.5" = 2018:2037,
  "GLW_2.0" = 2032:2051,
  "GLW_3.0" = 2055:2074,
  "GLW_4.0" = 2075:2094
)
all_years <- c(2018:2037, 2032:2051, 2055:2074, 2075:2094) %>% unique()

# Get country data function -------------------------------------------------------------------------------
# fnm <- fnms[1]
# id <- ids[1]
get_country_data <- function(fnms, ids) {
  annual_biomass <- map2_dfr(
    fnms, ids,
    function(fnm, id) {
      qs_read(fnm) %>% 
        tidy_stat("weight_scaled") %>% 
        mutate(farm_ID = id) %>% 
        slice_max(prod_day, n = 1) %>% 
        select(farm_ID, year, mean) %>% 
        filter(year %in% all_years)
    }
  ) %>% rename(biomass = mean)

  annual_food_prov <- map2_dfr(
    fnms, ids,
    function(fnm, id) {
      qs_read(fnm) %>% 
        tidy_stat("food_prov_scaled") %>% 
        mutate(farm_ID = id) %>% 
        group_by(farm_ID, year) %>% 
        reframe(food_prov = sumna(mean)) %>% 
        filter(year %in% all_years)
    }
  ) %>% 
    left_join(annual_biomass, by = join_by(farm_ID, year)) %>% 
    mutate(food_prov_biom = food_prov/biomass)

  left_join(
      annual_food_prov,
      annual_food_prov %>% filter(year == 2025) %>% select(-c(year, biomass)) %>% rename(food_prov_2025 = food_prov, food_prov_biom_2025 = food_prov_biom),
      by = join_by(farm_ID)
    ) %>% 
    left_join(
      annual_biomass %>% filter(year == 2025) %>% select(farm_ID, biomass) %>% rename(biomass_2025 = biomass),
      by = join_by(farm_ID)
    ) %>% 
    mutate(
      food_prov_change = food_prov - food_prov_2025,
      food_prov_biom_change = food_prov_biom - food_prov_biom_2025,
      biomass_change = biomass - biomass_2025
    )
}

# Get country data, save
# iso <- "NOR"
walk(
  unique(assigned_farms$ISO3_Code), 
  function(iso) {

    data_filename <- here("for-europe-seminars", "data", paste0(iso, ".qs"))
    af <- assigned_farms %>% filter(ISO3_Code == iso)

    if (
      # !file.exists(data_filename) & 
      any(af$prod_done)
    ) {

      df <- af %>%
        group_by(ISO3_Code, model_name) %>% 
        reframe(total_value = sum(value), value_done = sum(value[prod_done %in% TRUE]))

      fnms <- prod_files %>% str_subset(paste0(fix_int(unique(af$farm_ID), 8), collapse = "|"))
      ids <- as.integer(str_extract(basename(fnms), "(?<=farmID_)\\d+"))

      country_data <- get_country_data(fnms, ids) %>% 
        right_join(af, by = join_by(farm_ID))

      failed_farms <- nrow(country_data[is.na(country_data$biomass), ])
      country_data <- country_data %>% filter(!is.na(biomass))

      # Special for Japan - there's one farm that didn't seem to work but isn't being filtered out
      if (iso == "JPN") {
        weird_farms <- country_data$farm_ID[country_data$food_prov_biom_change > 5]
        country_data <- country_data %>% filter(!farm_ID %in% weird_farms)
        failed_farms <- failed_farms + length(weird_farms)
      }

      qs_save(country_data, data_filename)
    }
  },
  .progress = T
)

# Plot in terms of GLW change --------------------------------------------------------------------------------
fnms <- here("for-europe-seminars", "data") %>% 
  list.files(full.names = T) %>% 
  str_subset(paste0(
    c("ESP", "MMR"), # Chosen ocountries
    collapse = "|"
  ))

lims1 <- c(-0.65, 0.40)
breaks1 <- seq(-0.80, 0.40, 0.2)
lims2 <- c(-0.1, 0.2)
breaks2 <- seq(-1, 3, 0.05)
warming_cols <- c("1.5°C" = "#E9A23B", "2°C"   = "#DB5F32", "3°C"   = "#B3243B", "4°C"   = "#5C0F2E")
warming_levs <- c("1.5°C" = "GLW_1.5", "2°C" = "GLW_2.0", "3°C" = "GLW_3.0", "4°C" = "GLW_4.0")

walk(
  fnms,
  function(fnm) {
    df <- qs_read(fnm)
    all_prod <- df %>% 
      filter(year == 2025) %>% 
      group_by(model_name) %>% 
      reframe(biomass = sum(biomass))
    all_prod$perc_prod <- all_prod$biomass/sum(all_prod$biomass)
    all_prod <- all_prod %>% 
      mutate(label = paste0(model_name, " (", round(100*perc_prod, 2), "%)"))

    # yrs <- GLW_ls[[1]]
    # glw <- names(GLW_ls)[[1]]
    plot_df <- map2_dfr(
      GLW_ls,
      names(GLW_ls),
      function(yrs, glw) {
        df %>% 
          filter(year %in% yrs) %>% 
          mutate(
            biomass_change = biomass_change/biomass_2025,
            food_prov_biom_change = food_prov_biom_change/food_prov_biom_2025
          ) %>% 
          group_by(ISO3_Code, country, model_name) %>% 
          reframe(
            biomass_2025 = first(biomass_2025),
            biomass_change_mean = meanna(biomass_change),
            biomass_change_sd   = sdna(biomass_change),
            biomass_change_min  = minna(biomass_change),
            biomass_change_max  = maxna(biomass_change),
            food_prov_biom_change_mean = meanna(food_prov_biom_change),
            food_prov_biom_change_sd   = sdna(food_prov_biom_change),
            food_prov_biom_change_min  = minna(food_prov_biom_change),
            food_prov_biom_change_max  = maxna(food_prov_biom_change)
          ) %>% 
          mutate(GLW = factor(glw, levels = warming_levs, labels = names(warming_levs)))
      }
    )

    p1 <- plot_df %>% 
      left_join(all_prod, by = "model_name") %>% 
      ggplot(aes(
        x = label, y = biomass_change_mean, 
        ymin = biomass_change_mean - biomass_change_sd,
        ymax = biomass_change_mean + biomass_change_sd,
        colour = GLW
      )) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
      geom_errorbar(position = position_dodge(width = 0.5), width = 0.25) +
      scale_y_continuous(
        limits = lims1, breaks = breaks1,
        labels = scales::label_percent()
      ) +
      scale_colour_manual(values = warming_cols) +
      prettyplot() +
      theme(legend.position = "top") +
      labs(
        x = NULL, y = "Change in biomass",
        title = paste0("Projected change in farm production (", str_remove(basename(fnm), ".qs"), ")"),
        subtitle = "Mean and SD between farms"
      )
    
    ggsave(
      plot = p1,
      filename = here("for-europe-seminars", "plots", 
        paste0("biomass_change_", str_replace(basename(fnm), ".qs", ".png"))
      ),
      width = 6, height = 5
    )

    # plot_df %>% 
    #   group_by(year) %>%
    #   mutate(w = biomass_2025 / sum(biomass_2025)) %>%
    #   summarise(
    #     biomass_change_sd   = sqrt(sum(w^2 * biomass_change_sd^2)),
    #     biomass_change_mean = sum(w * biomass_change_mean),
    #     food_prov_biom_change_sd   = sqrt(sum(w^2 * food_prov_biom_change_sd^2)),
    #     food_prov_biom_change_mean = sum(w * food_prov_biom_change_mean),
    #     n_models    = n()
    #   )
    
    p2 <- plot_df %>% 
      left_join(all_prod, by = "model_name") %>% 
      ggplot(aes(
        x = label, y = food_prov_biom_change_mean, 
        ymin = food_prov_biom_change_mean - food_prov_biom_change_sd,
        ymax = food_prov_biom_change_mean + food_prov_biom_change_sd,
        colour = GLW
      )) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
      geom_errorbar(position = position_dodge(width = 0.5), width = 0.25) +
      scale_y_continuous(
        limits = lims2, breaks = breaks2,
        labels = scales::label_percent()
      ) +
      scale_colour_manual(values = warming_cols) +
      prettyplot() +
      theme(legend.position = "top") +
      labs(
        x = NULL, y = "Change in feed required (%)",
        title = paste0("Projected change in feed required (", str_remove(basename(fnm), ".qs"), ")"),
        subtitle = "Mean and SD between farms"
      )
    ggsave(
      plot = p2,
      filename = here("for-europe-seminars", "plots", 
        paste0("feedrequired_change_", str_replace(basename(fnm), ".qs", ".png"))
      ),
      width = 6, height = 5
    )
  },
  .progress = T
)
