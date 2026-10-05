library(tidyverse)
library(here)
library(qs2)

source("src/dirs.R")
source("src/model_functions.R")
source("src/other_functions.R")
source("src/species_colours.R")

prod_files <- file.path(outs_path, "data", "farm_production") %>% 
  list.files(recursive = T, full.names = T, pattern = "MC5000.qs")
farms_done <- as.integer(str_extract(basename(prod_files), "(?<=farmID_)\\d+"))

assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs_read() %>% 
  sf::st_drop_geometry() %>% 
  mutate(prod_done = farm_ID %in% farms_done)

# Get country data function
get_country_data <- function(fnms, ids) {
  annual_biomass <- map2_dfr(
    fnms, ids,
    function(fnm, id) {
      qs_read(fnm) %>% 
        tidy_stat("weight_scaled") %>% 
        mutate(farm_ID = id) %>% 
        slice_max(prod_day, n = 1) %>% 
        select(farm_ID, year, mean)
    }
  ) %>% rename(biomass = mean)

  annual_food_prov <- map2_dfr(
    fnms, ids,
    function(fnm, id) {
      qs_read(fnm) %>% 
        tidy_stat("food_prov_scaled") %>% 
        mutate(farm_ID = id) %>% 
        group_by(farm_ID, year) %>% 
        reframe(food_prov = sumna(mean))
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
  
# iso <- "BLZ"
walk(
  sort(unique(assigned_farms$ISO3_Code)), 
  function(iso) {

    plot_filename <- here("for-europe-seminars", "plots", paste0(iso, ".png"))
    af <- assigned_farms %>% filter(ISO3_Code == iso)

    if (
      !file.exists(plot_filename) & 
      any(af$prod_done) &
      sum(af$prod_done * af$value) > sum(af$value) * 0.5
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
     
      p <- country_data %>% 
        select(year, model_name, food_prov_change, food_prov_biom_change, biomass_change) %>% 
        mutate(
          food_prov_change = food_prov_change/1000000, 
          food_prov_biom_change = food_prov_biom_change * 100, 
          biomass_change = biomass_change/1000000
        ) %>% 
        pivot_longer(cols = c("food_prov_change", "food_prov_biom_change", "biomass_change"), names_to = "measure", values_to = "value") %>% 
        mutate(
          measure = factor(
            measure, 
            levels = c("biomass_change", "food_prov_change", "food_prov_biom_change"),
            labels = c("Change in biomass (t)", "Change in feed provided (t)", "Change in feed provided per biomass (%)")
          )
        ) %>% 
        group_by(year, model_name, measure) %>% 
        reframe(
          mean_value = mean(value),
          sd_value = sd(value)
        ) %>% 
        ggplot(aes(
          x = year, y = mean_value, 
          # ymin = mean_value-sd_value, ymax = mean_value+sd_value, 
          colour = model_name, fill = model_name
        )) +
        # geom_point(shape = 19, size = 1, alpha = 0.15) +
        # geom_smooth(alpha = 0.25) +
        geom_line(linewidth = 1) +
        # geom_ribbon(alpha = 0.1, linewidth = 0) +
        scale_colour_manual(values = models_pal) +
        scale_fill_manual(values = models_pal) +
        prettyplot() +
        theme(legend.position = "bottom") +
        facet_wrap(~measure, nrow = 2, scales = "free") +
        labs(
          title = paste(unique(af$country)),
          subtitle = paste0(round(100*sum(df$value_done)/sum(df$total_value), 1), "% of production done, ", failed_farms, " farms failed")
        )

      rm(country_data)

      ggsave(
        plot = p,
        filename = plot_filename,
        height = 7.5, width = 8.5
      )
    }
  },
  .progress = T
)

assigned_farms <- assigned_farms %>% 
  filter(ISO3_Code %in% c("USA", "PER", "NOR", "CHN", "MDG", "IDN")) %>% 
  mutate(country = factor(
    country, 
    levels = c("United States of America", "Norway", "China", "Peru", "Madagascar", "Indonesia")
  ))

p_data <- map_dfr(
  unique(assigned_farms$ISO3_Code),
  function(iso) {
  af <- assigned_farms %>% filter(ISO3_Code == iso)
  fnms <- prod_files %>% str_subset(paste0(fix_int(unique(af$farm_ID), 8), collapse = "|"))
  ids <- as.integer(str_extract(basename(fnms), "(?<=farmID_)\\d+"))
  get_country_data(fnms, ids) %>% 
    right_join(af, by = join_by(farm_ID))
  },
  .progress = T
)

# 6-country plots
p_biomass <- p_data %>% 
  select(country, year, model_name, biomass_change) %>% 
  mutate(biomass_change = biomass_change/1000000) %>% 
  group_by(country, year, model_name) %>% 
  reframe(mean_value = mean(biomass_change)) %>% 
  ggplot(aes(x = year, y = mean_value, colour = model_name, fill = model_name)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = models_pal) +
  scale_fill_manual(values = models_pal) +
  labs(x = "Year", y = "Change in biomass (t)") +
  prettyplot() +
  theme(legend.position = "right") +
  facet_wrap(~country, nrow = 2, scales = "free")

ggsave(plot = p_biomass, filename = "p_biomass.png", width = 12)

p_feed <- p_data %>% 
  select(country, year, model_name, food_prov_biom_change) %>% 
  mutate(food_prov_biom_change = food_prov_biom_change * 100) %>% 
  group_by(country, year, model_name) %>% 
  reframe(mean_value = mean(food_prov_biom_change)) %>% 
  ggplot(aes(x = year, y = mean_value, colour = model_name, fill = model_name)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = models_pal) +
  scale_fill_manual(values = models_pal) +
  labs(x = "Year", y = "Change in feed required (%)") +
  prettyplot() +
  theme(legend.position = "right") +
  facet_wrap(~country, nrow = 2, scales = "free")

ggsave(plot = p_feed, filename = "p_feed.png", width = 12)
