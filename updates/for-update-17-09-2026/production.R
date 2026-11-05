library(tidyverse)
library(here)
library(qs2)
library(sf)

source(here("src/dirs.R"))
source(here("src/other_functions.R"))
source(here("src/model_functions.R"))

production_files <- list.files(file.path(outs_path, "data/farm_production"), full.names = T, pattern = "MC5000.qs")
farms_done <- as.integer(str_extract(basename(production_files), "(?<=farmID_)\\d+"))

assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs_read() %>% 
  mutate(done = case_when(farm_ID %in% farms_done ~ T, T ~ F))

# Taiwan has about 60% of its production done over several species so let's go with that

# assigned_farms %>% 
#   st_drop_geometry() %>% 
#   group_by(country) %>% 
#   mutate(num_models = n_distinct(model_name)) %>% 
#   ungroup() %>% 
#   group_by(country, num_models, done) %>% 
#   reframe(n = n()) %>% 
#   pivot_wider(names_from = done, values_from = n) %>% 
#   mutate(perc_done = `TRUE`/(`TRUE` + `FALSE`)) %>% 
#   arrange(-perc_done)

assigned_farms %>% 
  st_drop_geometry() %>% 
  filter(country == "Taiwan") %>% 
  group_by(model_name, done) %>% 
  reframe(n = n()) %>% 
  arrange(-done, -n)

assigned_farms_TWN <- assigned_farms %>% 
  st_drop_geometry() %>% 
  filter(country == "Taiwan" & done)

# Showing results per country ------------------------------------------------------------

# First, I want to show how all farms have changed in their relative production over time.
farm_biomass <- map_dfr(
  1:nrow(assigned_farms_TWN),
  function(r) {
    prod_df <- production_files %>% 
      str_subset(fix_int(assigned_farms_TWN$farm_ID[r], 8)) %>% 
      qs_read()
    # names(prod_df[[1]][[2]])

    # Tidy farm biomass stat
    weight_scaled <- tidy_stat(prod_df, "weight_scaled") %>% 
      filter(prod_day == max(prod_day)) %>% 
      select(-prod_day, -day, -sd)
    weight_scaled <- weight_scaled %>% 
      mutate(
        mean_2025 = weight_scaled$mean[weight_scaled$year == 2025],
        mean_change = (mean - mean_2025)/mean_2025
      )
    
    # Add farm information
    weight_scaled %>% 
      mutate(
        farm_ID = assigned_farms_TWN$farm_ID[r],
        model_name = assigned_farms_TWN$model_name[r]
      )
  },
  .progress = T
)

# Plot for the whole country
nat_farm_biomass <- farm_biomass %>% 
  group_by(model_name, year) %>% 
  reframe(
    mean_nat_change = meanna(mean_change),
    sd_nat_change = sdna(mean_change),
    n = n()
  ) %>% 
  mutate()

labs <- nat_farm_biomass %>% 
  mutate(lab = paste0(str_to_sentence(str_replace(model_name, "_", " ")), " (", n, ")")) %>% 
  pull(lab) %>% 
  unique() %>% 
  setNames(unique(nat_farm_biomass$model_name))

nat_farm_biomass %>% 
  ggplot(
    aes(
      x = year, 
      y = mean_nat_change, 
      colour = model_name, 
      fill = model_name
    )
  ) +
  geom_point() +
  geom_smooth(alpha = 0.3) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(2025, 2095, 10)) +
  scale_colour_brewer(palette = "Set1", labels = labs) +
  scale_fill_brewer(palette = "Set1", labels = labs) +
  labs(
    title = "Country of Taiwan, ~60% of production represented",
    x = "Year",
    y = "Change in farm harvested biomass"
  ) +
  prettyplot() +
  theme(legend.position = "right")

ggsave(here("for-update-17-09-2026", "nat_biomass.png"), height = 6.5, width = 10)

## Weighted by 2017-2025 production ----------------------------------------------------------------
nat_farm_biomass_rel <- nat_farm_biomass %>% 
  left_join(
    assigned_farms_TWN %>% 
      distinct(model_name, value), # "value" is actually 2025 production (t)
    by = "model_name"
  ) %>% 
  group_by(year) %>% 
  reframe(
    mean = weighted.mean(mean_nat_change, value)
  )

nat_farm_biomass_rel %>% 
  ggplot(
    aes(
      x = year, 
      y = mean
    )
  ) +
  geom_point() +
  geom_smooth(alpha = 0.3, colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(2025, 2095, 10)) +
  labs(
    title = "Country of Taiwan, ~60% of production represented",
    x = "Year",
    y = "Change in farm harvested biomass, weighted by production"
  ) +
  prettyplot() +
  theme(legend.position = "right")

ggsave(here("for-update-17-09-2026", "nat_biomass_weighted.png"), height = 6.5, width = 10)

# Now, show how feed consumption per biomass changes ------------------------------------------
farm_rel_feed <- map_dfr(
  1:nrow(assigned_farms_TWN),
  function(r) {
    prod_df <- production_files %>% 
      str_subset(fix_int(assigned_farms_TWN$farm_ID[r], 8)) %>% 
      qs_read()
    # names(prod_df[[1]][[2]])

    # Tidy stats
    weight_scaled <- tidy_stat(prod_df, "weight_scaled") %>% 
      filter(prod_day == max(prod_day)) %>% 
      select(-prod_day, -day, -sd) %>% 
      rename(biomass = mean)
    food_prov <- tidy_stat(prod_df, "food_prov_scaled") %>% 
      group_by(year) %>% 
      reframe(mean = sumna(mean)) %>% 
      left_join(weight_scaled, by = "year") %>% 
      mutate(mean = mean/biomass) %>% 
      select(-biomass)
    food_prov <- food_prov %>% 
      mutate(
        mean_2025 = food_prov$mean[food_prov$year == 2025],
        mean_change = (mean - mean_2025)/mean_2025
      ) %>% 
      select(-mean)
    
    # Add farm information
    food_prov %>% 
      mutate(
        farm_ID = assigned_farms_TWN$farm_ID[r],
        model_name = assigned_farms_TWN$model_name[r]
      )
  },
  .progress = T
)

# Plot for the whole country
nat_farm_rel_feed <- farm_rel_feed %>% 
  group_by(model_name, year) %>% 
  reframe(
    mean_nat_change = meanna(mean_change),
    sd_nat_change = sdna(mean_change),
    n = n()
  ) %>% 
  mutate()

labs <- nat_farm_rel_feed %>% 
  mutate(lab = paste0(str_to_sentence(str_replace(model_name, "_", " ")), " (", n, ")")) %>% 
  pull(lab) %>% 
  unique() %>% 
  setNames(unique(nat_farm_rel_feed$model_name))

nat_farm_rel_feed %>% 
  ggplot(
    aes(
      x = year, 
      y = mean_nat_change, 
      colour = model_name, 
      fill = model_name
    )
  ) +
  geom_point() +
  geom_smooth(alpha = 0.3) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(2025, 2095, 10)) +
  scale_colour_brewer(palette = "Set1", labels = labs) +
  scale_fill_brewer(palette = "Set1", labels = labs) +
  labs(
    title = "Country of Taiwan, ~60% of production represented",
    x = "Year",
    y = "Change in feed provided per harvested biomass"
  ) +
  prettyplot() +
  theme(legend.position = "right")

ggsave(here("for-update-17-09-2026", "nat_feed.png"), height = 6.5, width = 10)

## Weighted by 2017-2025 production ----------------------------------------------------------------
nat_farm_rel_feed_rel <- nat_farm_rel_feed %>% 
  left_join(
    assigned_farms_TWN %>% 
      distinct(model_name, value), # "value" is actually 2025 production (t)
    by = "model_name"
  ) %>% 
  group_by(year) %>% 
  reframe(
    mean = weighted.mean(mean_nat_change, value)
  )

nat_farm_rel_feed_rel %>% 
  ggplot(
    aes(
      x = year, 
      y = mean
    )
  ) +
  geom_point() +
  geom_smooth(alpha = 0.3, colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(2025, 2095, 10)) +
  labs(
    title = "Country of Taiwan, ~60% of production represented",
    x = "Year",
    y = "Change in feed required per biomass, weighted by production"
  ) +
  prettyplot() +
  theme(legend.position = "right")

ggsave(here("for-update-17-09-2026", "nat_feed_weighted.png"), height = 6.5, width = 10)

# Species responses ---------------------------------------------------------------------------------
species_TWN <- unique(assigned_farms_TWN$model_name)

species_params <- file.path(prepdata_path, "species_parameters") %>% 
  list.files(full.names = T) %>% 
  str_subset(pattern = paste0(species_TWN, collapse = "|"))

species_params_nms <- species_params %>% basename() %>% str_remove_all(".qs|species_parameters_")

species_params <- map(species_params, qd_read) %>% 
  setNames(species_params_nms)

water_temp <- seq(-10, 40, 0.5)

temp_responses <- map(
  species_params, 
  function(sp) {
    FR <- feeding_rate(water_temp = water_temp, species_params = sp)
    FR[is.na(FR)] <- 0
    FR
  }
) %>% setNames(species_params_nms)

farm_temp_files <- file.path(prepdata_path, "assigned_farms_alltemps") %>% 
  list.files(full.names = T) %>% 
  str_subset(paste0(fix_int(assigned_farms_TWN$farm_ID, 8), collapse = "|"))

farm_temp_IDs <- as.integer(str_extract(basename(farm_temp_files), "(?<=farm_)\\d+"))

farm_temps <- map(
  farm_temp_files,
  function(fnm) {
    qs_read(fnm) %>% 
      mutate(year = year(date)) %>% 
      group_by(year) %>% 
      reframe(
        mean = mean(sst),
        min = min(sst),
        max = max(sst)
      ) %>% 
      filter(year %in% c(2025, 2099))
  }
) %>% set_names(farm_temp_IDs)

farm_temps <- map(
  names(temp_responses),
  function(mn) {
    mn_IDs <- unique(assigned_farms_TWN$farm_ID[assigned_farms_TWN$model_name == mn])
    mn_temps <- farm_temps[which(names(farm_temps) %in% mn_IDs)]
    mn_temps %>% 
      bind_rows() %>% 
      group_by(year) %>% 
      reframe(
        mean = mean(mean),
        min = mean(min),
        max = mean(max)
      )
  }
) %>% setNames(names(temp_responses))

plots <- map2(
  temp_responses,
  names(temp_responses),
  function(tr, nm) {
    dat_tr <- data.frame(
      water_temp = water_temp,
      resp = tr
    )
    temps <- farm_temps[[nm]] %>% 
      mutate(FR = feeding_rate(water_temp = mean, species_params = species_params[[nm]]))
    
    dat_tr %>% 
      ggplot(aes(x = water_temp, y = resp)) +
      geom_line() +
      geom_vline(xintercept = temps$mean[temps$year == 2025], colour = "steelblue") +
      geom_vline(xintercept = temps$mean[temps$year == 2099], colour = "darkorange") +
      geom_point(data = filter(temps, year == 2025), aes(x = mean, y = FR), colour = "steelblue", size = 2.5) +
      geom_text(
        label = round(temps$FR[temps$year == 2025], 2), x = temps$mean[temps$year == 2025], y = temps$FR[temps$year == 2025], 
        colour = "steelblue", size = 3.5, vjust = -0.5, hjust = -0.25
      ) +
      geom_point(data = filter(temps, year == 2099), aes(x = mean, y = FR), colour = "darkorange", size = 2.5) +
      geom_text(
        label = round(temps$FR[temps$year == 2099], 2), x = temps$mean[temps$year == 2099], y = temps$FR[temps$year == 2099], 
        colour = "darkorange", size = 3.5, vjust = -0.5, hjust = -0.25
      ) +
      scale_y_continuous(limits = c(0, 1.05)) +
      prettyplot() +
      labs(
        title = str_to_sentence(str_replace(nm, "_", " ")),
        x = expression("Mean farm SST ("*degree*"C)"), 
        y = "Feeding response"
      )
  }
)

wrap_plots(plots)
ggsave(here("for-update-17-09-2026", "FR.png"), height = 6.5, width = 11)

# Metabolic responses -------------------------------------------------------------------------------
feed_params <- file.path(prepdata_path, "feed_parameters") %>% 
  list.files(full.names = T)

feed_names <- feed_params %>% 
  basename() %>% 
  str_remove_all(".qs")

feed_params <- map(feed_params, qs_read) %>% 
  setNames(feed_names)

feed_params <- list(
  "japanese_amberjack" = feed_params[["japanese_amberjack-japan"]],
  "milkfish"           = feed_params[["milkfish-philippines"]],
  "seabream"           = feed_params[["seabream-china"]],
  "tiger_prawn"        = feed_params[["tiger_prawn-world"]],
  "whiteleg_shrimp"    = feed_params[["whiteleg_shrimp-china"]]
)

temp_responses_2 <- map2(
  species_params, 
  feed_params,
  function(sp, fp) {
    get_species_responses(water_temp = water_temp, species_params = sp, feed_params = fp, ref_weight = sp[["CS"]]/2) %>% 
      select(water_temp, metab)
  }
) %>% setNames(species_params_nms)

plots_2 <- map2(
  temp_responses_2,
  names(temp_responses_2),
  function(tr, nm) {
    temps <- farm_temps[[nm]] %>% 
      mutate(
        MR = get_species_responses(
          water_temp = mean, 
          species_params = species_params[[nm]], 
          feed_params = feed_params[[nm]], 
          ref_weight = species_params[[nm]]["CS"]/2
        )$metab/1000
      )
    
    tr %>% 
      mutate(metab = metab/1000) %>% 
      ggplot(aes(x = water_temp, y = metab)) +
      geom_line() +
      geom_vline(xintercept = temps$mean[temps$year == 2025], colour = "steelblue") +
      geom_vline(xintercept = temps$mean[temps$year == 2099], colour = "darkorange") +
      geom_point(data = filter(temps, year == 2025), aes(x = mean, y = MR), colour = "steelblue", size = 2.5) +
      geom_text(
        label = round(temps$MR[temps$year == 2025], 2), x = temps$mean[temps$year == 2025], y = temps$MR[temps$year == 2025], 
        colour = "steelblue", size = 3.5, vjust = -0.5, hjust = -0.25
      ) +
      geom_point(data = filter(temps, year == 2099), aes(x = mean, y = MR), colour = "darkorange", size = 2.5) +
      geom_text(
        label = round(temps$MR[temps$year == 2099], 2), x = temps$mean[temps$year == 2099], y = temps$MR[temps$year == 2099], 
        colour = "darkorange", size = 3.5, vjust = -0.5, hjust = -0.25
      ) +
      # scale_y_continuous(limits = c(0, 1.05)) +
      prettyplot() +
      labs(
        title = str_to_sentence(str_replace(nm, "_", " ")),
        x = expression("Mean farm SST ("*degree*"C)"), 
        y = "Metabolic rate [kJ/g]"
      )
  }
)

wrap_plots(plots_2)
ggsave(here("for-update-17-09-2026", "MR.png"), height = 6.5, width = 11)

# Scatterplot of temperature? ----------------------------------------------------------------------------------------------------
farm_temp_files <- file.path(prepdata_path, "assigned_farms_alltemps") %>% 
  list.files(full.names = T) %>% 
  str_subset(paste0(fix_int(assigned_farms_TWN$farm_ID, 8), collapse = "|"))

farm_temp_IDs <- as.integer(str_extract(basename(farm_temp_files), "(?<=farm_)\\d+"))

# fnm <- farm_temp_files[1]
farm_temps <- map_dfr(
  farm_temp_files,
  function(fnm) {
    qs_read(fnm) %>% 
      mutate(year = year(date)) %>% 
      group_by(farm_ID, year) %>% 
      reframe(
        mean = mean(sst),
        min = min(sst),
        max = max(sst)
      ) %>% 
      filter(year %in% c(2025, 2062, 2099))
  }
) %>% 
  left_join(assigned_farms_TWN, by = "farm_ID")

xy_data <- map_dfr(
  unique(farm_temps$farm_ID),
  function(f) {
    temp <- filter(farm_temps, farm_ID == f)
    nm <- unique(temp$model_name)
    resp <- get_species_responses(
      water_temp = temp$mean, 
      species_params = species_params[[nm]],
      feed_params = feed_params[[nm]],
      ref_weight = species_params[[nm]]["CS"]/2
    )
    temp %>% 
      mutate(FR_mean = resp$feeding_response)
  }
) %>% 
  select(farm_ID, year, model_name, FR_mean) %>% 
  pivot_wider(names_from = year, names_prefix = "yr_", values_from = FR_mean)

xy_data %>% 
  ggplot(aes(x = yr_2025, y = yr_2099, colour = model_name)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  geom_point(size = 3.5) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(limits = c(0.5,1)) +
  scale_y_continuous(limits = c(0.5,1)) +
  prettyplot() +
  theme(legend.position = "right")

