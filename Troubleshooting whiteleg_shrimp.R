# Troubleshooting whiteleg_shrimp

# This is a simplified version of the process detailed in documentation-qmds\05_01_run_production.qmd.
# I created it because whenever I tried to run the whiteleg_shrimp model line by line it runs fine, but running it in the full pipeline gives NAs for all years. There's an error somewhere, I just can't find it.
# This problem only seems to happen for whiteleg shrimp.

library(tidyverse)
library(qs2)
library(purrr)

source("src/model_functions.R")
source("src/other_functions.R")

focus_ID <- 23

# Get species_parameters
species_params <- "C:/Users/treimer/Documents/R-temp-files/climate_risks_to_fed_mariculture/data/prepped_data/species_parameters/species_parameters_whiteleg_shrimp.qs" %>% 
  qd_read()

# Get assigned farm information
assigned_farms <- "C:/Users/treimer/Documents/R-temp-files/climate_risks_to_fed_mariculture/data/prepped_data/assigned_farms.qs" %>% 
  qs_read() %>% 
  sf::st_drop_geometry() %>% 
  filter(farm_ID == focus_ID)

assigned_farms$filename <- "C:/Users/treimer/Documents/R-temp-files/climate_risks_to_fed_mariculture/data/prepped_data/assigned_farms_alltemps" %>% 
  list.files(full.names = T) %>% 
  str_subset(fix_int(focus_ID, 8))

# Load feed parameters
feed_params <- "C:/Users/treimer/Documents/R-temp-files/climate_risks_to_fed_mariculture/data/prepped_data/feed_parameters/whiteleg_shrimp-ecuador.qs" %>% 
  qs_read()

# Get the best "growing window" for this species
get_growing_window <- function(temperatures, grow_days, species_params) {
  temperatures$fr <- sapply(
    X = temperatures$sst, 
    FUN = feeding_rate, 
    species_params = species_params
  )
  temperatures$fr[is.na(temperatures$fr)] <- 0

  temperatures %>% 
    mutate(roll_fr = slide_dbl(fr, sum, .before = 0, .after = grow_days - 1)) %>%
    slice_max(roll_fr, n = 1) %>%
    mutate(end_doy = doy + grow_days - 1) %>% 
    select(start = doy, end = end_doy)
}

gd <- species_params[["gd"]] %>% as.integer()

# Get the base temperatures for this farm
temp_2025 <- assigned_farms$filename %>% 
  qs_read() %>% 
  filter(year(date) == 2025) %>% 
  mutate(doy = as.integer(format(date, "%j")))

# Repeat the base year (2025) to set growing conditions
temp_2025 <- rbind(
  temp_2025,
  temp_2025 %>% mutate(doy = doy + 365),
  temp_2025 %>% mutate(doy = doy + 730)
)

doys <- get_growing_window(
  temperatures = temp_2025,
  grow_days = gd,
  species_params = species_params
  ) %>% 
  slice_head(n = 1)

# Narrow down the temperature timeseries to only the growing days
temp_base <- temp_2025 %>% 
  filter(doy >= doys$start & doy <= doys$end) %>% 
  pull(sst)

# Grow a representative individual fish with this timeseries
rep_ind <- fish_growth(
  species_params = species_params, 
  water_temp = temp_base, 
  feed_params = feed_params, 
  times = c(t_start = doys$start, t_end = doys$end, dt = 1L), 
  init_weight = species_params[["meanW"]], 
  ingmax = species_params[["meanImax"]], 
  output_vars = c("weight")
)
rep_ind_weight <- unname(rep_ind[nrow(rep_ind), 2]) # get harvest weight (g)

# Get scaling factor from harvest weight (g) and per-farm production (t)
pop <- (assigned_farms$prod_perfarm * 10^6)/rep_ind_weight

# At the end, give the start and end days, the individual harvest size, and the final population (scaling factor)
assigned_farms <- assigned_farms %>% 
  mutate(
    start_doy = doys$start,
    end_doy = doys$end,
    harvest_size = rep_ind_weight,
    final_pop = pop
  )

# RUN PRODUCTION
# Get all temperatures at this farm (years 2025-2099)
temps_all <- assigned_farms$filename %>% 
  qs_read() %>% 
  mutate(
    doy = as.integer(format(date, "%j")),
    year = year(date)
  ) %>% 
  select(-date)

# Fill in 2024 to reconcile multi-year growing periods
temps_all <- rbind(
    temps_all %>% filter(year == 2025) %>% mutate(year = 2024),
    temps_all
  )

# Get correct species growtimes
grow_times <- c(
  t_start = assigned_farms$start_doy,
  t_end = assigned_farms$end_doy,
  dt = 1L
)

# Filter to actual growing window for each harvest year
# The previous year is added so that the harvest happens IN the year specified
harvest_years <- 2025:2099
temps_all <- map(
  harvest_years,
  function(yr) {
    if (grow_times[["t_end"]] > 365) {
      temp <- temps_all %>% filter(year %in% c(yr, yr-1))
      temp$doy <- 1:nrow(temp)
    } else {
      temp <- temps_all %>% filter(year == yr)
    }
    temp %>% 
      filter(doy >= grow_times[["t_start"]] & doy <= grow_times[["t_end"]]) %>% 
      mutate(act_year = year, year  = yr)
  }
) %>% 
  setNames(harvest_years)

# Generate the population (scaling factor)
pop <- generate_pop(harvest_n = assigned_farms$final_pop, mort = species_params[["mortmyt"]], times = grow_times)

# NEW - check full production
check_farm_growth_inputs(
  species_params = species_params, 
  feed_params = feed_params, 
  water_temp = temps_all[["2025"]]$sst, 
  times = grow_times,
  N_pop = pop, 
  use_MC_population = T, 
  MC_pop = 5000
)

# Run the full production
fg <- map(
  temps_all,
  function(temps) {
    farm_growth(
      species_params = species_params, 
      feed_params = feed_params, 
      water_temp = temps$sst, 
      times = grow_times,
      N_pop = pop, 
      use_MC_population = T, 
      MC_pop = 100
    )
  },
  .progress = T
)

# Have a quick look at the results
fg_weight <- tidy_stat(fg, "weight") #%>% 
  # mutate(
  #   year = case_when(prod_day <= 24 ~ year - 1, T ~ year),
  #   day  = case_when(prod_day <= 24 ~ day - 1, T ~ day - 365),
  #   date = as.Date(paste(year-1, day), format = "%Y %j")
  # )

ggplot(fg_weight, aes(x = day, y = mean)) + 
  geom_line() +
  facet_wrap(~year)
