library(tidyverse)
library(qs2)

assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs_read() %>% 
  filter(ISO3_Code == "FIN")
farm_IDs <- unique(assigned_farms$farm_ID) %>% sort()

species_params <- file.path(prepdata_path, "species_parameters", "species_parameters_rainbow_trout.qs") %>% qd_read()
feed_params <- file.path(prepdata_path, "feed_parameters", "rainbow_trout-norway.qs") %>% qs_read()

temp_files <- list.files(file.path(prepdata_path, "assigned_farms_alltemps"), full.names = T) %>% 
  str_subset(paste0(fix_int(farm_IDs, 8), collapse = "|")) %>% 
  sort() %>% 
  setNames(farm_IDs)

tf <- temp_files[1]
FR_df <- map_dfr(
  temp_files,
  function(tf) {
    tf %>% 
      qs_read() %>% 
      mutate(
        year = year(date),
        FR = feeding_rate(water_temp = sst, species_params = species_params)
      ) %>% 
      group_by(farm_ID, year) %>% 
      reframe(total_FR = sumna(FR)/365)
  },
  .progress = T
)

FR_df %>% 
  ggplot(aes(x = year, y = total_FR)) +
  geom_smooth()

MR_df <- map_dfr(
  temp_files,
  function(tf) {
    tdf <- tf %>% 
      qs_read() %>% 
      mutate(year = year(date))

    tdf$MR <- as.vector(unlist(map(.x = tdf$sst, .f = get_metabolism, species_params = sp, feed_params = feed_params,rep_size = sp[["CS"]])))

     tdf %>% 
      group_by(farm_ID, year) %>% 
      reframe(total_MR = sumna(MR)/365)
  },
  .progress = T
)

MR_df %>% 
  ggplot(aes(x = year, y = total_MR)) +
  geom_smooth()

