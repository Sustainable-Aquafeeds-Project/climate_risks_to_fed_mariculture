library(rnaturalearth)
library(rnaturalearthhires)
library(tidyverse)
library(qs2)
library(terra)

here("src", "dirs.R") %>% source()
here("src", "functions.R") %>% source()

world_land <- ne_countries(scale = "large", returnclass = "sf") %>% 
  select(name_long , adm0_a3)

EEZ_by_FAO_full <- file.path(prepdata_path, "EEZ_data_with_FAO_regions.qs") %>%
  qd_read() %>%
  st_transform(crs = 4326) %>%
  st_make_valid()

assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs_read() %>% 
  mutate(farm_ID = as.integer(farm_ID))

all_prod_files <- file.path(outs_path, "data", "farm_production") %>% 
  list.files(full.names = T, pattern = "MC5000\\.qs$")
all_farm_IDs <- as.integer(str_extract(basename(all_prod_files), "(?<=farmID_)\\d+"))

# For each iso
iso <- "NOR"

this_farm_locs <- assigned_farms %>% filter(ISO3_Code == iso)
this_EEZ_vect <- EEZ_by_FAO_full %>% filter(ISO3_Code == iso)
this_farm_IDs <- this_farm_locs %>% st_drop_geometry() %>% pull(farm_ID) %>% as.integer()
this_land <- world_land %>% filter(adm0_a3 == iso)

this_data <- all_prod_files[match(this_farm_IDs, all_farm_IDs)] %>% 
  map_dfr(
    function(fnm) {
      qs_read(fnm) %>% 
        tidy_stat("weight_scaled") %>% 
        mutate(farm_ID = as.integer(str_extract(basename(fnm), "(?<=farmID_)\\d+"))) %>%
        filter(prod_day == max(prod_day))
    }, 
  .progress = T
  )

this_data <- this_data %>% 
  filter(year == 2075) %>% 
  left_join(this_farm_locs, by = join_by(farm_ID))

ggplot() +
  geom_sf(data = this_land, fill = "grey80", color = NA) +
  geom_sf(data = this_EEZ_vect, aes(color = F_CODE), fill = NA, linewidth = 0.7) +
  geom_sf(data = st_as_sf(this_data), aes(fill = mean), shape = 21, color = "black") +
  facet_wrap(~ model_name) +
  coord_sf() +
  prettyplot()



