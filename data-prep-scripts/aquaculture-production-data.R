# Description ----------------------------------------------------------------------------------------------
# This document tidies up historical production data from the FAO, ready to be used in assessing which species will be chosen for future analysis.

library(tidyverse)
library(here)
library(qs2)
library(janitor)
library(units)
library(fuzzyjoin)

here("src", "dirs.R") %>% source()
here("src", "functions.R") %>% source()

# Get aquaculture quantity in tonnes live weight
quantity <- file.path(bigdata_path, "Aquaculture_2025.1.0", "Aquaculture_Quantity.csv") %>% 
  read.csv() %>% 
  clean_names() %>% 
  mutate(
    value = set_units(value, "t"),
    across(where(is.character), str_trim)
  ) %>% 
  select(-measure) # all values are quantity_liveweight

# Apply environment codes and filter by "marine"
environment_codes <- file.path(bigdata_path, "Aquaculture_2025.1.0", "CL_FI_PRODENVIRONMENT.csv") %>% 
  read.csv() %>% 
  distinct(Code, Name_En) %>% 
  mutate(across(where(is.character), str_trim))

quantity <- quantity %>% 
  left_join(environment_codes, by = c("environment_alpha_2_code" = "Code")) %>%
  rename(environment = Name_En) %>% 
  filter(environment == c("Marine", "Brackishwater"))

# Apply species codes — species_name and group_name are looked up from FAO_name_keys.csv
# using partial (grepl-style) matching on Scientific_Name
species_codes <- file.path(bigdata_path, "Aquaculture_2025.1.0", "CL_FI_SPECIES_GROUPS.csv") %>%
  read.csv() %>%
  distinct(X3A_Code, Scientific_Name) %>% 
  mutate(across(where(is.character), str_trim))

quantity <- quantity %>% 
  left_join(species_codes, by = c("species_alpha_3_code" = "X3A_Code"))

name_keys <- here("data", "raw_data", "FAO_name_keys.csv") %>%
  read.csv()

quantity <- quantity %>%
  regex_left_join(name_keys, by = c("Scientific_Name" = "Sci_Name_Match")) %>% 
  mutate(
    species_name = as.factor(species_name),
    group_name = as.factor(group_name)
  )

# Check that all species have a name and group
quantity %>% 
  filter(is.na(species_name)) %>% 
  arrange(-value) %>% 
  pull(Scientific_Name) %>% 
  unique()

quantity %>% 
  filter(is.na(group_name)) %>% 
  distinct(group_name, species_name)

# Apply status codes
status_codes <- file.path(bigdata_path, "Aquaculture_2025.1.0", "CL_FI_SYMBOL_SDMX.csv") %>% 
  read.csv() %>% 
  distinct(Symbol, Name_En, Description_En) %>% 
  mutate(across(where(is.character), str_trim))

quantity <- quantity %>% 
  left_join(status_codes, by = c("status" = "Symbol")) %>% 
  rename(status_nm = Name_En,
         status_description = Description_En) %>% 
  filter(!status %in% c("L", "N", "M", "O", "Q")) %>% 
  select(-status, -status_description)

# Apply country and area codes
country_codes <- file.path(bigdata_path, "Aquaculture_2025.1.0", "CL_FI_COUNTRY_GROUPS.csv") %>% 
  read.csv() %>% 
  select(UN_Code, Identifier, Name_En, ISO3_Code) %>% 
  mutate(across(where(is.character), str_trim))
area_codes <- file.path(bigdata_path, "Aquaculture_2025.1.0", "CL_FI_WATERAREA_GROUPS.csv") %>% 
  read.csv() %>% 
  select(Code, Name_En) %>% 
  mutate(across(where(is.character), str_trim))

quantity <- quantity %>% 
  left_join(country_codes, by = c("country_un_code" = "UN_Code")) %>% 
  left_join(area_codes, by = c("area_code" = "Code")) %>% 
  rename(
    country = Name_En.x, 
    fao_fishing_area = Name_En.y,
    country_id = Identifier
    )

# Final cleanup
quantity <- quantity %>% 
  mutate(
    environment = as.factor(environment),
    Scientific_Name = as.factor(Scientific_Name),
    status_nm = as.factor(status_nm),
    country = as.factor(country),
    fao_fishing_area = as.factor(fao_fishing_area),
    ISO3_Code = as.factor(ISO3_Code)
  ) %>% 
  select(-c(country_un_code, species_alpha_3_code, area_code, environment_alpha_2_code))

# Clean up countries that don't exist/have been renamed
quantity <- quantity %>% 
  filter(!ISO3_Code %in% c("YUG", "SUN", "SCG", "ANT")) %>%
  mutate(
    country = case_when(
      ISO3_Code == "EAZ" ~ "United Republic of Tanzania",
      T ~ country
    ),
    ISO3_Code = case_when(
      ISO3_Code == "EAZ" ~ "TZA",
      T ~ ISO3_Code
    ))

# Save two versions - one marine only, one marine + brackish
quantity %>% 
  filter(environment != "Brackishwater") %>% 
  qd_save(file.path(prepdata_path, "FAO_aquaculture_quantity_MA.qs"))

quantity %>% 
  qd_save(file.path(prepdata_path, "FAO_aquaculture_quantity_MA_BR.qs"))
