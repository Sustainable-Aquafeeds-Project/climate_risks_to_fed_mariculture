library(qs2)
library(tidyverse)
library(here)

here("src", "dirs.R") %>% source()
here("src", "functions.R") %>% source()
here("src", "model_functions.R") %>% source()

# The purpose of this script is to identify and fix any discrepancies between the farm production files that have already been run and any changes to farm numbers/model assignments when the location assignment protocol is run again. Here's why this might be needed.

# When the location assignment protocol is run, it feeds the temperature response curves of each species into its optimisation algorithm and assigns existing ("mapped") farms to species as needed. If mapped farms don't exist for a certain species, it creates them based on thermal niche suitability ("extracted").

# If a species model's temperature response parameters change, this could change the model_names assigned to farm_IDs within any country where that species is grown. Additionally, if any changes are made to "extracted" farms assigned to that species, new farm_ID numbers will have been generated. 

assigned_farms_old <- file.path(prepdata_path, "assigned_farms_old.qs") %>% qs_read() %>% 
  mutate(model_name = case_when(model_name == "japanese_seabass" ~ "asian_seabass", T ~ model_name))
assigned_farms_new <- file.path(prepdata_path, "assigned_farms.qs") %>% qs_read()

existing_prod_files <- file.path(outs_path, "data", "farm_production") %>% 
  list.files(recursive = T, full.names = T) %>% 
  str_subset("MC5000.qs")
existing_prod_IDs <- existing_prod_files %>% 
  basename() %>% 
  str_extract("\\d+(?=\\_MC5000.qs$)") %>% 
  as.integer()

# Do "mapped" locations first
assigned_farms_mapped <- full_join(
  assigned_farms_new %>% 
    filter(source == "mapped") %>% 
    st_drop_geometry() %>% 
    select(ISO3_Code, F_CODE, farm_ID, model_name) %>% 
    rename(model_name_new = model_name),
  assigned_farms_old %>% 
    filter(source == "mapped") %>% 
    st_drop_geometry() %>% 
    select(ISO3_Code, F_CODE, farm_ID, model_name) %>% 
    rename(model_name_old = model_name),
  by = join_by(ISO3_Code, F_CODE, farm_ID)
)

# First, we're going to check if any of the production files that already exist have changed in their assigned status or assigned species.
# If the farm has NA in the model_name_new column, it is no longer assigned and any existing files should be deleted.
to_del_IDs <- assigned_farms_mapped$farm_ID[is.na(assigned_farms_mapped$model_name_new)]
to_del_ints <- which(existing_prod_IDs %in% to_del_IDs)
file.remove(existing_prod_files[to_del_ints])
# If the farm has NA in the model_name_old column, it is new and can be ignored, although the production will need to be run again.
assigned_farms_mapped %>% 
  filter(is.na(model_name_old)) %>% 
  pull(model_name_new) %>% 
  unique()

# If the values in the model_name_old and model_name_new columns do not match, any existing file should be deleted.
to_del_IDs <- assigned_farms_mapped %>% 
  mutate(match = model_name_old == model_name_new) %>% 
  filter(!match) %>% 
  pull(farm_ID)

to_del_ints <- which(existing_prod_IDs %in% to_del_IDs)
file.remove(existing_prod_files[to_del_ints])

# Mapped locations have a maximum farm_ID of 95443. There are only 309 non-mapped farms.
assigned_farms_new %>% 
  filter(source != "mapped") %>% 
  nrow()
# Rather than trying to match locations, I'm instead going to make the decision that any existing production files from non-mapped farms should be deleted.
to_del_ints <- which(existing_prod_IDs > 95443)
file.remove(existing_prod_files[to_del_ints])

# CHECK - make sure there are no existing production files where:
# 1. The farm_ID does not exist at all in assigned_farms_new
# 2. The model_name in assigned_farms_old is different to that in assigned_farms_new
# 3. The farm_ID is greater than 95443

existing_prod_files <- file.path(outs_path, "data", "farm_production") %>% 
  list.files(recursive = T, full.names = T) %>% 
  str_subset("MC5000.qs")
existing_prod_IDs <- existing_prod_files %>% 
  basename() %>% 
  str_extract("\\d+(?=\\_MC5000.qs$)") %>% 
  as.integer()

any(!existing_prod_IDs %in% assigned_farms_new$farm_ID) # expecting FALSE
assigned_farms_mapped %>% 
  filter(farm_ID %in% existing_prod_IDs) %>% 
  filter(model_name_new != model_name_old)              # expecting 0 rows
existing_prod_IDs[existing_prod_IDs > 95443]            # expecting 0


