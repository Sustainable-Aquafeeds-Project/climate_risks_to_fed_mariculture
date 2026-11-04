library(qs2)
library(tidyverse)
library(here)

here("src", "dirs.R") %>% source()
here("src", "other_functions.R") %>% source()
here("src", "model_functions.R") %>% source()

# The purpose of this script is to identify and fix any discrepancies between the farm production files that have already been run and any changes to farm numbers/model assignments when the location assignment protocol is run again. Here's why this might be needed.

# When the location assignment protocol is run, it feeds the temperature response curves of each species into its optimisation algorithm and assigns existing ("mapped") farms to species as needed. If mapped farms don't exist for a certain species, it creates them based on thermal niche suitability ("extracted").

# If a species model's temperature response parameters change, this could change the model_names assigned to farm_IDs within any country where that species is grown. Additionally, if any changes are made to "extracted" farms assigned to that species, new farm_ID numbers will have been generated. 

assigned_farms_old <- file.path(prepdata_path, "assigned_farms_OLD.qs") %>% qs_read()
assigned_farms_new <- file.path(prepdata_path, "assigned_farms.qs") %>% qs_read()

# Get all existing production files - if you want to re-do production files for any reason other than assignment changes, make sure they're deleted (or moved) before this step.
existing_prod_files <- file.path(outs_path, "data", "farm_production") %>% 
  list.files(recursive = T, full.names = T) %>% 
  str_subset("MC5000.qs") %>% 
  str_subset("RESCALED", negate = T)

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

sprintf("Files to remove, no longer assigned: %d", length(to_del_ints))
# file.remove(existing_prod_files[to_del_ints])
file.rename(existing_prod_files[to_del_ints], file.path(outs_path, "data", "farm_production", "temp_hold", basename(existing_prod_files[to_del_ints])))

# If the farm has NA in the model_name_old column, it is new and can be ignored, although the production will need to be run again if it exists already.
new_mapped <- assigned_farms_mapped %>% 
  filter(is.na(model_name_old)) %>% 
  pull(model_name_new) %>% 
  unique()
sprintf("New files to run: %d", length(new_mapped))

# If the values in the model_name_old and model_name_new columns do not match, any existing file should be deleted.
to_del_IDs <- assigned_farms_mapped %>% 
  mutate(match = model_name_old == model_name_new) %>% 
  filter(!match) %>% 
  pull(farm_ID)
to_del_ints <- which(existing_prod_IDs %in% to_del_IDs)

sprintf("Files to remove, new species assigned: %d", length(to_del_ints))
# file.remove(existing_prod_files[to_del_ints])
file.rename(existing_prod_files[to_del_ints], file.path(outs_path, "data", "farm_production", "temp_hold", basename(existing_prod_files[to_del_ints])))

# Mapped locations have a maximum farm_ID of 95443. Any farm with a higher number is extracted.
assigned_farms_new %>% 
  filter(source != "mapped") %>% 
  nrow()
# Rather than trying to match locations, I'm instead going to make the decision that any existing production files from non-mapped farms should be deleted.
to_del_ints <- which(existing_prod_IDs > 95443)
sprintf("Files to remove, extracted location: %d", length(to_del_ints))
file.rename(existing_prod_files[to_del_ints], file.path(outs_path, "data", "farm_production", "temp_hold", basename(existing_prod_files[to_del_ints])))

# CHECK - make sure there are no existing production files where:
# 1. The farm_ID does not exist at all in assigned_farms_new
# 2. The model_name in assigned_farms_old is different to that in assigned_farms_new
# 3. The farm_ID is greater than 95443

existing_prod_files <- file.path(outs_path, "data", "farm_production") %>% 
  list.files(recursive = T, full.names = T) %>% 
  str_subset("MC5000.qs") %>% 
  str_subset("testing", negate = T) %>% 
  str_subset("temp_hold", negate = T) %>% 
  str_subset("RESCALED", negate = T)
existing_prod_IDs <- existing_prod_files %>% 
  basename() %>% 
  str_extract("\\d+(?=\\_MC5000.qs$)") %>% 
  as.integer()

any(!existing_prod_IDs %in% assigned_farms_new$farm_ID) # expecting FALSE
assigned_farms_mapped %>% 
  filter(farm_ID %in% existing_prod_IDs) %>% 
  filter(model_name_new != model_name_old)              # expecting 0 rows
existing_prod_IDs[existing_prod_IDs > 95443]            # expecting 0

# NEXT STEP - Now we need to make sure that the remaining files (which have the correct species, location, etc) have the correct scaling factor. This is because the number of farms and the production they're trying to match may have changed.

existing_prod <- data.frame(
  farm_ID = existing_prod_IDs,
  filename = existing_prod_files
)
assigned_farms_new <- right_join(assigned_farms_new, existing_prod, by = "farm_ID") # Only need to rescale farms where production file exists
assigned_farms_new <- right_join(
  assigned_farms_old %>% 
    dplyr::select(farm_ID, prod_perfarm) %>% 
    st_drop_geometry() %>% 
    rename(prod_perfarm_old = prod_perfarm),
  assigned_farms_new,
  by = "farm_ID"
) %>% 
  filter(prod_perfarm != prod_perfarm_old) %>%                                       # Only need to rescale farms where prod_perfarm has changed
  mutate(prod_perfarm_change = prod_perfarm/prod_perfarm_old)

walk(
  1:nrow(assigned_farms_new),
  function(i) {
    
    data_all  <- qs_read(assigned_farms_new$filename[1])
    vars <- c("weight_scaled", "dw_scaled", "ing_act_scaled", "total_excr_scaled", "food_prov_scaled")
    prod_perfarm_change <- assigned_farms_new$prod_perfarm_change[i]
    
    # dat <- data_all[[2]]
    data_all_rescaled <- map(
      data_all, 
      function(dat) {
        for (var in vars) {
          vardat     <- dat[["stats"]][[var]]                # get both columns
          vardat[,2] <- vardat[,2]/vardat[,1]                # express sd as proportion
          vardat[,1] <- vardat[,1] * prod_perfarm_change     # rescale mean
          vardat[,2] <- vardat[,2] * prod_perfarm_change     # rescale sd
          vardat[,2] <- vardat[,2] * vardat[,1]              # re-express sd as absolute
          dat[["stats"]][[var]] <- vardat                    # put data back
        }
        dat
      }
    )

    qs_save(data_all_rescaled, assigned_farms_new$filename[1])

  },
  .progress = T
)

