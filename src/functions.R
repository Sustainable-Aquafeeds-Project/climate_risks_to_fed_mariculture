library(tidyverse)

sumna <- function(x) {sum(x, na.rm = T)}
maxna <- function(x) {max(x, na.rm = T)}
minna <- function(x) {min(x, na.rm = T)}
meanna <- function(x) {mean(x, na.rm = T)}
medianna <- function(x) {median(x, na.rm = T)}
sdna <- function(x) {sd(x, na.rm = T)}
rangena <- function(x) {range(x, na.rm = T)}

# Adds leading zeros to integers for nicer filenames
fix_int <- function(n, digits = 4) {
  vapply(n, function(x) {
    stringr::str_flatten(c(rep("0", digits-nchar(as.character(x))), as.character(x)))
  }, character(1))
}

# This takes a list of files and groups them by DOY (month-day)
by_doy <- function(files) {
  tibble(filenames = files) %>%
    mutate(
      date = str_extract(filenames, "\\d{4}-\\d{2}-\\d{2}") %>% as.Date(),
      doy = yday(date) %>% as.integer(),
      month_day = format(date, "%m-%d"),
      leap_year = leap_year(date)
    ) %>%
    mutate(
      # For leap years, Feb 29 gets converted to Feb 28
      month_day = if_else(month_day == "02-29", "02-28", month_day),
      doy_adjusted = case_when(
        leap_year & doy >= 60 ~ doy-1,
        T ~ doy
      ) %>% as.integer()
    ) %>%
    group_by(doy_adjusted, month_day) %>%
    summarise(
      filenamess = list(filenames),
      .groups = "drop"
    ) %>%
    arrange(doy_adjusted)
}

# Helper: returns the meanyear layer name(s) and weight(s) for a given date.
# For Feb 29: returns DOY_059 (Feb 28) and DOY_060 (Mar 1) with weight 0.5 each.
# For all other dates in a leap year after Feb 28: shifts DOY back by 1.
# For all non-leap-year dates: direct lookup.
get_doy_lookup <- function(date) {
  
  if (leap_year(date) && month(date) == 2 && day(date) == 29) {
    return(list(
      layers  = c("DOY_059", "DOY_060"),
      weights = c(0.5, 0.5)
    ))
  }
  
  doy <- yday(date)
  
  # In a leap year, Mar 1 onwards is DOY 61+; shift back by 1 to align with
  # the 365-day meanyear (where Mar 1 = DOY_060, Dec 31 = DOY_365)
  if (leap_year(date) && doy > 60) doy <- doy - 1
  
  list(
    layers  = paste0("DOY_", str_pad(doy, 3, pad = "0")),
    weights = 1
  )
}

# The function extracting a single location is intended for mapped farms, which have IDs. These will be saved for use later.
extract_point_location_temperature <- function(
  location_df,               # dataframe containing the geometry (sf_point) for the location needed
  global_raster_stack,       # raster stack of global daily temperatures in a typical year
  save_file = T,
  out_path = NA,
  overwrite = F
) {

  out_file <- file.path(out_path, paste0("meanyear_farmID_", fix_int(location_df$farm_ID, 8), ".qs"))
  
  if (!file.exists(out_file) | overwrite | !save_file) {
    # Get coordinates and resolve cell number directly - avoids terra::vect() overhead
    coords   <- sf::st_coordinates(location_df)
    cell_num <- terra::cellFromXY(global_raster_stack, coords)

    # Extract at the exact cell (cell-number extraction returns values only, no ID column)
    temps <- terra::extract(global_raster_stack, cell_num) %>% as.vector() %>% unlist()

    # If any layers are NA (e.g. point falls on a land mask), fall back to the mean of surrounding cells
    if (any(is.na(temps))) {
      cell_nums <- c(cell_num, terra::adjacent(global_raster_stack, cell_num, directions = 8))

      # Get mean of all cell temps
      temps <- terra::extract(global_raster_stack, cell_nums)
      temps <- colMeans(temps, na.rm = T)

      # If this isn't enough, expand out in the longitude direction to minimise latitude impacts on temperature (up to ~0.5 degrees away)
      if (all(is.na(temps))) {
        # message("\n Had to extend search to 0.5 degrees away for farm_ID ", fix_int(location_df$farm_ID, 5), ".")
        # Get row and col numbers for 9 focus cells
        rc   <- terra::rowColFromCell(global_raster_stack, cell_nums)
        cols <- unique(rc[, 2])
        rows <- unique(rc[, 1])
        # Extend to 1 more cell in E and W directions (9 + 6 = 15)
        new_rows    <- unique(c(rows, rows - 1, rows + 1))
        new_rowcols <- expand.grid(row = new_rows, col = cols) %>% distinct()
        new_cells   <- terra::cellFromRowCol(global_raster_stack, new_rowcols$row, new_rowcols$col)

        # Get temps of expanded cell list
        cell_nums <- c(cell_nums, new_cells) %>% unique()
        temps     <- terra::extract(global_raster_stack, cell_nums)
        temps     <- colMeans(temps, na.rm = T)

        # If this STILL isn't enough, expand out in the longitude direction up to 1 degree away
        if (all(is.na(temps))) {
          # message("\n Had to extend search to 1.0 degree away for farm_ID ", fix_int(location_df$farm_ID, 5), ".")
          # Get row and col numbers for expanded cells
          new_rows    <- unique(c(new_rowcols$row, new_rowcols$row - 1, new_rowcols$row + 1, new_rowcols$row - 2, new_rowcols$row + 2))
          new_rowcols <- expand.grid(row = new_rows, col = cols) %>% distinct()
          new_cells   <- terra::cellFromRowCol(global_raster_stack, new_rowcols$row, new_rowcols$col)

          # Get temps of expanded cell list
          cell_nums <- c(cell_nums, new_cells) %>% unique()
          temps     <- terra::extract(global_raster_stack, cell_nums)
          temps     <- colMeans(temps, na.rm = T)

          # If this STILL isn't enough, just report it
          if (all(is.na(temps))) {
            message("\n WARNING: temperatures for farm_ID ", fix_int(location_df$farm_ID, 5), " still blank with 1-degree buffer.")
          }
        }
      }
    }

    temps <- tibble(doy = 1:length(temps), sst = temps)
    if (save_file) {
      qd_save(temps, out_file)
    } else {
      return(temps)
    }
  }
}

# For the thermal niche suitability
get_species_responses <- function(water_temp, species_params, feed_params, ref_weight) {
  fr             <- sapply(water_temp, feeding_rate, species_params = species_params)
  T_response     <- exp(species_params['pk'] * water_temp)
  feed_ingested  <- unname(species_params['meanImax'] * (ref_weight^species_params['m']) * fr)

  assim <- vapply(
    names(feed_params),
    function(nm) {
      fp <- feed_params[[nm]]
      apportion_feed_v(feed_ingested, feed_ingested, fp$proportion, fp$macro, fp$digestibility)[['assimilated']]
    },
    numeric(length(feed_ingested))
  )
  protein       <- assim[, "protein"] * species_params['epsprot']
  lipid         <- assim[, "lipid"]   * species_params['epslip']
  carb          <- assim[, "carb"]    * species_params['epscarb']
  E_assim_total <- as.numeric(protein + lipid + carb)
  
  anab      <- E_assim_total * (1 - species_params['alpha'])
  catab     <- species_params['epsO2'] * species_params['k0'] * T_response * (ref_weight^species_params['n']) * species_params['omega'] # catabolic response
  E_somat   <- rep(species_params['a'] * ref_weight^species_params['k'], length(water_temp))       # energy content of somatic tissue at half harvest weight
  dw        <- (anab - catab) / E_somat
  
  data.frame(
    water_temp = water_temp,
    feeding_response = fr,
    feed_ingested = feed_ingested/ref_weight,
    E_assim = E_assim_total, # total energy assimilated
    anabolism = anab, # energy gained from feeding (minus the energy it took to process the feed)
    catabolism = catab,
    metab = anab - catab,
    somatic_energy = E_somat, # energy content of somatic tissue at half harvest weight
    weight_change = dw/ref_weight # weight change at this temperature when at half harvest weight and feed to satiety
  )
}

# Wrapper to use the above function on a raster stack of temperatures instead of a single vector
get_species_responses_r <- function(
  temp_stack, species_params, feed_params, ref_weight, 
  out_cols = c("feeding_response", "feed_ingested", "E_assim", "anabolism", "catabolism", "metab", "somatic_energy", "weight_change")
) {
  ncell <- terra::ncell(temp_stack)
  nlyr  <- terra::nlyr(temp_stack)
  temp_vec <- as.vector(terra::values(temp_stack))

  # Only run the model on real temperatures; cells masked outside the EEZ/land are NA
  ok   <- is.finite(temp_vec)
  resp <- get_species_responses(temp_vec[ok], species_params, feed_params, ref_weight)

  result <- lapply(out_cols, function(col) {
    full <- rep(NA_real_, length(temp_vec))      # scatter results back to full length
    full[ok] <- resp[[col]]
    r <- terra::setValues(temp_stack,            # reuse the stack purely as a geometry template
                          matrix(full, nrow = ncell, ncol = nlyr))
    names(r) <- names(temp_stack)                # keep the layer names
    r
  })
  names(result) <- out_cols
  result
}


