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
