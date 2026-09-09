library(tidyverse)
library(tidyverse)
library(ozmaps)
library(paletteer)
library(patchwork)

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
  fr             <- feeding_rate(water_temp, species_params)  # already vectorised
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

# Shortcut plotting functions
prettyplot <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "serif", size = 12, colour = "black"),
      legend.position = "none",
      axis.title.y = element_text(vjust = 1.5),
      axis.title.x = element_text(vjust = 1.5),
      legend.title = element_blank()
    )
}

rm_y_axis <- function() {
  theme(
    axis.title.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.line = element_blank()
  )
}

rm_x_axis <- function() {
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.line = element_blank()
  )
}

plot_sim_results <- function(df, CS = NA) {  
  # --- helper: base theme ---
  base_theme <- function(d) {
    ggplot(d, aes(days, value)) +
      theme_classic(base_size = 10) +
      theme(
        legend.position  = "bottom", 
        legend.title     = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_blank()
      ) +
      labs(x = "Day")
  }
  
  # --- single-variable panels ---
  p_weight <- base_theme(filter(df, output_var == "weight")) +
    geom_line() + labs(title = "Individiaul weight (g)", y = NULL) 
  
  if (!is.na(CS)) {
    p_weight <- p_weight +
      geom_hline(yintercept = CS, linewidth = 0.75, colour = "red", linetype = "dashed") +
      annotate("text", x = min(df$days), y = CS, label = sprintf("CS: %d", CS), hjust = -0.1, vjust = -0.5, colour = "red")
  }
  
  p_dw <- base_theme(filter(df, output_var == "dw")) +
    geom_line() + labs(title = "Daily weight change (g/d)", y = NULL)
  
  p_sgr <- base_theme(filter(df, output_var == "SGR")) +
    geom_line() + labs(title = "SGR (%)", y = NULL)
  
  # --- energy: overlaid lines ---
  p_energy <- base_theme(
    df %>% 
      filter(output_var %in% c("E_assim", "E_somat", "metab"))
    ) +
    geom_line(aes(colour = output_var)) +
    labs(title = "Energy (J/g/d)", y = NULL)
  
  # --- environment: facetted, free y ---
  p_env <- base_theme(
      df %>% 
        filter(output_var %in% c("water_temp", "T_response", "rel_feeding"))
    ) +
    geom_line(aes(colour = output_var)) +
    facet_wrap(~output_var, scales = "free_y", ncol = 1) +
    labs(title = "Temperature & feeding", y = NULL)
  
  # --- food / ingestion: overlaid ---
  p_food <- base_theme(
      df %>% 
        filter(output_var %in% c("food_prov", "food_enc", "ing_pot", "ing_act", "weight")) %>% 
        pivot_wider(names_from = output_var, values_from = value) %>% 
        mutate(
          food_prov = food_prov/weight, 
          food_enc = food_enc/weight, 
          ing_pot = ing_pot/weight, 
          ing_act = ing_act/weight
        ) %>% 
        select(-weight) %>% 
        pivot_longer(cols = -days, names_to = "output_var", values_to = "value")
    ) +
    geom_line(aes(colour = output_var)) +
    labs(title = "Feeding & ingestion (g/g/d)", y = NULL)
  
  # --- O2 & NH4: overlaid ---
  p_o2_nh4 <- base_theme(filter(df, output_var %in% c("O2", "NH4"))) +
    geom_line(aes(colour = output_var)) +
    labs(title = "O2 & NH4 (g/g/d)", y = NULL)
  
  # --- excretion / uneaten: totals black, components coloured ---
  d_excr <- df %>% 
    filter(!str_detect(output_var, "carbon|nitrogen")) %>%
    filter(str_detect(output_var, "excr|uneat")) %>%
    mutate(
      type = str_split_i(output_var, "_", 1) %>% as.factor(),
      out = str_split_i(output_var, "_", 2) %>% as.factor()
    )
  
  p_excr <- base_theme(d_excr) +
    geom_line(aes(colour = type, linetype = out)) +
    scale_colour_manual(values = c("total" = "black", setNames(scales::hue_pal()(3), c("P", "L", "C")))) +
    labs(title = "Excretion & uneaten feed (g/g/d)", y = NULL)
  
  # --- elemental budgets: overlaid ---
  df_elem <- df %>% 
    filter(str_detect(output_var, "carbon|nitrogen")) %>%
    mutate(
      type = case_when(str_detect(output_var, "carbon") ~ "carbon", T ~ "nitrogen") %>% as.factor(),
      out = case_when(str_detect(output_var, "total") ~ "total", str_detect(output_var, "excr") ~ "excr", T ~ "uneat") %>% as.factor()
    )

  p_elem <- base_theme(df_elem) +
    geom_line(aes(colour = type, linetype = out)) +
    scale_linetype_manual(values = c("total" = "solid", "excr" = "dashed", "uneat" = "dotted")) +
    labs(title = "Carbon & nitrogen (g/g/d)", y = NULL)
  
  # --- assemble with patchwork ---
  (p_weight | p_dw | p_sgr) /
    (p_energy | p_env | p_food) /
    (p_o2_nh4 | p_excr | p_elem)
}

#' plot_sim_results(), trimmed for production-file troubleshooting
#'
#' Identical to `plot_sim_results()` except it drops every panel/variable
#' that the saved production files don't retain: SGR, E_assim, E_somat,
#' water_temp, food_enc, ing_pot, O2, and NH4.
plot_sim_results_troubleshooting <- function(df, CS = NA) {
  # --- helper: base theme ---
  base_theme <- function(d) {
    ggplot(d, aes(days, value)) +
      theme_classic(base_size = 10) +
      theme(
        legend.position  = "bottom", 
        legend.title     = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_blank()
      ) +
      labs(x = "Day")
  }
  
  # --- single-variable panels ---
  p_weight <- base_theme(filter(df, output_var == "weight")) +
    geom_line() + labs(title = "Individiaul weight (g)", y = NULL) 
  
  if (!is.na(CS)) {
    p_weight <- p_weight +
      geom_hline(yintercept = CS, linewidth = 0.75, colour = "red", linetype = "dashed") +
      annotate("text", x = min(df$days), y = CS, label = sprintf("CS: %d", CS), hjust = -0.1, vjust = -0.5, colour = "red")
  }
  
  p_dw <- base_theme(filter(df, output_var == "dw")) +
    geom_line() + labs(title = "Daily weight change (g/d)", y = NULL)
  
  # --- energy: overlaid lines (E_assim, E_somat not saved) ---
  p_energy <- base_theme(
    df %>% 
      filter(output_var %in% c("metab"))
    ) +
    geom_line(aes(colour = output_var)) +
    labs(title = "Energy (J/g/d)", y = NULL)
  
  # --- environment: facetted, free y (water_temp not saved) ---
  p_env <- base_theme(
      df %>% 
        filter(output_var %in% c("T_response", "rel_feeding"))
    ) +
    geom_line(aes(colour = output_var)) +
    facet_wrap(~output_var, scales = "free_y", ncol = 1) +
    labs(title = "Temperature & feeding", y = NULL)
  
  # --- food / ingestion: overlaid (food_enc, ing_pot not saved) ---
  p_food <- base_theme(
      df %>% 
        filter(output_var %in% c("food_prov", "ing_act", "weight")) %>% 
        pivot_wider(names_from = output_var, values_from = value) %>% 
        mutate(
          food_prov = food_prov/weight, 
          ing_act = ing_act/weight
        ) %>% 
        select(-weight) %>% 
        pivot_longer(cols = -days, names_to = "output_var", values_to = "value")
    ) +
    geom_line(aes(colour = output_var)) +
    labs(title = "Feeding & ingestion (g/g/d)", y = NULL)
  
  # --- excretion / uneaten: totals black, components coloured ---
  d_excr <- df %>% 
    filter(!str_detect(output_var, "carbon|nitrogen")) %>%
    filter(str_detect(output_var, "excr|uneat")) %>%
    mutate(
      type = str_split_i(output_var, "_", 1) %>% as.factor(),
      out = str_split_i(output_var, "_", 2) %>% as.factor()
    )
  
  p_excr <- base_theme(d_excr) +
    geom_line(aes(colour = type, linetype = out)) +
    scale_colour_manual(values = c("total" = "black", setNames(scales::hue_pal()(3), c("P", "L", "C")))) +
    labs(title = "Excretion & uneaten feed (g/g/d)", y = NULL)
  
  # --- elemental budgets: overlaid ---
  df_elem <- df %>% 
    filter(str_detect(output_var, "carbon|nitrogen")) %>%
    mutate(
      type = case_when(str_detect(output_var, "carbon") ~ "carbon", T ~ "nitrogen") %>% as.factor(),
      out = case_when(str_detect(output_var, "total") ~ "total", str_detect(output_var, "excr") ~ "excr", T ~ "uneat") %>% as.factor()
    )

  p_elem <- base_theme(df_elem) +
    geom_line(aes(colour = type, linetype = out)) +
    scale_linetype_manual(values = c("total" = "solid", "excr" = "dashed", "uneat" = "dotted")) +
    labs(title = "Carbon & nitrogen (g/g/d)", y = NULL)
  
  # --- assemble with patchwork (SGR and O2/NH4 panels dropped) ---
  (p_weight | p_dw) /
    (p_energy | p_env | p_food) /
    (p_excr | p_elem)
}


