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

# Feeding response with temperature raster only
get_feeding_responses_r <- function(
  temp_stack, species_params,
  out_cols = c("feeding_response")
) {
  ncell <- terra::ncell(temp_stack)
  nlyr  <- terra::nlyr(temp_stack)
  temp_vec <- as.vector(terra::values(temp_stack))

  # Only run the model on real temperatures; cells masked outside the EEZ/land are NA
  ok   <- is.finite(temp_vec)
  resp <- data.frame(feeding_response = feeding_rate(temp_vec[ok], species_params))

  result <- lapply("feeding_response", function(col) {
    full <- rep(NA_real_, length(temp_vec))      # scatter results back to full length
    full[ok] <- resp[[col]]
    r <- terra::setValues(temp_stack,            # reuse the stack purely as a geometry template
                          matrix(full, nrow = ncell, ncol = nlyr))
    names(r) <- names(temp_stack)                # keep the layer names
    r
  })
  names(result) <- "feeding_response"
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

country_summary <- function(iso, digits = 1) {
  d <- file.path(prepdata_path, "FAO_with_models.qs") %>% 
    qs_read() %>% 
    filter(period %in% 2014:2023 & ISO3_Code == iso) %>%
    mutate(model_name = case_when(
      model_name %in% c("drum", "large_yellow_croaker", "whitemouth_croaker", "meagre") ~ "drum_croaker", 
      model_name %in% c("red_seabream", "gilthead_seabream", "blackhead_seabream") ~ "seabream", 
      TRUE ~ model_name)
    ) %>% 
    group_by(ISO3_Code, country, model_name, Scientific_Name, period) %>% 
    reframe(value = sum(value)) %>% 
    group_by(ISO3_Code, country, model_name, Scientific_Name) %>% 
    reframe(value = mean(value))

  d <- d |>
    mutate(value = units::drop_units(value))  # drop [t] units

  if (nrow(d) == 0) stop("No rows found for ISO3_Code '", iso, "'")

  total <- sum(d$value, na.rm = TRUE)

  # % of country total per model
  models <- d |>
    summarise(value = sum(value, na.rm = TRUE), .by = model_name) |>
    mutate(pct = 100 * value / total) |>
    arrange(desc(value))

  # % within each model per species
  species <- d |>
    summarise(value = sum(value, na.rm = TRUE), .by = c(model_name, Scientific_Name)) |>
    mutate(pct = 100 * value / sum(value), .by = model_name) |>
    arrange(model_name, desc(value))

  fmt_t <- \(x) format(round(x), big.mark = ",", trim = TRUE)
  fmt_p <- \(x) sprintf("%.*f%%", digits, x)

  cat(sprintf("%s (%s)\n", d$country[1], iso))
  cat(sprintf("Total production: %s t\n", fmt_t(total)))
  cat(sprintf("Number of models: %d\n", nrow(models)))

  for (i in seq_len(nrow(models))) {
    m <- models$model_name[i]
    cat(sprintf("\n  %s: %s t (%s of country total)\n",
                m, fmt_t(models$value[i]), fmt_p(models$pct[i])))

    sp <- filter(species, model_name == m)
    cat(sprintf("    - %s: %s t (%s of model)\n",
                sp$Scientific_Name, fmt_t(sp$value), fmt_p(sp$pct)), sep = "")
  }

  invisible(list(total = total, models = models, species = species))
}

# country_summary("DZA")

get_country_data <- function(fnms, ids) {
  annual_biomass <- map2_dfr(
    fnms, ids,
    function(fnm, id) {
      qs_read(fnm) %>% 
        tidy_stat("weight_scaled") %>% 
        mutate(farm_ID = id) %>% 
        slice_max(prod_day, n = 1) %>% 
        select(farm_ID, year, mean)
    }
  ) %>% rename(biomass = mean)

  annual_food_prov <- map2_dfr(
    fnms, ids,
    function(fnm, id) {
      qs_read(fnm) %>% 
        tidy_stat("food_prov_scaled") %>% 
        mutate(farm_ID = id) %>% 
        group_by(farm_ID, year) %>% 
        reframe(food_prov = sumna(mean))
    }
  ) %>% 
    left_join(annual_biomass, by = join_by(farm_ID, year)) %>% 
    mutate(food_prov_biom = food_prov/biomass)

  annual_food_prov_2025 <- annual_food_prov %>% 
    filter(year == 2025) %>% 
    select(-c(year)) %>% 
    rename(
      biomass_2025        = biomass,
      food_prov_2025      = food_prov, 
      food_prov_biom_2025 = food_prov_biom
    )

  left_join(
    annual_food_prov,
    annual_food_prov_2025,
    by = join_by(farm_ID)
  ) %>% 
  mutate(
    food_prov_change           = food_prov - food_prov_2025,
    food_prov_change_pct       = food_prov_change/food_prov_2025,
    food_prov_biom_change      = food_prov_biom - food_prov_biom_2025,
    food_prov_biom_change_pct  = food_prov_biom_change/food_prov_biom_2025,
    biomass_change             = biomass - biomass_2025,
    biomass_change_pct         = biomass_change/biomass_2025
  )
}

get_country_data2 <- function(fnms) {  
  model_name  <- str_extract(fnms, "(?<=_).*(?=\\.qs$)")
  cons_data   <- map2_dfr(fnms, function(fn, mn) {
      qs_read(fnm) %>% mutate(model_name = mn)
  })

  cons_data_2025 <- cons_data %>% 
    filter(year == 2025) %>% 
    select(-c(year)) %>% 
    rename(
      biomass_2025        = biomass,
      feed_prov_2025      = feed_prov, 
      feed_prov_per_biomass_2025 = feed_prov_per_biomass
    )

  left_join(
    cons_data,
    cons_data_2025,
    by = join_by(farm_ID, model_name)
  ) %>% 
  mutate(
    feed_prov_change           = feed_prov - feed_prov_2025,
    feed_prov_change_pct       = feed_prov_change/feed_prov_2025,
    feed_prov_biom_change      = feed_prov_biom - feed_prov_per_biomass_2025,
    feed_prov_biom_change_pct  = feed_prov_biom_change/feed_prov_per_biomass_2025,
    biomass_change             = biomass - biomass_2025,
    biomass_change_pct         = biomass_change/biomass_2025
  )
}

country_coord_sf <- function(bbox_data) {
  coord_sf(
    xlim = c(bbox_data[["xmin"]], bbox_data[["xmax"]]), 
    ylim = c(bbox_data[["ymin"]], bbox_data[["ymax"]])
  ) 
}



# Percentage y-axis helpers ---------------------------------------------------
#
# Limits are rounded outward to a "round" percentage just outside the data
# range, and breaks run between them in 10% steps, or 5% steps if 10% would
# give too few breaks. Values are assumed to be proportions (0.85 = 85%).
#
# pct_limits() and pct_breaks() can be passed straight to
# scale_y_continuous(limits = , breaks = ), or use scale_y_pct() for both.

# Internal: choose step and compute rounded limits from a numeric vector
.pct_axis <- function(x,
                      steps = c(0.10, 0.05),
                      min_breaks = 5,
                      strict = FALSE,
                      clamp = c(0, 1)) {
  rng <- range(x, na.rm = TRUE, finite = TRUE)
  if (!all(is.finite(rng))) stop("`x` has no finite values.")
  steps <- sort(steps, decreasing = TRUE)

  round_out <- function(step) {
    # round() guards against floating point, e.g. 0.85 / 0.05 = 16.9999...
    lo <- floor(round(rng[1] / step, 8)) * step
    hi <- ceiling(round(rng[2] / step, 8)) * step
    if (strict) {
      if (isTRUE(all.equal(lo, rng[1]))) lo <- lo - step
      if (isTRUE(all.equal(hi, rng[2]))) hi <- hi + step
    }
    if (hi == lo) { lo <- lo - step; hi <- hi + step }  # constant data
    # Don't let rounding push past the clamp (e.g. 0% / 100%), unless the
    # data are already outside it on that side
    if (!is.null(clamp)) {
      if (rng[1] >= clamp[1]) lo <- max(lo, clamp[1])
      if (rng[2] <= clamp[2]) hi <- min(hi, clamp[2])
    }
    c(lo, hi)
  }

  # Use the largest step that still gives at least `min_breaks` breaks;
  # fall back to the smallest step
  for (step in steps) {
    lim <- round_out(step)
    n <- round((lim[2] - lim[1]) / step) + 1
    if (n >= min_breaks) break
  }
  list(limits = round(lim, 10), step = step)
}

# Rounded limits. Works as `limits =` in scale_y_continuous() (it receives
# the data range) or on a raw data vector.
pct_limits <- function(x, ...) {
  .pct_axis(x, ...)$limits
}

# Breaks between the rounded limits. Works as `breaks =` in
# scale_y_continuous() (it receives the limits) or on a raw data vector.
pct_breaks <- function(x, ...) {
  ax <- .pct_axis(x, ...)
  lim <- ax$limits
  # Breaks on multiples of the step, in case clamping left an odd limit
  first <- ceiling(round(lim[1] / ax$step, 8)) * ax$step
  round(seq(first, lim[2] + 1e-9, by = ax$step), 10)
}

# Convenience wrapper: percentage scale with rounded limits and breaks.
# `...` goes to scale_y_continuous(); pass axis options via `axis_args`,
# e.g. axis_args = list(steps = c(0.1, 0.05), strict = TRUE).
scale_y_pct <- function(...,
                        axis_args = list(),
                        accuracy = 1,
                        expand = c(0, 0)) {
  # expand = c(0, 0) so the axis ends exactly on the rounded limits. ggplot
  # passes the *expanded* range to the breaks function, so any expansion can
  # make pct_breaks() pick a different step than pct_limits() did.
  scale_y_continuous(
    limits = function(x) do.call(pct_limits, c(list(x), axis_args)),
    breaks = function(x) do.call(pct_breaks, c(list(x), axis_args)),
    labels = scales::label_percent(accuracy = accuracy),
    expand = expand,
    ...
  )
}
