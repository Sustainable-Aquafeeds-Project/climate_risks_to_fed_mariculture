# Make sure you've run "3_regional-and-maps-of-HS-change.R" first - it defines get_HS_change_rast()

library(tidyverse)
library(here)
library(qs2)
library(rnaturalearth)
library(terra)

source("src/dirs.R")
source("src/model_functions.R")
source("src/other_functions.R")

stopifnot("Run 3_regional-and-maps-of-HS-change.R first" = exists("get_HS_change_rast"))

# Files and regional information -----------------------------------------------------------------------------------------
land <- ne_countries(scale = "large", returnclass = "sf")

# Coastal gap-filling ----------------------------------------------------------------------------------------------------
# Band around the coast of a single country, within which NA cells are allowed to be filled
get_coast_zone <- function(iso, width_m = 62500) {
  this_land <- land %>% filter(iso_a3 == iso | adm0_a3 == iso)
  if (nrow(this_land) == 0) {
    warning("No land polygon found for ", iso, "; skipping coastal gap-filling")
    return(NULL)
  }
  buffer(vect(this_land), width_m)
}

fill_coastal_gaps <- function(r, coast_zone, w = 3) {
  if (is.null(coast_zone)) return(r)
  filled <- focal(r, w, "mean", na.rm = TRUE, na.policy = "only")
  cover(r, mask(filled, coast_zone))       # only fill NAs that fall inside the coastal band
}

dir.create(here("for-europe-seminars", "plots", "HS_change"), showWarnings = FALSE)
dir.create(here("for-europe-seminars", "plots", "HS_change_map"), showWarnings = FALSE)

hs_files <- file.path(outs_path, "data", "habitat_suitability_decadal") %>% 
  list.files(recursive = T, full.names = T, pattern = "GLW|ref") %>% 
  str_subset("MR") #%>%                                     # Use only the metabolism HS
  # str_subset(paste0(c("ESP", "MMR"), collapse = "|"))      # Picking out some example countries

isos <- hs_files %>% dirname() %>% dirname() %>% basename()
spec <- hs_files %>% dirname() %>% basename()

hs <- data.frame(
  iso = isos, 
  spec = spec, 
  fnm = hs_files
)

assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs_read() %>% 
  sf::st_drop_geometry() %>% 
  # filter(ISO3_Code %in% c("ESP", "MMR")) %>% 
  group_by(ISO3_Code, country, model_name) %>% 
  reframe(value = sum(value))

assigned_farms_points <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs_read() # %>% filter(ISO3_Code %in% c("ESP", "MMR"))

# Get change in HS per country (EEZ) --------------------------------------------------------------------------------------
# i <- unique(isos)[1]
HS_data <- map_dfr(
  unique(isos),
  function(i) {
    # Filter to only the species that are actually grown there
    af_iso <- assigned_farms %>% filter(ISO3_Code == i)
    hsi_spec <- hs %>% filter(iso == i & spec %in% af_iso$model_name)
    hsi_spec <- split(hsi_spec, hsi_spec$spec)
    coast_zone <- get_coast_zone(i)

    # hsis <- hsi_spec[[1]]
    # nm <- names(hsi_spec)[1]
    map2_dfr(
      hsi_spec,
      names(hsi_spec),
      function(hsis, nm) {
        r <- get_HS_change_rast(hsis$fnm) #%>% fill_coastal_gaps(coast_zone)

        # Get mean change in habitat suitability - whole EEZ
        df1 <- rbind(
          global(r[["GLW_1.5"]], "mean", na.rm = T),
          global(r[["GLW_2.0"]], "mean", na.rm = T),
          global(r[["GLW_3.0"]], "mean", na.rm = T),
          global(r[["GLW_4.0"]], "mean", na.rm = T)
        ) %>% rownames_to_column()
        df2 <- rbind(
          global(r[["GLW_1.5"]], "sd", na.rm = T),
          global(r[["GLW_2.0"]], "sd", na.rm = T),
          global(r[["GLW_3.0"]], "sd", na.rm = T),
          global(r[["GLW_4.0"]], "sd", na.rm = T)
        ) %>% rownames_to_column()

        left_join(df1, df2, by = "rowname") %>% 
          mutate(spec = nm)
      }
    ) %>% 
    rename(GLW = rowname, mean_HS_change = mean, sd_HS_change = sd) %>% 
    left_join(
      af_iso,
      by = c("spec" = "model_name")
    )
  },
  .progress = T
)

warming_cols <- c("1.5°C" = "#E9A23B", "2°C"   = "#DB5F32", "3°C"   = "#B3243B", "4°C"   = "#5C0F2E")
warming_levs <- c("1.5°C" = "GLW_1.5", "2°C" = "GLW_2.0", "3°C" = "GLW_3.0", "4°C" = "GLW_4.0")
lims1 <- c(-0.55, 0.4)
breaks1 <- seq(-0.5, 1, 0.1)

walk(
  unique(isos),
  function(i) {
    p1 <- HS_data %>% 
      filter(ISO3_Code == i) %>% 
      mutate(GLW = factor(GLW, levels = warming_levs, labels = names(warming_levs))) %>% 
      ggplot(aes(
        x = spec, y = mean_HS_change, 
        ymin = mean_HS_change - sd_HS_change,
        ymax = mean_HS_change + sd_HS_change,
        colour = GLW
      )) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
      geom_errorbar(position = position_dodge(width = 0.5), width = 0.25) +
      scale_y_continuous(
        # limits = lims1, breaks = breaks1,
        labels = scales::label_percent()
      ) +
      scale_colour_manual(values = warming_cols) +
      prettyplot() +
      theme(legend.position = "top") +
      labs(
        x = NULL, y = "Change in habitat suitability",
        title = paste0("Projected change in habitat suitability (", i, ")"),
        subtitle = "Mean and SD of all cells within the country EEZ"
      )

    ggsave(
      plot = p1,
      filename = here("for-europe-seminars", "plots", "HS_change", paste0(i, ".png")),
      width = 6, height = 5
    )
  }
)

# Get change in HS per country (EEZ) as a map -----------------------------------------------------------------------------
lims <- c(-0.5, 0.5)

# i <- unique(isos)[1]
walk(
  unique(isos),
  function(i) {
    # Filter to only the species that are actually grown there
    af_iso <- assigned_farms %>% filter(ISO3_Code == i)
    af_points <- assigned_farms_points %>% filter(ISO3_Code == i)
    hsi_spec <- hs %>% filter(iso == i & spec %in% af_iso$model_name)
    hsi_spec <- split(hsi_spec, hsi_spec$spec)
    coast_zone <- get_coast_zone(i)

    # hsis <- hsi_spec[[1]]
    hsi_rasts <- map(hsi_spec, \(hsis) get_HS_change_rast(hsis$fnm) %>% 
      fill_coastal_gaps(coast_zone)) %>% 
      setNames(names(hsi_spec))

    w <- setNames(af_iso$value, af_iso$model_name)
    w <- w[names(hsi_rasts)]   # drops drum_croaker and tuna, and matches the list order
    glw <- names(hsi_rasts[[1]])

    hsi_wm <- lapply(glw, \(g) {
      s <- rast(lapply(hsi_rasts, \(r) r[[g]]))   # one layer per species
      weighted.mean(s, w, na.rm = TRUE)
      }) |>
        rast() |>
        setNames(glw)
    
    country_bbox <- bbox(i, crs = 4326)
    
    p1 <- ggplot() +
      geom_spatraster(data = hsi_wm) +
      geom_sf(data = land, fill = "grey85", colour = "grey55", linewidth = 0.2) + 
      geom_sf(data = af_points, shape = 4, size = 1.5, colour = "black") + 
      coord_sf(
        xlim = c(country_bbox[["xmin"]], country_bbox[["xmax"]]), 
        ylim = c(country_bbox[["ymin"]], country_bbox[["ymax"]])
      ) +
      scale_fill_viridis_c(
        option = "viridis",
        # limits = lims, 
        labels = scales::label_percent(),
        na.value = "transparent"
      ) +
      facet_wrap(~lyr) +
      labs(
        title = paste0("Change in habitat suitability (", i, ")"), 
        subtitle = "Weighted mean for all species",
        fill = NULL, x = "Longitude", y = "Latitude"
      ) +
      guides(fill = guide_colourbar(
        theme = theme(
          legend.key.height     = unit(14, "lines"),
          legend.key.width      = unit(1, "lines"),
          legend.ticks          = element_line(colour = "white", linewidth = 0.4),
          legend.ticks.length   = unit(0.2, "lines"),
          legend.frame          = element_rect(colour = "grey30", linewidth = 0.3),
          legend.title.position = "top"
        )
      )) +
      theme_classic() +
      theme(
        legend.position = "right",
        legend.title    = element_text(hjust = 0.5)
      )
    
    ggsave(
      plot = p1,
      filename = here("for-europe-seminars", "plots", "HS_change_map", paste0(i, ".png")),
      width = 9.5, height = 6.5
    )    
  },
  .progress = T
)


