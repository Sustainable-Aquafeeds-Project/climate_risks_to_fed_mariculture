# Make sure you've run "country-production.R" first

library(tidyverse)
library(here)
library(qs2)
library(rnaturalearth)

source("src/dirs.R")
source("src/model_functions.R")
source("src/other_functions.R")

# Files and regional information -----------------------------------------------------------------------------------------

prod_files <- here("for-europe-seminars", "data") %>% 
  list.files(full.names = T, pattern = ".qs")
isos <- prod_files %>% basename() %>% str_remove(".qs")
prod <- data.frame(iso = isos, fnm = prod_files)

regions <- ne_countries(scale = "medium", returnclass = "sf") |>
  sf::st_drop_geometry() |>
  select(name, iso_a3, iso_a3_eh, adm0_a3, continent, region_un, subregion, region_wb) %>% 
  filter(adm0_a3 %in% isos)
regions_PSE <- regions %>% # Need to add Palestine
  filter(name == "Israel") %>% 
  mutate(
    name = "Palestine",
    iso_a3 = "PSE",
    iso_a3_eh = "PSE",
    adm0_a3 = "PSE"
  )
regions <- rbind(regions, regions_PSE)


GLW_ls <- list(
  "GLW_1.5" = 2018:2037,
  "GLW_2.0" = 2032:2051,
  "GLW_3.0" = 2055:2074,
  "GLW_4.0" = 2075:2094
)
all_years <- c(2018:2037, 2032:2051, 2055:2074, 2075:2094) %>% unique()

# Get change in biomass per region ---------------------------------------------------------------------------------------
# reg <- unique(regions$region_wb)[1]
biomass_data <- map_dfr(
  unique(regions$region_wb),
  function(reg) {
    isos <- regions$adm0_a3[regions$region_wb == reg]
    fnms <- prod %>% filter(iso %in% isos) %>% pull(fnm)

    df <- map_dfr(
      fnms,
      function(fnm) {
        qs_read(fnm) %>% 
          filter(year %in% all_years) %>% 
          mutate(biomass_change = biomass_change/biomass_2025) %>% # biomass change is already (biomass - biomass_2025)
          select(year, biomass_2025, biomass_change)
      }
    ) %>% 
      mutate(region = reg)
  },
  .progress = T
)

regional_biomass_data <- biomass_data %>% 
  group_by(region, year) %>% 
  reframe(mean_biomass_change = weighted.mean(biomass_change, w = biomass_2025))

# yrs <- GLW_ls[[1]]
# glw <- names(GLW_ls)[1]
regional_biomass_data <- map2_dfr(
  GLW_ls,
  names(GLW_ls),
  function(yrs, glw) {
    regional_biomass_data %>% 
      filter(year %in% yrs) %>% 
      group_by(region) %>% 
      reframe(
        median = median(mean_biomass_change),
        q25 = quantile(mean_biomass_change, 0.25),
        q75 = quantile(mean_biomass_change, 0.75)
      ) %>% 
      mutate(GLW = as.character(str_remove(glw, "GLW_")))
  },
  .progress = T
) %>% 
  mutate(GLW = factor(GLW, levels = c("1.5", "2.0", "3.0", "4.0"), labels = c("1.5°C", "2°C", "3°C", "4°C")))

# Put it into a plot like Gage's -----------------------------------------------------------------------------------------

warming_cols <- c("1.5°C" = "#E9A23B", "2°C"   = "#DB5F32", "3°C"   = "#B3243B", "4°C"   = "#5C0F2E")
offsets <- c("1.5°C" = -0.24, "2°C" = -0.08, "3°C" = 0.08, "4°C" = 0.24)
region_order <- regional_biomass_data |>
  filter(GLW == "4°C") |>
  arrange(median) |>   # bottom-to-top
  pull(region)

p <- regional_biomass_data %>% 
  mutate(
    region = factor(region, levels = region_order),
    y = as.numeric(region) + offsets[as.character(GLW)]
  ) %>% 
  ggplot(aes(x = y, y = median, colour = GLW)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.6) +
  geom_segment(aes(y = q25, yend = q75), linewidth = 0.8) +
  geom_point(aes(y = median), size = 2.3) +
  scale_colour_manual(values = warming_cols, name = NULL) +
  scale_y_continuous(
    breaks = seq(-0.25, 0.75, 0.25),
    labels = scales::label_percent(),
    limits = c(-0.35, 0.35),
    expand = expansion(0)
  ) +
  scale_x_continuous(
    breaks = seq_along(region_order),
    labels = region_order,
    expand = expansion(add = 0.5)
  ) +
  coord_flip() +
  labs(
    title    = "Projected changes in farm production",
    subtitle = "Median and interquartile range across countries within each region",
    y        = "Projected change (%) weighted by current (2025) production",
    x        = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 0, size = 3))) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title            = element_text(face = "bold", size = 14),
    plot.subtitle         = element_text(size = 10, margin = margin(b = 6)),
    plot.title.position   = "plot",
    plot.caption.position = "plot",
    plot.caption          = element_text(hjust = 0, colour = "grey45", size = 9),
    legend.position       = "top",
    legend.justification  = "center",
    legend.text           = element_text(size = 10),
    legend.margin         = margin(b = -4),
    panel.grid.major.y    = element_blank(),
    panel.grid.minor      = element_blank(),
    panel.grid.major.x    = element_line(colour = "grey92", linewidth = 0.4),
    axis.text.y           = element_text(colour = "grey25", size = 10),
    axis.title.x          = element_text(margin = margin(t = 6)),
    plot.margin           = margin(10, 15, 10, 10)
  )

  ggsave(
    plot = p,
    filename = here("for-europe-seminars", "plots", "regional_biomass_change.png"),
    width = 10,
    height = 5.5
  )

# Get change in food required per region ---------------------------------------------------------------------------------
# reg <- unique(regions$region_wb)[1]
feed_data <- map_dfr(
  unique(regions$region_wb),
  function(reg) {
    isos <- regions$adm0_a3[regions$region_wb == reg]
    fnms <- prod %>% filter(iso %in% isos) %>% pull(fnm)

    df <- map_dfr(
      fnms,
      function(fnm) {
        qs_read(fnm) %>% 
          filter(year %in% all_years) %>% 
          mutate(food_prov_biom_change = food_prov_biom_change/food_prov_biom_2025) %>% 
          select(year, biomass_2025, food_prov_biom_change)
      }
    ) %>% 
      mutate(region = reg)
  },
  .progress = T
)

regional_feed_data <- feed_data %>% 
  group_by(region, year) %>% 
  reframe(mean_food_prov_biom_change = weighted.mean(food_prov_biom_change, w = biomass_2025))

# yrs <- GLW_ls[[1]]
# glw <- names(GLW_ls)[1]
regional_feed_data <- map2_dfr(
  GLW_ls,
  names(GLW_ls),
  function(yrs, glw) {
    regional_feed_data %>% 
      filter(year %in% yrs) %>% 
      group_by(region) %>% 
      reframe(
        median = median(mean_food_prov_biom_change),
        q25 = quantile(mean_food_prov_biom_change, 0.25),
        q75 = quantile(mean_food_prov_biom_change, 0.75)
      ) %>% 
      mutate(GLW = as.character(str_remove(glw, "GLW_")))
  },
  .progress = T
) %>% 
  mutate(GLW = factor(GLW, levels = c("1.5", "2.0", "3.0", "4.0"), labels = c("1.5°C", "2°C", "3°C", "4°C")))

# Put it into a plot like Gage's -----------------------------------------------------------------------------------------

warming_cols <- c("1.5°C" = "#E9A23B", "2°C"   = "#DB5F32", "3°C"   = "#B3243B", "4°C"   = "#5C0F2E")
offsets <- c("1.5°C" = -0.24, "2°C" = -0.08, "3°C" = 0.08, "4°C" = 0.24)

p <- regional_feed_data %>% 
  mutate(
    region = factor(region, levels = region_order),
    y = as.numeric(region) + offsets[as.character(GLW)]
  ) %>% 
  ggplot(aes(x = y, y = median, colour = GLW)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.6) +
  geom_segment(aes(y = q25, yend = q75), linewidth = 0.8) +
  geom_point(aes(y = median), size = 2.3) +
  scale_colour_manual(values = warming_cols, name = NULL) +
  scale_y_continuous(
    breaks = seq(-0.25, 0.75, 0.05),
    labels = scales::label_percent(),
    limits = c(-0.05, 0.16),
    expand = expansion(0)
  ) +
  scale_x_continuous(
    breaks = seq_along(region_order),
    labels = region_order,
    expand = expansion(add = 0.5)
  ) +
  coord_flip() +
  labs(
    title    = "Projected changes in the feed required per unit biomass",
    subtitle = "Median and interquartile range across countries within each region",
    y        = "Projected change (%) weighted by current (2025) production",
    x        = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 0, size = 3))) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title            = element_text(face = "bold", size = 14),
    plot.subtitle         = element_text(size = 10, margin = margin(b = 6)),
    plot.title.position   = "plot",
    plot.caption.position = "plot",
    plot.caption          = element_text(hjust = 0, colour = "grey45", size = 9),
    legend.position       = "top",
    legend.justification  = "center",
    legend.text           = element_text(size = 10),
    legend.margin         = margin(b = -4),
    panel.grid.major.y    = element_blank(),
    panel.grid.minor      = element_blank(),
    panel.grid.major.x    = element_line(colour = "grey92", linewidth = 0.4),
    axis.text.y           = element_text(colour = "grey25", size = 10),
    axis.title.x          = element_text(margin = margin(t = 6)),
    plot.margin           = margin(10, 15, 10, 10)
  )

  ggsave(
    plot = p,
    filename = here("for-europe-seminars", "plots", "regional_feedrequired_change.png"),
    width = 10,
    height = 5.5
  )

# Country-level change maps, facetted by GLW -----------------------------------------------------------------------------

# Per-country, per-year change (weighted by 2025 biomass), then mean across the years in each GLW window
country_data <- map_dfr(
  prod$fnm,
  function(fnm) {
    qs_read(fnm) %>% 
      filter(year %in% all_years) %>% 
      group_by(iso = ISO3_Code, year) %>% 
      reframe(
        biomass_change = weighted.mean(biomass_change/biomass_2025, w = biomass_2025),
        food_prov_biom_change = weighted.mean(food_prov_biom_change/food_prov_biom_2025, w = biomass_2025)
      )
  },
  .progress = T
)

country_GLW_data <- map2_dfr(
  GLW_ls,
  names(GLW_ls),
  function(yrs, glw) {
    country_data %>% 
      filter(year %in% yrs) %>% 
      group_by(iso) %>% 
      reframe(
        biomass_change = mean(biomass_change),
        food_prov_biom_change = mean(food_prov_biom_change)
      ) %>% 
      mutate(GLW = as.character(str_remove(glw, "GLW_")))
  }
) %>% 
  mutate(GLW = factor(GLW, levels = c("1.5", "2.0", "3.0", "4.0"), labels = c("1.5°C", "2°C", "3°C", "4°C")))

world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  # Natural Earth codes Palestine as PSX; match the ISO code used in the production files
  mutate(iso = if_else(iso_a3 == "PSE", "PSE", adm0_a3)) %>% 
  filter(name != "Antarctica") %>% 
  select(iso, geometry)

country_GLW_sf <- world %>% 
  inner_join(country_GLW_data, by = "iso")

# Limits default to the full data range (asymmetric), with white fixed at 0
plot_GLW_map <- function(
  var, title, 
  limits = range(country_GLW_sf[[var]]), 
  midpoint = 0
) {
  ggplot() +
    geom_sf(data = world, fill = "grey80", colour = "grey45", linewidth = 0.1) +
    geom_sf(data = country_GLW_sf, aes(fill = .data[[var]]), colour = "grey30", linewidth = 0.1) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC", midpoint = midpoint,
      limits = limits,
      # Values outside the limits are shown at the extreme colour rather than dropped
      oob = scales::squish,
      labels = scales::label_percent(),
      name = "Change (%)"
    ) +
    coord_sf(crs = "+proj=robin", expand = FALSE) +
    facet_wrap(~GLW, ncol = 2) +
    labs(
      title = title,
      subtitle = "Mean of each GLW across species, weighted by 2025 production"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title            = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle         = element_text(hjust = 0.5),
      strip.text            = element_text(face = "bold", size = 11),
      legend.position       = "bottom",
      legend.title.position = "bottom",
      legend.title          = element_text(hjust = 0.5),
      legend.frame          = element_rect(colour = "grey50", linewidth = 0.4),
      legend.key.width      = unit(1.5, "cm"),
      panel.grid            = element_blank(),
      axis.text             = element_blank()
    )
}

p <- plot_GLW_map(
  var = "biomass_change",
  title = "Projected changes in farm production by country",
  limits = unname(quantile(country_GLW_sf$biomass_change, c(0.05, 0.95)))
)
ggsave(
  plot = p,
  filename = here("for-europe-seminars", "plots", "country_biomass_change_map.png"),
  width = 12,
  height = 7
)

p <- plot_GLW_map(
  var = "food_prov_biom_change",
  title = "Projected changes in the feed required per unit biomass by country",
  # A few countries (e.g. BHR, TUR) are far above the rest, so clip to the 5th-95th percentiles
  limits = unname(quantile(country_GLW_sf$food_prov_biom_change, c(0.05, 0.95)))
)

ggsave(
  plot = p,
  filename = here("for-europe-seminars", "plots", "country_feedrequired_change_map.png"),
  width = 12,
  height = 7
)
