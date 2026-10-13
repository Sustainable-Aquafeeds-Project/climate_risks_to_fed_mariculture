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
