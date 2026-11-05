library(tidyverse)
library(here)
library(qs2)
library(rnaturalearth)
library(terra)

source("src/dirs.R")
source("src/model_functions.R")
source("src/other_functions.R")

# Files and regional information -----------------------------------------------------------------------------------------
hs_files <- file.path(outs_path, "data", "habitat_suitability_decadal") %>% 
  list.files(recursive = T, full.names = T, pattern = "GLW|ref") %>% 
  str_subset("MR")
isos <- hs_files %>% dirname() %>% dirname() %>% basename()
spec <- hs_files %>% dirname() %>% basename()
hs <- data.frame(
  iso = isos, 
  spec = spec, 
  fnm = hs_files
) #%>% filter(!spec %in% c("golden_pompano", "mud_crab", "nile_tilapia", "drum_croaker"))

assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs_read() %>% 
  sf::st_drop_geometry() %>% 
  group_by(ISO3_Code, country, model_name) %>% 
  reframe(value = sum(value))

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

# HS change function (also used in 4_per-country-HS-change.R) ------------------------------------------------------------
# Takes the file paths for one country x species (ref + GLW rasters) and returns a 4-layer SpatRaster of
# proportional change in HS relative to the reference period. Positive = HS increases under warming.
GLW_names <- c("GLW_1.5", "GLW_2.0", "GLW_3.0", "GLW_4.0")

get_HS_change_rast <- function(fnms) {
  r <- rast(fnms)/365
  names(r) <- str_remove_all(basename(fnms), "MR_|.tif")
  ref <- r[["ref_2015_2024"]]
  r <- (r[[GLW_names]] - ref)/ref
  names(r) <- GLW_names
  r
}

# Get change in HS per region (EEZ) --------------------------------------------------------------------------------------
# reg <- unique(regions$region_wb)[1]
HS_data <- map_dfr(
  unique(regions$region_wb),
  function(reg) {
    isos <- regions$adm0_a3[regions$region_wb == reg]
    hs_iso <- hs %>% filter(iso %in% isos)
    hs_iso <- split(hs_iso, hs_iso$iso)

    # hsi <- hs_iso[[1]]
    df_iso <- map_dfr(
      hs_iso,
      function(hsi) {
        # Filter to only the species that are actually grown there
        af_iso <- assigned_farms %>% filter(ISO3_Code == unique(hsi$iso))
        hsi_spec <- hsi %>% filter(spec %in% af_iso$model_name)
        hsi_spec <- split(hsi_spec, hsi_spec$spec)

        if (length(hsi_spec) != 0) {
        # hsis <- hsi_spec[[1]]
        map2_dfr(
          hsi_spec,
          names(hsi_spec),
          function(hsis, nm) {
            r <- get_HS_change_rast(hsis$fnm)

            # Get mean change in habitat suitability - whole EEZ
            rbind(
              global(r[["GLW_1.5"]], "mean", na.rm = T),
              global(r[["GLW_2.0"]], "mean", na.rm = T),
              global(r[["GLW_3.0"]], "mean", na.rm = T),
              global(r[["GLW_4.0"]], "mean", na.rm = T)
            ) %>% rownames_to_column() %>% 
              mutate(spec = nm)
          }
        ) %>% 
          rename(GLW = rowname, mean_HS_change = mean) %>% 
          left_join(
            af_iso,
            by = c("spec" = "model_name")
          )
        }
      }
    ) %>% 
      mutate(region = reg)

    # Get average per country, weighted by production
    df_iso <- df_iso %>% 
      group_by(region, ISO3_Code, country, GLW) %>% 
      reframe(mean_HS_change = weighted.mean(mean_HS_change, value))
  },
  .progress = T
)

# Get stats per region
regional_HS_data <- HS_data %>% 
  group_by(region, GLW) %>% 
  reframe(
    median = median(mean_HS_change, na.rm = T),
    q25 = quantile(mean_HS_change, 0.25, na.rm = T),
    q75 = quantile(mean_HS_change, 0.75, na.rm = T)
  ) %>% 
  mutate(GLW = factor(
    GLW, 
    levels = c("GLW_1.5", "GLW_2.0", "GLW_3.0", "GLW_4.0"), 
    labels = c("1.5°C", "2°C", "3°C", "4°C")
  ))

# Put it into a plot like Gage's -----------------------------------------------------------------------------------------

warming_cols <- c("1.5°C" = "#E9A23B", "2°C"   = "#DB5F32", "3°C"   = "#B3243B", "4°C"   = "#5C0F2E")
offsets <- c("1.5°C" = -0.24, "2°C" = -0.08, "3°C" = 0.08, "4°C" = 0.24)
region_order <- c("Sub-Saharan Africa", "South Asia", "East Asia & Pacific", "Middle East & North Africa", "North America", "Latin America & Caribbean", "Europe & Central Asia")

p <- regional_HS_data %>% 
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
    breaks = seq(-1, 0.75, 0.1),
    labels = scales::label_percent(),
    # limits = c(-0.25, 0.41),
    expand = expansion(0)
  ) +
  scale_x_continuous(
    breaks = seq_along(region_order),
    labels = region_order,
    expand = expansion(add = 0.5)
  ) +
  coord_flip() +
  labs(
    title    = "Projected changes in habitat suitability within EEZs",
    subtitle = "Median and interquartile range across countries within each region",
    y        = "Projected change (%) weighted by current (2025) production of fed species",
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
    filename = here("for-europe-seminars", "plots", "regional_HS_change.png"),
    width = 10,
    height = 5.5
  )

# Country-level HS change map, facetted by GLW ---------------------------------------------------------------------------

country_HS_data <- HS_data %>% 
  select(iso = ISO3_Code, GLW, mean_HS_change) %>% 
  filter(!is.na(mean_HS_change)) %>% 
  mutate(GLW = factor(
    GLW, 
    levels = c("GLW_1.5", "GLW_2.0", "GLW_3.0", "GLW_4.0"), 
    labels = c("1.5°C", "2°C", "3°C", "4°C")
  ))

world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  # Natural Earth codes Palestine as PSX; match the ISO code used in the HS files
  mutate(iso = if_else(iso_a3 == "PSE", "PSE", adm0_a3)) %>% 
  filter(name != "Antarctica") %>% 
  select(iso, geometry)

country_HS_sf <- world %>% 
  inner_join(country_HS_data, by = "iso")

p <- ggplot() +
  geom_sf(data = world, fill = "grey80", colour = "grey45", linewidth = 0.1) +
  geom_sf(data = country_HS_sf, aes(fill = mean_HS_change), colour = "grey30", linewidth = 0.1) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
    # A few countries (e.g. NOR, PRK) have extreme values from near-zero reference HS, so clip the scale
    limits = c(-1, 1),
    oob = scales::squish,
    labels = scales::label_percent(),
    name = "Change (%)"
  ) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  facet_wrap(~GLW, ncol = 2) +
  labs(
    title = "Projected changes in habitat suitability within EEZs by country",
    subtitle = "Mean of each GLW across species, weighted by 2025 production of fed species"
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

ggsave(
  plot = p,
  filename = here("for-europe-seminars", "plots", "country_HS_change_map.png"),
  width = 12,
  height = 7
)
