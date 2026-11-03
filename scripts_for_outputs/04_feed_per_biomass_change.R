# Dot-and-line plot of projected change in feed required per unit harvested
# biomass per model at each global warming level (GLW), one png per focus country.
# Assumes scripts_for_outputs/00_run_all.R has already been run (provides
# `assigned_farms`, `focus_countries`, `GLW_ls`, `GLW_labels`, `report_theme()`,
# and paths from src/dirs.R).

library(sf)
library(here)

here("src", "species_colours.R") %>% source()    # models_lab
here("src", "other_functions.R") %>% source()    # get_country_data2(), pct_limits(), pct_breaks()
here("src", "country_bboxes.R")  %>% source()    # country_info()

warming_cols <- c("1.5°C" = "#E9A23B", "2°C" = "#DB5F32", "3°C" = "#B3243B", "4°C" = "#5C0F2E")
offsets      <- c("1.5°C" = -0.24, "2°C" = -0.08, "3°C" = 0.08, "4°C" = 0.24)

feed_per_biomass_change_plot <- function(this_iso, out_dir) {
  these_farms <- assigned_farms %>%
    st_drop_geometry() %>%
    filter(ISO3_Code == this_iso)

  production_files <- file.path(
    outs_path, "data", "farm_production_consolidated",
    paste0(this_iso, "_", unique(these_farms$model_name), ".qs")
  )
  missing <- production_files[!file.exists(production_files)]
  if (length(missing) > 0) {
    warning(this_iso, ": no consolidated production file for ", paste(basename(missing), collapse = ", "))
    production_files <- setdiff(production_files, missing)
  }
  if (length(production_files) == 0) {
    warning(this_iso, ": no production data, skipping plot")
    return(invisible(NULL))
  }

  farm_data <- get_country_data2(production_files) %>%
    semi_join(these_farms, by = join_by(farm_ID))

  # Models ordered by mean 2025 biomass, largest first (same order as the biomass plot)
  model_order <- farm_data %>%
    group_by(model_name) %>%
    reframe(biomass_2025 = mean(biomass_2025)) %>%
    arrange(-biomass_2025) %>%
    pull(model_name)

  # Mean change per farm within each GLW window, then median + IQR across farms per model
  plt_data <- map2_dfr(
    GLW_ls, names(GLW_ls),
    function(yrs, nm) {
      farm_data %>%
        filter(year %in% yrs) %>%
        group_by(farm_ID, model_name) %>%
        reframe(feed_prov_biom_change_pct = mean(feed_prov_biom_change_pct, na.rm = T)) %>%
        mutate(GLW = nm)
    }
  ) %>%
    mutate(GLW = factor(GLW, levels = names(GLW_labels), labels = GLW_labels)) %>%
    group_by(model_name, GLW) %>%
    reframe(
      median = median(feed_prov_biom_change_pct, na.rm = T),
      q25 = quantile(feed_prov_biom_change_pct, 0.25, na.rm = T),
      q75 = quantile(feed_prov_biom_change_pct, 0.75, na.rm = T)
    ) %>%
    mutate(
      model_name = factor(model_name, levels = rev(model_order)),
      y = as.numeric(model_name) + offsets[as.character(GLW)]
    )

  plt_limits <- pct_limits(c(plt_data$q25, plt_data$q75))
  plt_breaks <- pct_breaks(c(plt_data$q25, plt_data$q75))

  p <- plt_data %>%
    ggplot(aes(x = y, y = median, colour = GLW)) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.5, linetype = "dashed") +
    geom_segment(aes(y = q25, yend = q75), linewidth = 0.8) +
    geom_point(aes(y = median), size = 2.3) +
    scale_colour_manual(values = warming_cols, name = "Global warming level") +
    scale_y_continuous(
      breaks = plt_breaks,
      labels = scales::label_percent(),
      limits = plt_limits,
      expand = expansion(0)
    ) +
    scale_x_continuous(
      breaks = seq_along(rev(model_order)),
      labels = models_lab[rev(model_order)],
      expand = expansion(add = 0.5)
    ) +
    coord_flip() +
    labs(
      # title = paste0("Projected change in feed per unit biomass: ", country_info(this_iso)$name, " (", this_iso, ")"),
      y = "Projected change (%)",
      x = NULL
    ) +
    guides(colour = guide_legend(override.aes = list(linewidth = 0, size = 3))) +
    report_theme()

  ggsave(
    filename = file.path(out_dir, paste0("feed_per_biomass_change_GLW_", this_iso, ".png")),
    plot = p, width = 8.5, height = GLW_plot_height(length(model_order)), dpi = 300, bg = "white"
  )
}

walk2(focus_countries$ISO3_Code, focus_countries$dir_name, feed_per_biomass_change_plot)
