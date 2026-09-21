# Correcting temperatures
# =======================
# Uses the change in temperature per grid cell from CMIP data to inform climate change impacts, starting from contemporary NASA SST data.
# Bias correction uses the decadal mean from 2015-2024 for both datasets.
#
# Outputs:
#   - bigdata_path/ESGF/EC-Earth3P-HR/projected_adjusted/   (bias-corrected rasters)
#   - data_path/intermediates/comparing-corrected-temps.qs  (diagnostic comparison)
#   - prepdata_path/global-temperature-changes.qs           (global mean temps)

library(tidyverse)
library(terra)
library(qs2)
library(here)

here("src", "dirs.R") %>% source()
here("src", "other_functions.R") %>% source()

# ---- Meanyear file listings ----

nasa_fnms <- file.path(bigdata_path, "MUR25-JPL-L4-GLOB-v04.2_4.2", "meanyears") %>%
  list.files(full.names = T, pattern = "2015_2024")
cmip_fnms <- file.path(bigdata_path, "ESGF", "EC-Earth3P-HR", "meanyears") %>%
  list.files(full.names = T, pattern = "2015_2024")

# ---- NASA meanyear temperatures (2015-2024) ----

nasa_nms <- nasa_fnms %>%
  basename() %>%
  str_remove_all(".tif") %>%
  str_split_i("_", 5) %>%
  paste0("2015-", .) %>%
  as.Date() %>% format("%j") %>% paste0("DOY_", .)

nasa_meanyear <- rast(nasa_fnms)
names(nasa_meanyear) <- nasa_nms

# Quick check - plot DOY 100
plot(rast(nasa_fnms[100]))

# ---- CMIP meanyear temperatures (2015-2024) ----

cmip_nms <- cmip_fnms %>%
  basename() %>%
  str_remove_all(".tif") %>%
  str_split_i("_", 10) %>%
  paste0("2015-", .) %>%
  as.Date() %>% format("%j") %>% paste0("DOY_", .)

cmip_meanyear <- rast(cmip_fnms)
names(cmip_meanyear) <- cmip_nms

# Quick check - plot DOY 100 (should be gapfilled)
plot(rast(cmip_fnms[100]))

# ---- Apply delta correction to all future CMIP projections ----
# For each future date, the corrected temperature is:
#   corrected = CMIP_future - CMIP_meanyear + NASA_meanyear

cmip_fnms <- file.path(bigdata_path, "ESGF", "EC-Earth3P-HR", "rasts") %>%
  list.files(full.names = T) # %>%
  # str_subset(paste0("_gn_", 2015:2024, collapse = "|"), negate = T)

# All future dates needing to be adjusted
future_dates <- as.Date(
  str_extract(cmip_fnms, "\\d{4}-\\d{2}-\\d{2}"),
  format = "%Y-%m-%d"
)

date_df <- tibble(
  date  = future_dates,
  fnm   = cmip_fnms,
  feb29 = leap_year(future_dates) & month(future_dates) == 2 & day(future_dates) == 29
) |>
  mutate(
    doy = yday(date),
    # Shift post-Feb-28 leap year dates back by 1 to align with 365-day meanyear
    doy = if_else(!feb29 & leap_year(date) & doy > 60, doy - 1L, doy),
    # Feb 29 gets a special group key; all others group by adjusted DOY
    group_key = if_else(feb29, "feb29", str_pad(doy, 3, pad = "0"))
  )
date_ls <- split(date_df, date_df$group_key)

# Where to save files
out_dir <- file.path(bigdata_path, "ESGF", "EC-Earth3P-HR", "projected_adjusted")
dir.create(out_dir, showWarnings = F)

# Process all dates sharing the same climatology layer together
# grp <- date_ls[[100]]
walk(
  date_ls,
  function(grp) {
    cmip_fnames <- file.path(out_dir, basename(grp$fnm))
    if (any(!file.exists(cmip_fnames))) {
      if (grp$feb29[1]) {
        cmip_clim_lyr <- weighted.mean(cmip_meanyear[[ c("DOY_059", "DOY_060") ]], c(0.5, 0.5))
        nasa_clim_lyr <- weighted.mean(nasa_meanyear[[  c("DOY_059", "DOY_060") ]], c(0.5, 0.5))
      } else {
        lyr_name      <- paste0("DOY_", grp$group_key[1])
        cmip_clim_lyr <- cmip_meanyear[[ lyr_name ]]
        nasa_clim_lyr <- nasa_meanyear[[ lyr_name ]]
      }

      # Filter to the files that actually need to be created
      cmip_fnames_todo <- cmip_fnames[!file.exists(cmip_fnames)]
      grp_todo <- grp[!file.exists(cmip_fnames), ]

      # Stack all CMIP rasters for this DOY and apply delta in one operation
      cmip_stack        <- rast(grp_todo$fnm)
      projected         <- cmip_stack - cmip_clim_lyr + nasa_clim_lyr
      names(cmip_stack) <- grp_todo$date

      # Write out individual files, preserving original filenames
      walk2(
        seq_len(nlyr(projected)),
        cmip_fnames_todo,
        function(j, fnm) {
          writeRaster(projected[[j]], fnm, overwrite = TRUE)
        })
    }
  },
  .progress = T
)

# ---- Compare corrected vs. uncorrected global temperatures (diagnostic) ----
corr_files <- out_dir %>%
  list.files(full.names = T)
uncorr_files <- file.path(bigdata_path, "ESGF", "EC-Earth3P-HR", "rasts") %>%
  list.files(full.names = T) # %>%
  # str_subset(paste0(2015:2024, collapse = "|"), negate = T)

dates <- str_extract(uncorr_files, "\\d{4}-\\d{2}-\\d{2}")

file_list <- setNames(
  Map(list, uncorr = uncorr_files, corr = corr_files),
  dates
)

diffs <- map2_dfr(
  file_list,
  names(file_list),
  function(ls, nm) {
    uncorr <- rast(ls[["uncorr"]])
    corr   <- rast(ls[["corr"]])
    diff   <- uncorr - corr
    diff_vals <- values(diff)
    diff_vals <- diff_vals[!is.na(diff_vals)]
    data.frame(
      date = nm,
      mean = mean(diff_vals),
      min  = min(diff_vals),
      max  = max(diff_vals),
      med  = median(diff_vals)
    )
  },
  .progress = T
)

diffs <- diffs %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"))

dir.create(file.path(data_path, "intermediates"), showWarnings = FALSE)
qs_save(diffs, file.path(data_path, "intermediates", "comparing-corrected-temps.qs"))

# Plot the comparison
diffs %>%
  pivot_longer(cols = -date, names_to = "measure", values_to = "value") %>%
  ggplot(aes(x = date, y = value)) +
  geom_line(linewidth = 0.75) +
  facet_wrap(~measure, scales = "free", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

# ---- Global mean temperature summary across corrected projections ----

glob_temp <- map2_dfr(
  corr_files,
  dates,
  function(fnm, nm) {
    r <- rast(fnm)
    vals <- values(r)
    vals <- vals[!is.na(vals)]
    data.frame(
      date = nm,
      mean = mean(vals),
      sd   = sd(vals),
      min  = min(vals),
      max  = max(vals),
      med  = median(vals)
    )
  },
  .progress = T
)

qs_save(glob_temp, file.path(prepdata_path, "global-temperature-changes.qs"))

# Plot global mean temperature change relative to 2025
glob_temp <- qs_read(file.path(prepdata_path, "global-temperature-changes.qs")) %>%
  mutate(date = as.Date(date)) %>%
  mutate(yr = year(date)) %>%
  group_by(yr) %>%
  reframe(mean = mean(mean))

glob_temp$diff_mean <- glob_temp$mean - glob_temp$mean[glob_temp$yr == 2025]

glob_temp %>%
  ggplot(aes(x = yr, y = diff_mean)) +
  geom_line() +
  theme_classic()
