# Description ---------------------------------------------------------------------------------------------------------
# Exploring Earth System Grid Federation data - to assess how it can be used to modify temperature with climate

library(tidyverse)
library(terra)
library(ncdf4)
library(sf)
library(here)

here("src", "dirs.R") %>% source()
here("src", "other_functions.R") %>% source()

# Extracting data -----------------------------------------------------------------------------------------------------
# Takes raw NC files from Earth Consortium and saves one raster per date. Rasters are gapfilled to compensate for the gaps reated by comforming the curvilinear grid to a flat raster. This means that there are some "land" cells which have SST values, but as long as you don't try to extract an SST value from a location on land it'll be fine.
# This process takes 10-12 minutes per file, there are 85 files (2015-2099), so ~14 hours in total.
tifs_nm <- "tos_Oday_EC-Earth3P-HR_highres-future_r3i1p2f1_gn_"
ncs_list <- file.path(bigdata_path, "ESGF", "EC-Earth3P-HR", "raw_ncs") %>% 
  list.files(full.names = T, recursive = T, pattern = tifs_nm) %>% 
  sort()
overwrite <- F

# Inspect the files
# highres_nc <- nc_open(ncs_list[1])
# highres_nc
# nc_close(highres_nc)

for (nc_file in ncs_list) {
  highres_nc <- nc_open(nc_file)

  # Get times
  time_raw  <- ncvar_get(highres_nc, "time")
  time_units <- ncatt_get(highres_nc, "time", "units")$value  
  # "days since 1900-01-01 00:00:00"
  origin <- as_datetime(sub("days since ", "", time_units))
  dates  <- as_date(origin + duration(time_raw, "days"))

  tif_nms <- file.path(bigdata_path, "ESGF", "EC-Earth3P-HR", "rasts", paste0(tifs_nm, dates, ".tif"))
  lyrs_todo <- if (overwrite) {
    seq_along(tif_nms)
  } else {
    which(!file.exists(tif_nms))
  }

  if (length(lyrs_todo) > 0) {
    message(sprintf("\nStarting file %s...\n", basename(nc_file)))

    lon  <- ncvar_get(highres_nc, "longitude")   # [1442, 1050]
    lat  <- ncvar_get(highres_nc, "latitude")
    # all_tos <- ncvar_get(highres_nc, "tos")
    # nc_close(highres_nc)

    template <- rast(resolution = 0.25, xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = "EPSG:4326")

    walk(
      lyrs_todo,
      function(lyr) { 
        df <- data.frame(
          lon = as.vector(lon),
          lat = as.vector(lat),
          tos = as.vector(ncvar_get(highres_nc, "tos", start = c(1,1,lyr), count = c(-1,-1,1)))
        )
        df$lon <- ifelse(df$lon > 180, df$lon - 360, df$lon)
        pts <- vect(df, geom = c("lon", "lat"), crs = "EPSG:4326")
        r <- rasterize(pts, template, field = "tos", fun = "mean")
        names(r) <- dates[lyr]
        # plot(r)

        # Gapfill the raster - note that this will give some (coastal) land cells SST values!
        gfr <- focal(
          r, 
          w = matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3), 
          fun = "mean", 
          na.policy = "only", 
          na.rm = T
        )
        # plot(gfr)

        writeRaster(gfr, tif_nms[lyr], overwrite = T)
      },
      .progress = T
    )

    nc_close(highres_nc)

  } else {
    message(sprintf("\nSkipping file %s, raster files already exist.\n", basename(nc_file)))
    nc_close(highres_nc)
  }
}

# Mean year 2015-2024 -------------------------------------------------------------------------------------------------
files <- file.path(bigdata_path, "ESGF", "EC-Earth3P-HR", "rasts") %>% 
  list.files(full.names = TRUE, recursive = T) %>% 
  sort() %>% 
  str_subset(paste0(2015:2024, collapse = "|"))

files_by_doy <- by_doy(files)

# Make sure directory exists
destpath <- file.path(bigdata_path, "ESGF", "EC-Earth3P-HR", "meanyears")
dir.create(destpath, showWarnings = F)

# Loop through each DOY and create average (gapfilled) SST raster
for (i in seq_len(nrow(files_by_doy))) {
  current_doy <- files_by_doy$doy_adjusted[i]
  current_date <- files_by_doy$month_day[i]
  current_files <- files_by_doy$filenamess[[i]]
  
  cat("\nProcessing DOY", current_doy, "(", current_date, ") - ", length(current_files), "files\n")
  
  raster_stack <- rast(current_files)
  mean_raster <- app(raster_stack, "mean", na.rm = T)

  # # Gapfill raster - WARNING: will give some land cells SST!
  # mean_gf_raster <- focal(
  #   mean_raster, 
  #   w = matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3), 
  #   fun = "mean", 
  #   na.policy = "only", 
  #   na.rm = T
  # )

  writeRaster(
    mean_raster, 
    file.path(
      destpath, 
      paste0("tos_Oday_EC-Earth3P-HR_highres-future_r3i1p2f1_gn_meanyear_2015_2024_", current_date, ".tif")
    ),
    datatype = "FLT4S",
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=3"), # compression for space
    overwrite = TRUE
  )
}
