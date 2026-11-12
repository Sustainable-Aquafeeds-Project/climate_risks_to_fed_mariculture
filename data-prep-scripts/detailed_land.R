# The purpose of this script is to take detailed polygons of coastlines and save them into more manageable files for each country.
# Unlike most of the data prep scripts, this one was conceived of quite late in the process so requires that many other variables/scripts already exist.

library(sf)
library(here)
library(tidyverse)
library(qs2)
library(future)
library(furrr)

here("src/dirs.R") %>% source()
here("src/country_bboxes.R") %>% source()

land_path <- file.path(prepdata_path, "detailed_land_polys")
coast_path <- file.path(prepdata_path, "coastal_zone_polys")
coast2_path <- file.path(prepdata_path, "coastal_zone_50_polys")
dir.create(land_path, showWarnings = F)
dir.create(coast_path, showWarnings = F)
dir.create(coast2_path, showWarnings = F)

assigned_farms_final <- file.path(prepdata_path, "assigned_farms_final.qs") %>% qs_read()
isos <- assigned_farms_final %>% pull(ISO3_Code) %>% unique() %>% sort()

land_files <- file.path(land_path, paste0(isos, ".qs"))
coast_files <- file.path(coast_path, paste0(isos, ".qs"))
coast2_files <- file.path(coast2_path, paste0(isos, ".qs"))

# Only create files that don't already exist (delete the files to regenerate)
land_isos <- isos[!file.exists(land_files)]
coast_isos <- isos[!file.exists(coast_files)]
coast2_isos <- isos[!file.exists(coast2_files)]

rm(assigned_farms_final)
gc()

plan(multisession, workers = 10)

# 1. Detailed land polygons ----
plan(multisession, workers = 10)
future_walk(
  land_isos,
  function(this_iso) {
    sf_use_s2(FALSE)
    this_bbox <- bbox_lims(this_iso)
    this_bbox <- st_bbox(c(
      xmin = this_bbox[['xlim']][1], xmax = this_bbox[['xlim']][2], 
      ymin = this_bbox[['ylim']][1], ymax = this_bbox[['ylim']][2]
    ), crs = st_crs(4326))

    this_land <- st_read(
      file.path(rawdata_path, "land-polygons-split-4326", "land_polygons.shp"),
      wkt_filter = st_as_text(st_as_sfc(this_bbox)),
      quiet = T
    )
    this_land <- st_make_valid(this_land)
    # crop before union so the union only works on what we keep
    this_land <- st_crop(this_land, this_bbox)
    this_land <- st_union(this_land)          # merge the tiles the file is split into
    this_land <- st_make_valid(this_land)

    qs_save(this_land, file.path(land_path, paste0(this_iso, ".qs")))
    sf_use_s2(TRUE)
  },
  .progress = T,
  .options = furrr_options(seed = T)
)

# 2. Coastal zone polygons ----
# Buffer around the land polygons (in degrees, since s2 is off), read from the files saved above 
make_coast_1 <- function(this_iso, dist = 0.75, tol = dist / 10) {
  sf_use_s2(FALSE)
  on.exit(sf_use_s2(TRUE), add = TRUE)

  land <- file.path(land_path, paste0(this_iso, ".qs")) |>
    qs_read() |>
    st_geometry() |>
    st_collection_extract("POLYGON") %>% 
    st_cast("POLYGON")                       # one geometry per landmass/island

  # per-polygon bbox: xmin, ymin, xmax, ymax
  bbs   <- t(vapply(land, \(p) as.numeric(st_bbox(p)), numeric(4)))
  diag  <- sqrt((bbs[, 3] - bbs[, 1])^2 + (bbs[, 4] - bbs[, 2])^2)
  small <- diag < tol

  # big polygons: simplify as before
  big <- st_simplify(land[!small], dTolerance = tol, preserveTopology = TRUE)

  # small islands: bbox centre, snapped to the tolerance grid, deduplicated
  pieces <- big
  if (any(small)) {
    ctr <- cbind((bbs[small, 1] + bbs[small, 3]) / 2,
                 (bbs[small, 2] + bbs[small, 4]) / 2)
    ctr <- unique(round(ctr / tol) * tol)
    pts <- st_cast(st_sfc(st_multipoint(ctr), crs = st_crs(land)), "POINT")
    pieces <- c(big, pts)
  }
  rm(land, big); gc()

  this_coast <- pieces |>
    st_buffer(dist, nQuadSegs = 5) |>        # many small, independent buffers
    st_union() |>                            # STRtree-based cascaded union
    st_make_valid()

  qs_save(this_coast, file.path(coast_path, paste0(this_iso, ".qs")))
  invisible(NULL)
}

# biggest countries sequentially, the rest in parallel
sizes     <- file.size(land_files)
big_lands <- isos[sizes > quantile(sizes, 0.9)]
big_isos  <- coast_isos[coast_isos %in% big_lands]

purrr::walk(big_isos, make_coast_1, .progress = TRUE)
future_walk(setdiff(coast_isos, big_isos), make_coast_1, .progress = TRUE)

# 3. Coastal zone 2 polygons ----
# Buffer ~62.5 km around the land polygons. A fixed buffer in degrees is a different ground distance at every latitude (1 deg longitude ~ 111 km * cos(lat)), so instead project each country to a local azimuthal equidistant CRS centred on it, buffer in metres, then transform back to lon/lat.

# After st_transform() back to 4326, put longitudes in [lon0-180, lon0+180]
# so nothing jumps across the antimeridian, then split at ±180 and shift back.
unwrap_lon <- function(x, lon0) {
  fix <- function(m) { m[, 1] <- lon0 + ((m[, 1] - lon0 + 180) %% 360) - 180; m }
  x <- st_sfc(
    lapply(x, \(g) structure(rapply(unclass(g), fix, how = "replace"), class = class(g))),
    crs = 4326
  )
  bb <- st_bbox(x)
  if (bb[["xmin"]] >= -180 && bb[["xmax"]] <= 180) return(x)

  parts <- lapply(c(-360, 0, 360), \(dx) {
    win <- st_as_sfc(st_bbox(c(xmin = -180 + dx, ymin = -90, xmax = 180 + dx, ymax = 90),
                             crs = st_crs(4326)))
    p <- suppressMessages(st_intersection(x, win))
    p <- p[!st_is_empty(p)]
    if (length(p) == 0) return(NULL)
    st_set_crs(p - c(dx, 0), 4326)
  })
  do.call(c, Filter(Negate(is.null), parts))
}

make_coast_m <- function(this_iso, dist_m = 62000, tol_m = dist_m / 10, cell_deg = 5, nQuadSegs = 5) {
  sf_use_s2(FALSE)
  on.exit(sf_use_s2(TRUE), add = TRUE)

  land <- file.path(land_path, paste0(this_iso, ".qs")) |>
    qs_read() |>
    st_geometry() |>
    st_collection_extract("POLYGON") %>% 
    st_cast("POLYGON") |>
    # a tolerance of tol_m / 111320 degrees can never move a vertex more than tol_m,
    # at any latitude, so this is a safe, cheap pre-thin before cutting
    st_simplify(dTolerance = tol_m / 111320, preserveTopology = TRUE)

  # --- cut into grid cells -------------------------------------------------
  bb   <- st_bbox(land)
  grid <- st_make_grid(land, cellsize = cell_deg,
                       offset = floor(bb[c("xmin", "ymin")] / cell_deg) * cell_deg)
  land <- suppressWarnings(suppressMessages(
    st_intersection(st_sf(geometry = land), st_sf(cell = seq_along(grid), geometry = grid))
  ))
  land <- land[st_dimension(land) == 2 & !st_is_empty(land), ]
  land <- suppressWarnings(st_collection_extract(land, "POLYGON"))

  # --- buffer each cell in its own local projection ------------------------
  cell_bufs <- lapply(split(st_geometry(land), land$cell), function(g) {
    i    <- as.integer(attr(g, "cell_id"))  # placeholder, set below
    g
  })
  cell_ids <- as.integer(names(cell_bufs))

  cell_bufs <- Map(function(g, i) {
    cb   <- st_bbox(grid[i])
    lon0 <- (cb[["xmin"]] + cb[["xmax"]]) / 2
    lat0 <- (cb[["ymin"]] + cb[["ymax"]]) / 2
    crs_loc <- st_crs(sprintf("+proj=aeqd +lat_0=%.6f +lon_0=%.6f +datum=WGS84 +units=m",
                              lat0, lon0))

    g <- g |> st_transform(crs_loc) |> st_cast("MULTIPOLYGON") |> st_cast("POLYGON")

    # small islands -> snapped points (in metres now)
    bbs   <- t(vapply(g, \(p) as.numeric(st_bbox(p)), numeric(4)))
    diag  <- sqrt((bbs[, 3] - bbs[, 1])^2 + (bbs[, 4] - bbs[, 2])^2)
    small <- diag < tol_m

    pieces <- st_simplify(g[!small], dTolerance = tol_m, preserveTopology = TRUE)
    if (any(small)) {
      ctr <- cbind((bbs[small, 1, drop = FALSE] + bbs[small, 3, drop = FALSE]) / 2,
                   (bbs[small, 2, drop = FALSE] + bbs[small, 4, drop = FALSE]) / 2)
      ctr <- unique(round(ctr / tol_m) * tol_m)
      pts <- st_cast(st_sfc(st_multipoint(ctr), crs = crs_loc), "POINT")
      pieces <- c(pieces, pts)
    }

    pieces |>
      st_buffer(dist_m, nQuadSegs = nQuadSegs) |>
      st_union() |>
      st_transform(4326) |>
      unwrap_lon(lon0)
  }, cell_bufs, cell_ids)

  rm(land); gc()

  this_coast <- do.call(c, cell_bufs) |>
    st_union() |>
    st_make_valid()

  qs_save(this_coast, file.path(coast2_path, paste0(this_iso, ".qs")))
  invisible(NULL)
}

big2_isos  <- coast2_isos[coast2_isos %in% big_lands]

purrr::walk(big2_isos, make_coast_m, .progress = TRUE)
future_walk(setdiff(coast2_isos, big2_isos), make_coast_m, .progress = TRUE)
plan(sequential)
