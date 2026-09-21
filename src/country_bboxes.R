# ---------------------------------------------------------------------------
# country_bboxes.R
# Get bounding boxes (and a place to hang other per-country metadata) keyed by ISO 3166-1 alpha-3 code.
#
#   source("country_bboxes.R")
#   bbox("LBN", crs = 4326)
#
# HOW THE NUMBERS WERE MADE
# Extents come from Natural Earth 1:10m Admin 0 countries, padded by roughly 100 km on every side (0.9 deg of latitude; the longitude pad is widened by 1/cos(lat) so it is also ~100 km on the ground, capped at 20 deg near the poles).
#
# `bbox` covers the country's main landmass plus nearby islands, which is almost always what you want for a map. `bbox_full` is only present where the whole territory is meaningfully bigger (Chile + Easter Island, Ecuador + Galapagos, Spain + Canaries, USA + Alaska/Hawaii, ...).
#
# `label_position` is Natural Earth's hand-placed LABEL_X/LABEL_Y point, i.e. where a cartographer decided the country name should sit. For the 10 entries with no usable point of their own (the carved-out territories below, plus SGS and UMI) it is the centroid of the largest landmass; the bbox centre is the last resort and is currently unused. Note that these are label anchors, not centroids, so a few sit just offshore of long thin countries.
#
# Overseas territories with their own ISO3 code are split out into their own entries (GUF, GLP, MTQ, REU, MYT out of France; SJM, BVT out of Norway; BES out of the Netherlands), so FRA is metropolitan France.
#
# ANTIMERIDIAN
# For countries that straddle 180 deg (RUS, FJI, NZL, ...) xmax is allowed to exceed 180 rather than wrapping, e.g. Fiji is 173.61 -> 182.76. This keeps the box a real interval, which is what st_crop()/coord_sf() need; it is also why you should not assume xmax <= 180 anywhere downstream.
#
# EDITING
# These are hand-tweakable on purpose. Change a number, add a field, add a country; nothing here is derived at load time. Extra fields are free:
#   country_bbox_data$LBN$display_name <- "Lebanon"
#   country_bbox_data$LBN$rank         <- 3
# or just edit the list entry in place below.
#
# NOTE: `bbox()` will mask sp::bbox() / raster::bbox() if you use those.
# ---------------------------------------------------------------------------

.need_sf <- function() {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for bbox().", call. = FALSE)
  }
}

.bbox_order <- c("xmin", "ymin", "xmax", "ymax")

#' Pad a bare bbox vector by an approximate distance in km (negative shrinks)
pad_bbox <- function(b, km) {
  if (km == 0) return(b)
  dlat <- km / 111.32
  ymin <- max(-90, b[["ymin"]] - dlat)
  ymax <- min(90,  b[["ymax"]] + dlat)
  latref <- min(89, max(abs(ymin), abs(ymax)))
  dlon <- min(20, km / (111.32 * max(0.05, cos(latref * pi / 180))))
  xmin <- b[["xmin"]] - dlon
  xmax <- b[["xmax"]] + dlon
  if (xmax - xmin >= 360) {
    xmin <- -180
    xmax <- 180
  }
  c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
}

#' Everything stored for a country (or a named list of them)
#'
#' @param iso3 One or more ISO 3166-1 alpha-3 codes, case-insensitive.
country_info <- function(iso3, data = country_bbox_data) {
  iso3 <- toupper(iso3)
  unknown <- setdiff(iso3, names(data))
  if (length(unknown)) {
    stop("No entry for ISO3 code(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  if (length(iso3) == 1L) data[[iso3]] else data[iso3]
}

#' Bounding box for one or more countries
#'
#' @param iso3   ISO3 code(s). Several codes give the box that covers them all.
#' @param crs    Target CRS: an EPSG number, a proj/WKT string, or an
#'               `sf::st_crs()` object. The box is densified before it is
#'               reprojected, so the result is the true extent of the region
#'               rather than the corners joined by straight lines.
#' @param which  "bbox" (main landmass, the default) or "full" (whole
#'               territory, falling back to "bbox" where they are the same).
#'
#' @return An `sf::st_bbox` object.
#'
#' @examples
#'   bbox("LBN", crs = 4326)
#'   bbox("NOR", crs = 25833)                 # ETRS89 / UTM 33N, in metres
#'   bbox(c("LBN", "SYR", "ISR"))
#'   bbox("CHL", which = "full")              # includes Easter Island
bbox <- function(iso3,
                 crs    = 4326,
                 which  = c("bbox", "full"),
                 data   = country_bbox_data) {
  .need_sf()
  which <- match.arg(which)
  info  <- country_info(iso3, data = data)
  if (length(iso3) == 1L) info <- list(info)

  boxes <- lapply(info, function(x) {
    b <- if (which == "full" && !is.null(x$bbox_full)) x$bbox_full else x$bbox
    b[.bbox_order]
  })
  b <- c(
    xmin = min(vapply(boxes, `[[`, numeric(1), "xmin")),
    ymin = min(vapply(boxes, `[[`, numeric(1), "ymin")),
    xmax = max(vapply(boxes, `[[`, numeric(1), "xmax")),
    ymax = max(vapply(boxes, `[[`, numeric(1), "ymax"))
  )
  b <- pad_bbox(b, info$pad_km)

  target <- sf::st_crs(crs)
  if (is.na(target) || target == sf::st_crs(4326)) {
    return(sf::st_bbox(b, crs = sf::st_crs(4326)))
  }
  sf::st_bbox(sf::st_transform(.bbox_sfc(b), target))
}

#' The box as a (densified) polygon, e.g. for st_crop() or st_intersection()
bbox_poly <- function(iso3, crs = 4326, which = c("bbox", "full"), n = 50, data = country_bbox_data) {
  .need_sf()
  b <- bbox(iso3, crs = 4326, which = match.arg(which), data = data)
  p <- .bbox_sfc(b, n = n)
  target <- sf::st_crs(crs)
  if (is.na(target) || target == sf::st_crs(4326)) p else sf::st_transform(p, target)
}

#' xlim/ylim ready to splice into coord_sf()
#'
#'   ggplot(x) + geom_sf() + do.call(coord_sf, bbox_lims("LBN"))
bbox_lims <- function(iso3, crs = 4326, which = c("bbox", "full"), data = country_bbox_data) {
  b <- bbox(iso3, crs = crs, which = match.arg(which), data = data)
  list(xlim = c(b[["xmin"]], b[["xmax"]]),
       ylim = c(b[["ymin"]], b[["ymax"]]),
       expand = FALSE)
}

#' Label anchor point(s) as an sf data frame
#'
#' Columns iso3, name and geometry, so it drops straight into a plot:
#'
#'   ggplot(map) +
#'     geom_sf() +
#'     geom_sf_text(data = label_point(c("LBN", "SYR")), aes(label = name))
#'
#' @param crs Target CRS; the points are reprojected from EPSG:4326.
label_point <- function(iso3, crs = 4326, data = country_bbox_data) {
  .need_sf()
  iso3 <- toupper(iso3)
  info <- country_info(iso3, data = data)
  if (length(iso3) == 1L) info <- list(info)

  pts <- sf::st_sfc(
    lapply(info, function(x) sf::st_point(unname(x$label_position[c("lon", "lat")]))),
    crs = sf::st_crs(4326)
  )
  target <- sf::st_crs(crs)
  if (!is.na(target) && target != sf::st_crs(4326)) pts <- sf::st_transform(pts, target)

  sf::st_sf(
    iso3 = iso3,
    name = vapply(info, function(x) x$name, character(1)),
    geometry = pts
  )
}

# Rectangle in EPSG:4326 with n points along each edge, so that reprojecting it
# follows the curved graticule instead of cutting the corners.
.bbox_sfc <- function(b, n = 50) {
  .need_sf()
  x <- seq(b[["xmin"]], b[["xmax"]], length.out = n)
  y <- seq(b[["ymin"]], b[["ymax"]], length.out = n)
  ring <- rbind(
    cbind(head(x, -1), b[["ymin"]]),
    cbind(b[["xmax"]], head(y, -1)),
    cbind(head(rev(x), -1), b[["ymax"]]),
    cbind(b[["xmin"]], head(rev(y), -1))
  )
  ring <- rbind(ring, ring[1, ])
  sf::st_sfc(sf::st_polygon(list(unname(ring))), crs = sf::st_crs(4326))
}

# ---------------------------------------------------------------------------
# Data: 247 entries, alphabetical by ISO3.
# xmin/ymin/xmax/ymax are decimal degrees, EPSG:4326.
# ---------------------------------------------------------------------------
#' Load country bounding-box data from CSV into the original nested list
#'
#' Any columns beyond the standard set are attached to each country's entry
#' under their own column name, e.g. a `custom_name` column becomes
#' `result[["ABW"]][["custom_name"]]`.
#'
#' @param path Path to the CSV.
#' @return A named list keyed by ISO3 code.
load_country_bbox_data <- function(
  path = file.path(rawdata_path, "country_bbox_data.csv")
) {
  bbox_parts <- c("xmin", "ymin", "xmax", "ymax")
  bbox_cols <- paste0("bbox_", bbox_parts)
  bbox_full_cols <- paste0("bbox_full_", bbox_parts)
  label_cols <- c("label_lon", "label_lat")

  required_cols <- c("iso3", "name", "iso2", bbox_cols, label_cols)
  core_cols <- c(required_cols, bbox_full_cols)

  # na = "" so Namibia's ISO2 code "NA" is kept as a string
  df <- readr::read_csv(
    path,
    na = "",
    col_types = readr::cols(
      iso3 = readr::col_character(),
      name = readr::col_character(),
      iso2 = readr::col_character(),
      !!!rlang::set_names(
        rep(list(readr::col_double()), length(c(bbox_cols, bbox_full_cols, label_cols))),
        c(bbox_cols, bbox_full_cols, label_cols)
      ),
      .default = readr::col_guess()
    ),
    show_col_types = FALSE
  )

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("CSV is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  if (anyDuplicated(df$iso3)) {
    stop("Duplicate iso3 codes: ", paste(unique(df$iso3[duplicated(df$iso3)]), collapse = ", "))
  }

  has_full_cols <- all(bbox_full_cols %in% names(df))
  extra_cols <- setdiff(names(df), core_cols)

  get_vec <- function(i, cols, nms) {
    stats::setNames(vapply(cols, \(col) df[[col]][[i]], numeric(1)), nms)
  }

  entries <- lapply(seq_len(nrow(df)), function(i) {
    entry <- list(
      name = df$name[[i]],
      iso2 = df$iso2[[i]],
      bbox = get_vec(i, bbox_cols, bbox_parts)
    )

    if (has_full_cols) {
      bbox_full <- get_vec(i, bbox_full_cols, bbox_parts)
      if (!all(is.na(bbox_full))) entry$bbox_full <- bbox_full
    }

    entry$label_position <- get_vec(i, label_cols, c("lon", "lat"))

    for (col in extra_cols) {
      entry[[col]] <- df[[col]][[i]]
    }

    entry
  })

  stats::setNames(entries, df$iso3)
}

country_bbox_data <- load_country_bbox_data()
