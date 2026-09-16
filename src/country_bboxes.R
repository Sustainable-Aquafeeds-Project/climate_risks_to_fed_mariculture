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
#' @param pad_km Extra padding in km on top of the ~100 km already baked in.
#'
#' @return An `sf::st_bbox` object.
#'
#' @examples
#'   bbox("LBN", crs = 4326)
#'   bbox("NOR", crs = 25833)                 # ETRS89 / UTM 33N, in metres
#'   bbox(c("LBN", "SYR", "ISR"), pad_km = 50)
#'   bbox("CHL", which = "full")              # includes Easter Island
bbox <- function(iso3,
                 crs    = 4326,
                 which  = c("bbox", "full"),
                 pad_km = 0,
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
  b <- pad_bbox(b, pad_km)

  target <- sf::st_crs(crs)
  if (is.na(target) || target == sf::st_crs(4326)) {
    return(sf::st_bbox(b, crs = sf::st_crs(4326)))
  }
  sf::st_bbox(sf::st_transform(.bbox_sfc(b), target))
}

#' The box as a (densified) polygon, e.g. for st_crop() or st_intersection()
bbox_poly <- function(iso3, crs = 4326, which = c("bbox", "full"),
                      pad_km = 0, n = 50, data = country_bbox_data) {
  .need_sf()
  b <- bbox(iso3, crs = 4326, which = match.arg(which), pad_km = pad_km,
            data = data)
  p <- .bbox_sfc(b, n = n)
  target <- sf::st_crs(crs)
  if (is.na(target) || target == sf::st_crs(4326)) p else sf::st_transform(p, target)
}

#' xlim/ylim ready to splice into coord_sf()
#'
#'   ggplot(x) + geom_sf() + do.call(coord_sf, bbox_lims("LBN"))
bbox_lims <- function(iso3, crs = 4326, which = c("bbox", "full"),
                      pad_km = 0, data = country_bbox_data) {
  b <- bbox(iso3, crs = crs, which = match.arg(which), pad_km = pad_km,
            data = data)
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

country_bbox_data <- list(
  ABW = list(
    name = "Aruba", iso2 = "AW",
    bbox = c(xmin = -70.99, ymin = 11.51, xmax = -68.95, ymax = 13.54),
    label_position = c(lon = -69.9728, lat = 12.5174)
  ),
  AFG = list(
    name = "Afghanistan", iso2 = "AF",
    bbox = c(xmin = 59.32, ymin = 28.48, xmax = 76.06, ymax = 39.38),
    label_position = c(lon = 66.4966, lat = 34.1643)
  ),
  AGO = list(
    name = "Angola", iso2 = "AO",
    bbox = c(xmin = 10.71, ymin = -18.93, xmax = 25.02, ymax = -3.49),
    label_position = c(lon = 17.9842, lat = -12.1828)
  ),
  AIA = list(
    name = "Anguilla", iso2 = "AI",
    bbox = c(xmin = -64.39, ymin = 17.27, xmax = -61.97, ymax = 19.5),
    label_position = c(lon = -63.0264, lat = 18.243)
  ),
  ALA = list(
    name = "Aland", iso2 = "AX",
    bbox = c(xmin = 17.63, ymin = 59, xmax = 22.98, ymax = 61.38),
    label_position = c(lon = 19.8697, lat = 60.1565)
  ),
  ALB = list(
    name = "Albania", iso2 = "AL",
    bbox = c(xmin = 18.03, ymin = 38.73, xmax = 22.28, ymax = 43.56),
    label_position = c(lon = 20.1138, lat = 40.6549)
  ),
  AND = list(
    name = "Andorra", iso2 = "AD",
    bbox = c(xmin = 0.16, ymin = 41.53, xmax = 3.01, ymax = 43.55),
    label_position = c(lon = 1.5394, lat = 42.5476)
  ),
  ARE = list(
    name = "United Arab Emirates", iso2 = "AE",
    bbox = c(xmin = 50.56, ymin = 21.72, xmax = 57.4, ymax = 26.98),
    label_position = c(lon = 54.5473, lat = 23.4663)
  ),
  ARG = list(
    name = "Argentina", iso2 = "AR",
    bbox = c(xmin = -75.18, ymin = -55.96, xmax = -52.05, ymax = -20.88),
    label_position = c(lon = -64.1733, lat = -33.5012)
  ),
  ARM = list(
    name = "Armenia", iso2 = "AM",
    bbox = c(xmin = 42.22, ymin = 37.96, xmax = 47.82, ymax = 42.19),
    label_position = c(lon = 44.8006, lat = 40.4591)
  ),
  ASM = list(
    name = "American Samoa", iso2 = "AS",
    bbox = c(xmin = -172.02, ymin = -15.44, xmax = -167.22, ymax = -10.15),
    label_position = c(lon = -170.7472, lat = -14.3267)
  ),
  ATA = list(
    name = "Antarctica", iso2 = "AQ",
    bbox = c(xmin = -180, ymin = -90, xmax = 180, ymax = -59.61),
    label_position = c(lon = 35.8855, lat = -79.8432)
  ),
  ATF = list(
    name = "French Southern and Antarctic Lands", iso2 = "TF",
    bbox = c(xmin = 67.19, ymin = -50.62, xmax = 71.99, ymax = -47.66),
    bbox_full = c(xmin = 38.31, ymin = -50.62, xmax = 79.01, ymax = -10.65),
    label_position = c(lon = 69.1221, lat = -49.3037)
  ),
  ATG = list(
    name = "Antigua and Barbuda", iso2 = "AG",
    bbox = c(xmin = -63.3, ymin = 16.03, xmax = -60.71, ymax = 18.63),
    label_position = c(lon = -61.7906, lat = 17.3522)
  ),
  AUS = list(
    name = "Australia", iso2 = "AU",
    bbox = c(xmin = 111.65, ymin = -44.54, xmax = 154.9, ymax = -8.34),
    bbox_full = c(xmin = 111.32, ymin = -55.65, xmax = 160.7, ymax = -8.34),
    label_position = c(lon = 134.0497, lat = -24.1295)
  ),
  AUT = list(
    name = "Austria", iso2 = "AT",
    bbox = c(xmin = 8.12, ymin = 45.48, xmax = 18.55, ymax = 49.91),
    label_position = c(lon = 14.1305, lat = 47.5189)
  ),
  AZE = list(
    name = "Azerbaijan", iso2 = "AZ",
    bbox = c(xmin = 43.55, ymin = 37.49, xmax = 51.85, ymax = 42.79),
    label_position = c(lon = 47.211, lat = 40.4024)
  ),
  BDI = list(
    name = "Burundi", iso2 = "BI",
    bbox = c(xmin = 28.08, ymin = -5.37, xmax = 31.74, ymax = -1.4),
    label_position = c(lon = 29.9171, lat = -3.3328)
  ),
  BEL = list(
    name = "Belgium", iso2 = "BE",
    bbox = c(xmin = 1.04, ymin = 48.59, xmax = 7.85, ymax = 52.4),
    label_position = c(lon = 4.8004, lat = 50.7854)
  ),
  BEN = list(
    name = "Benin", iso2 = "BJ",
    bbox = c(xmin = -0.17, ymin = 5.31, xmax = 4.77, ymax = 13.3),
    label_position = c(lon = 2.352, lat = 10.3248)
  ),
  BES = list(
    name = "Caribbean Netherlands", iso2 = "BQ",
    bbox = c(xmin = -69.35, ymin = 11.12, xmax = -67.26, ymax = 13.21),
    bbox_full = c(xmin = -69.37, ymin = 11.12, xmax = -61.99, ymax = 18.55),
    label_position = c(lon = -68.2876, lat = 12.1865)
  ),
  BFA = list(
    name = "Burkina Faso", iso2 = "BF",
    bbox = c(xmin = -6.46, ymin = 8.49, xmax = 3.33, ymax = 15.98),
    label_position = c(lon = -1.3639, lat = 12.673)
  ),
  BGD = list(
    name = "Bangladesh", iso2 = "BD",
    bbox = c(xmin = 87, ymin = 19.84, xmax = 93.66, ymax = 27.53),
    label_position = c(lon = 89.685, lat = 24.215)
  ),
  BGR = list(
    name = "Bulgaria", iso2 = "BG",
    bbox = c(xmin = 21.07, ymin = 40.33, xmax = 29.88, ymax = 45.13),
    label_position = c(lon = 25.1571, lat = 42.5088)
  ),
  BHR = list(
    name = "Bahrain", iso2 = "BH",
    bbox = c(xmin = 49.37, ymin = 24.68, xmax = 51.83, ymax = 27.19),
    label_position = c(lon = 50.5548, lat = 26.056)
  ),
  BHS = list(
    name = "The Bahamas", iso2 = "BS",
    bbox = c(xmin = -80.62, ymin = 20.01, xmax = -71.73, ymax = 27.83),
    label_position = c(lon = -77.1467, lat = 26.4018)
  ),
  BIH = list(
    name = "Bosnia and Herzegovina", iso2 = "BA",
    bbox = c(xmin = 14.41, ymin = 41.66, xmax = 20.92, ymax = 46.19),
    label_position = c(lon = 18.0684, lat = 44.0911)
  ),
  BLM = list(
    name = "Saint Barthelemy", iso2 = "BL",
    bbox = c(xmin = -63.82, ymin = 16.98, xmax = -61.84, ymax = 18.83),
    label_position = c(lon = -62.8332, lat = 17.902)
  ),
  BLR = list(
    name = "Belarus", iso2 = "BY",
    bbox = c(xmin = 21.51, ymin = 50.33, xmax = 34.38, ymax = 57.06),
    label_position = c(lon = 28.4177, lat = 53.8219)
  ),
  BLZ = list(
    name = "Belize", iso2 = "BZ",
    bbox = c(xmin = -90.19, ymin = 14.98, xmax = -86.83, ymax = 19.39),
    label_position = c(lon = -88.713, lat = 17.2021)
  ),
  BMU = list(
    name = "Bermuda", iso2 = "BM",
    bbox = c(xmin = -65.97, ymin = 31.34, xmax = -63.57, ymax = 33.29),
    label_position = c(lon = -64.7636, lat = 32.2966)
  ),
  BOL = list(
    name = "Bolivia", iso2 = "BO",
    bbox = c(xmin = -70.65, ymin = -23.8, xmax = -56.48, ymax = -8.78),
    label_position = c(lon = -64.5934, lat = -16.666)
  ),
  BRA = list(
    name = "Brazil", iso2 = "BR",
    bbox = c(xmin = -75.12, ymin = -34.65, xmax = -27.78, ymax = 6.17),
    label_position = c(lon = -49.5594, lat = -12.0987)
  ),
  BRB = list(
    name = "Barbados", iso2 = "BB",
    bbox = c(xmin = -60.59, ymin = 12.15, xmax = -58.5, ymax = 14.25),
    label_position = c(lon = -59.569, lat = 13.1637)
  ),
  BRN = list(
    name = "Brunei", iso2 = "BN",
    bbox = c(xmin = 113.09, ymin = 3.11, xmax = 116.27, ymax = 5.96),
    label_position = c(lon = 114.5519, lat = 4.4483)
  ),
  BTN = list(
    name = "Bhutan", iso2 = "BT",
    bbox = c(xmin = 87.7, ymin = 25.79, xmax = 93.12, ymax = 29.26),
    label_position = c(lon = 90.0403, lat = 27.5367)
  ),
  BVT = list(
    name = "Bouvet Island", iso2 = "BV",
    bbox = c(xmin = 1.76, ymin = -55.37, xmax = 5.07, ymax = -53.48),
    label_position = c(lon = 3.4131, lat = -54.4193)
  ),
  BWA = list(
    name = "Botswana", iso2 = "BW",
    bbox = c(xmin = 18.96, ymin = -27.8, xmax = 30.37, ymax = -16.88),
    label_position = c(lon = 24.1792, lat = -22.1026)
  ),
  CAF = list(
    name = "Central African Republic", iso2 = "CF",
    bbox = c(xmin = 13.46, ymin = 1.33, xmax = 28.36, ymax = 11.9),
    label_position = c(lon = 20.9069, lat = 6.9897)
  ),
  CAN = list(
    name = "Canada", iso2 = "CA",
    bbox = c(xmin = -149.63, ymin = 40.77, xmax = -44, ymax = 84.02),
    label_position = c(lon = -101.9107, lat = 60.3243)
  ),
  CHE = list(
    name = "Switzerland", iso2 = "CH",
    bbox = c(xmin = 4.59, ymin = 44.92, xmax = 11.83, ymax = 48.7),
    label_position = c(lon = 7.464, lat = 46.7191)
  ),
  CHL = list(
    name = "Chile", iso2 = "CL",
    bbox = c(xmin = -82.42, ymin = -56.82, xmax = -64.77, ymax = -16.6),
    bbox_full = c(xmin = -111.1, ymin = -56.82, xmax = -64.77, ymax = -16.6),
    label_position = c(lon = -72.3189, lat = -38.1518)
  ),
  CHN = list(
    name = "China", iso2 = "CN",
    bbox = c(xmin = 72.05, ymin = 14.87, xmax = 136.32, ymax = 54.47),
    label_position = c(lon = 106.3373, lat = 32.4982)
  ),
  CIV = list(
    name = "Ivory Coast", iso2 = "CI",
    bbox = c(xmin = -9.54, ymin = 3.44, xmax = -1.58, ymax = 11.63),
    label_position = c(lon = -5.5686, lat = 7.4914)
  ),
  CMR = list(
    name = "Cameroon", iso2 = "CM",
    bbox = c(xmin = 7.57, ymin = 0.75, xmax = 17.14, ymax = 13.98),
    label_position = c(lon = 12.4735, lat = 4.585)
  ),
  COD = list(
    name = "Democratic Republic of the Congo", iso2 = "CD",
    bbox = c(xmin = 11.28, ymin = -14.36, xmax = 32.21, ymax = 6.28),
    label_position = c(lon = 23.4588, lat = -1.8582)
  ),
  COG = list(
    name = "Republic of the Congo", iso2 = "CG",
    bbox = c(xmin = 10.21, ymin = -5.92, xmax = 19.55, ymax = 4.61),
    label_position = c(lon = 15.9005, lat = 0.1423)
  ),
  COK = list(
    name = "Cook Islands", iso2 = "CK",
    bbox = c(xmin = -160.83, ymin = -22.84, xmax = -156.33, ymax = -17.93),
    bbox_full = c(xmin = -166.8, ymin = -22.84, xmax = -156.33, ymax = -8.04),
    label_position = c(lon = -159.7857, lat = -21.216)
  ),
  COL = list(
    name = "Colombia", iso2 = "CO",
    bbox = c(xmin = -82.66, ymin = -5.14, xmax = -65.94, ymax = 14.48),
    label_position = c(lon = -73.1743, lat = 3.3731)
  ),
  COM = list(
    name = "Comoros", iso2 = "KM",
    bbox = c(xmin = 42.29, ymin = -13.28, xmax = 45.46, ymax = -10.46),
    label_position = c(lon = 43.3181, lat = -11.7277)
  ),
  CPV = list(
    name = "Cabo Verde", iso2 = "CV",
    bbox = c(xmin = -26.31, ymin = 13.9, xmax = -21.72, ymax = 18.1),
    label_position = c(lon = -23.6394, lat = 15.0748)
  ),
  CRI = list(
    name = "Costa Rica", iso2 = "CR",
    bbox = c(xmin = -88.04, ymin = 4.61, xmax = -81.64, ymax = 12.11),
    label_position = c(lon = -84.0779, lat = 10.0651)
  ),
  CUB = list(
    name = "Cuba", iso2 = "CU",
    bbox = c(xmin = -85.94, ymin = 18.92, xmax = -73.14, ymax = 24.17),
    label_position = c(lon = -77.9759, lat = 21.334)
  ),
  CUW = list(
    name = "Curaçao", iso2 = "CW",
    bbox = c(xmin = -70.1, ymin = 11.14, xmax = -67.81, ymax = 13.29),
    label_position = c(lon = -68.9206, lat = 12.145)
  ),
  CYM = list(
    name = "Cayman Islands", iso2 = "KY",
    bbox = c(xmin = -82.38, ymin = 18.36, xmax = -78.76, ymax = 20.66),
    label_position = c(lon = -81.2405, lat = 19.3199)
  ),
  CYP = list(
    name = "Cyprus", iso2 = "CY",
    bbox = c(xmin = 31.16, ymin = 33.72, xmax = 35.22, ymax = 36.09),
    label_position = c(lon = 33.0842, lat = 34.9133)
  ),
  CZE = list(
    name = "Czechia", iso2 = "CZ",
    bbox = c(xmin = 10.61, ymin = 47.65, xmax = 20.3, ymax = 51.94),
    label_position = c(lon = 15.3776, lat = 49.8824)
  ),
  DEU = list(
    name = "Germany", iso2 = "DE",
    bbox = c(xmin = 4.24, ymin = 46.37, xmax = 16.63, ymax = 55.97),
    label_position = c(lon = 9.6783, lat = 50.9617)
  ),
  DJI = list(
    name = "Djibouti", iso2 = "DJ",
    bbox = c(xmin = 40.82, ymin = 10.03, xmax = 44.35, ymax = 13.61),
    label_position = c(lon = 42.4988, lat = 11.9763)
  ),
  DMA = list(
    name = "Dominica", iso2 = "DM",
    bbox = c(xmin = -62.43, ymin = 14.3, xmax = -60.31, ymax = 16.54),
    label_position = c(lon = -61.345, lat = 15.4588)
  ),
  DNK = list(
    name = "Denmark", iso2 = "DK",
    bbox = c(xmin = 6.36, ymin = 53.67, xmax = 16.88, ymax = 58.65),
    label_position = c(lon = 9.0182, lat = 55.967)
  ),
  DOM = list(
    name = "Dominican Republic", iso2 = "DO",
    bbox = c(xmin = -72.98, ymin = 16.64, xmax = -67.36, ymax = 20.84),
    label_position = c(lon = -70.654, lat = 19.1041)
  ),
  DZA = list(
    name = "Algeria", iso2 = "DZ",
    bbox = c(xmin = -9.83, ymin = 18.07, xmax = 13.11, ymax = 38),
    label_position = c(lon = 2.8082, lat = 27.3974)
  ),
  ECU = list(
    name = "Ecuador", iso2 = "EC",
    bbox = c(xmin = -81.92, ymin = -5.91, xmax = -74.32, ymax = 2.34),
    bbox_full = c(xmin = -92.92, ymin = -5.91, xmax = -74.32, ymax = 2.57),
    label_position = c(lon = -78.1884, lat = -1.2591)
  ),
  EGY = list(
    name = "Egypt", iso2 = "EG",
    bbox = c(xmin = 23.62, ymin = 21.09, xmax = 37.97, ymax = 32.56),
    label_position = c(lon = 29.4458, lat = 26.1862)
  ),
  ERI = list(
    name = "Eritrea", iso2 = "ER",
    bbox = c(xmin = 35.47, ymin = 11.46, xmax = 44.08, ymax = 18.91),
    label_position = c(lon = 38.2856, lat = 15.7874)
  ),
  ESH = list(
    name = "Western Sahara", iso2 = "EH",
    bbox = c(xmin = -18.13, ymin = 19.86, xmax = -7.65, ymax = 28.56),
    label_position = c(lon = -12.6303, lat = 23.9676)
  ),
  ESP = list(
    name = "Spain", iso2 = "ES",
    bbox = c(xmin = -10.56, ymin = 34.27, xmax = 5.61, ymax = 44.7),
    bbox_full = c(xmin = -19.44, ymin = 26.74, xmax = 5.61, ymax = 44.7),
    label_position = c(lon = -3.4647, lat = 40.091)
  ),
  EST = list(
    name = "Estonia", iso2 = "EE",
    bbox = c(xmin = 20, ymin = 56.61, xmax = 30.02, ymax = 60.57),
    label_position = c(lon = 25.8671, lat = 58.7249)
  ),
  ETH = list(
    name = "Ethiopia", iso2 = "ET",
    bbox = c(xmin = 32.05, ymin = 2.5, xmax = 48.92, ymax = 15.78),
    label_position = c(lon = 39.0886, lat = 8.0328)
  ),
  FIN = list(
    name = "Finland", iso2 = "FI",
    bbox = c(xmin = 17.86, ymin = 58.91, xmax = 34.33, ymax = 70.98),
    label_position = c(lon = 27.2764, lat = 63.2524)
  ),
  FJI = list(
    name = "Fiji", iso2 = "FJ",
    bbox = c(xmin = 173.61, ymin = -22.61, xmax = 182.76, ymax = -11.57),
    label_position = c(lon = 177.9754, lat = -17.8261)
  ),
  FLK = list(
    name = "Falkland Islands", iso2 = "FK",
    bbox = c(xmin = -62.83, ymin = -53.31, xmax = -56.23, ymax = -50.12),
    label_position = c(lon = -58.7386, lat = -51.6089)
  ),
  FRA = list(
    name = "France", iso2 = NA_character_,
    bbox = c(xmin = -6.6, ymin = 40.46, xmax = 11.02, ymax = 51.99),
    label_position = c(lon = 2.5523, lat = 46.6961)
  ),
  FRO = list(
    name = "Faroe Islands", iso2 = "FO",
    bbox = c(xmin = -9.65, ymin = 60.49, xmax = -4.27, ymax = 63.3),
    label_position = c(lon = -7.0584, lat = 62.1856)
  ),
  FSM = list(
    name = "Federated States of Micronesia", iso2 = "FM",
    bbox = c(xmin = 137.14, ymin = 0.01, xmax = 163.97, ymax = 10.68),
    label_position = c(lon = 158.234, lat = 6.8876)
  ),
  GAB = list(
    name = "Gabon", iso2 = "GA",
    bbox = c(xmin = 7.79, ymin = -4.84, xmax = 15.41, ymax = 3.23),
    label_position = c(lon = 11.8359, lat = -0.4377)
  ),
  GBR = list(
    name = "United Kingdom", iso2 = "GB",
    bbox = c(xmin = -10.52, ymin = 49.01, xmax = 3.67, ymax = 61.75),
    bbox_full = c(xmin = -15.59, ymin = 49.01, xmax = 3.67, ymax = 61.75),
    label_position = c(lon = -2.1163, lat = 54.4027)
  ),
  GEO = list(
    name = "Georgia", iso2 = "GE",
    bbox = c(xmin = 38.72, ymin = 40.14, xmax = 47.96, ymax = 44.48),
    label_position = c(lon = 43.7357, lat = 41.8701)
  ),
  GGY = list(
    name = "Guernsey", iso2 = "GG",
    bbox = c(xmin = -4.09, ymin = 48.51, xmax = -0.75, ymax = 50.63),
    label_position = c(lon = -2.5617, lat = 49.4635)
  ),
  GHA = list(
    name = "Ghana", iso2 = "GH",
    bbox = c(xmin = -4.19, ymin = 3.83, xmax = 2.11, ymax = 12.07),
    label_position = c(lon = -1.0369, lat = 7.7176)
  ),
  GIB = list(
    name = "Gibraltar", iso2 = "GI",
    bbox = c(xmin = -6.49, ymin = 35.21, xmax = -4.21, ymax = 37.04),
    label_position = c(lon = -5.3467, lat = 36.1294)
  ),
  GIN = list(
    name = "Guinea", iso2 = "GN",
    bbox = c(xmin = -16.01, ymin = 6.29, xmax = -6.73, ymax = 13.58),
    label_position = c(lon = -10.0164, lat = 10.6185)
  ),
  GLP = list(
    name = "Guadeloupe", iso2 = "GP",
    bbox = c(xmin = -62.74, ymin = 14.94, xmax = -60.04, ymax = 17.42),
    label_position = c(lon = -61.6745, lat = 16.1669)
  ),
  GMB = list(
    name = "Gambia", iso2 = "GM",
    bbox = c(xmin = -17.76, ymin = 12.16, xmax = -12.88, ymax = 14.72),
    label_position = c(lon = -14.9983, lat = 13.6417)
  ),
  GNB = list(
    name = "Guinea-Bissau", iso2 = "GW",
    bbox = c(xmin = -17.66, ymin = 10.02, xmax = -12.73, ymax = 13.58),
    label_position = c(lon = -14.5241, lat = 12.1637)
  ),
  GNQ = list(
    name = "Equatorial Guinea", iso2 = "GQ",
    bbox = c(xmin = 4.71, ymin = -2.38, xmax = 12.24, ymax = 4.68),
    label_position = c(lon = 8.9902, lat = 2.333)
  ),
  GRC = list(
    name = "Greece", iso2 = "GR",
    bbox = c(xmin = 18.4, ymin = 33.91, xmax = 29.47, ymax = 42.65),
    label_position = c(lon = 21.7257, lat = 39.4928)
  ),
  GRD = list(
    name = "Grenada", iso2 = "GD",
    bbox = c(xmin = -62.72, ymin = 11.1, xmax = -60.49, ymax = 13.43),
    label_position = c(lon = -61.6805, lat = 12.1132)
  ),
  GRL = list(
    name = "Greenland", iso2 = "GL",
    bbox = c(xmin = -82.49, ymin = 58.89, xmax = -1.94, ymax = 84.54),
    label_position = c(lon = -39.3353, lat = 74.3194)
  ),
  GTM = list(
    name = "Guatemala", iso2 = "GT",
    bbox = c(xmin = -93.2, ymin = 12.83, xmax = -87.27, ymax = 18.72),
    label_position = c(lon = -90.4971, lat = 14.9821)
  ),
  GUF = list(
    name = "French Guiana", iso2 = "GF",
    bbox = c(xmin = -55.52, ymin = 1.21, xmax = -50.74, ymax = 6.65),
    label_position = c(lon = -53.2446, lat = 3.9225)
  ),
  GUM = list(
    name = "Guam", iso2 = "GU",
    bbox = c(xmin = 143.69, ymin = 12.34, xmax = 145.89, ymax = 14.56),
    label_position = c(lon = 144.7036, lat = 13.3542)
  ),
  GUY = list(
    name = "Guyana", iso2 = "GY",
    bbox = c(xmin = -62.31, ymin = 0.28, xmax = -55.57, ymax = 9.46),
    label_position = c(lon = -58.9426, lat = 5.1243)
  ),
  HKG = list(
    name = "Hong Kong S.A.R.", iso2 = "HK",
    bbox = c(xmin = 112.85, ymin = 21.27, xmax = 115.39, ymax = 23.47),
    label_position = c(lon = 114.0978, lat = 22.4488)
  ),
  HMD = list(
    name = "Heard Island and McDonald Islands", iso2 = "HM",
    bbox = c(xmin = 71.7, ymin = -54.1, xmax = 75.35, ymax = -52.06),
    label_position = c(lon = 73.5052, lat = -53.1035)
  ),
  HND = list(
    name = "Honduras", iso2 = "HN",
    bbox = c(xmin = -90.32, ymin = 12.08, xmax = -82.18, ymax = 18.32),
    label_position = c(lon = -86.8876, lat = 14.7948)
  ),
  HRV = list(
    name = "Croatia", iso2 = "HR",
    bbox = c(xmin = 12.17, ymin = 41.51, xmax = 20.74, ymax = 47.45),
    label_position = c(lon = 16.3724, lat = 45.8058)
  ),
  HTI = list(
    name = "Haiti", iso2 = "HT",
    bbox = c(xmin = -75.46, ymin = 17.12, xmax = -70.67, ymax = 20.99),
    label_position = c(lon = -72.2241, lat = 19.2638)
  ),
  HUN = list(
    name = "Hungary", iso2 = "HU",
    bbox = c(xmin = 14.71, ymin = 44.84, xmax = 24.26, ymax = 49.47),
    label_position = c(lon = 19.4479, lat = 47.0868)
  ),
  IDN = list(
    name = "Indonesia", iso2 = "ID",
    bbox = c(xmin = 94.09, ymin = -11.83, xmax = 141.9, ymax = 6.81),
    label_position = c(lon = 101.8929, lat = -0.9544)
  ),
  IMN = list(
    name = "Isle of Man", iso2 = "IM",
    bbox = c(xmin = -6.37, ymin = 53.15, xmax = -2.73, ymax = 55.32),
    label_position = c(lon = -4.5301, lat = 54.2208)
  ),
  IND = list(
    name = "India", iso2 = "IN",
    bbox = c(xmin = 67.02, ymin = 5.84, xmax = 98.48, ymax = 36.4),
    label_position = c(lon = 79.3581, lat = 22.6869)
  ),
  IOT = list(
    name = "British Indian Ocean Territory", iso2 = "IO",
    bbox = c(xmin = 70.35, ymin = -8.34, xmax = 73.41, ymax = -4.32),
    label_position = c(lon = 71.3483, lat = -6.1908)
  ),
  IRL = list(
    name = "Ireland", iso2 = "IE",
    bbox = c(xmin = -12.1, ymin = 50.54, xmax = -4.37, ymax = 56.29),
    label_position = c(lon = -7.7986, lat = 53.0787)
  ),
  IRN = list(
    name = "Iran", iso2 = "IR",
    bbox = c(xmin = 42.83, ymin = 24.16, xmax = 64.51, ymax = 40.67),
    label_position = c(lon = 54.9315, lat = 32.1662)
  ),
  IRQ = list(
    name = "Iraq", iso2 = "IQ",
    bbox = c(xmin = 37.63, ymin = 28.16, xmax = 49.71, ymax = 38.28),
    label_position = c(lon = 43.2618, lat = 33.094)
  ),
  ISL = list(
    name = "Iceland", iso2 = "IS",
    bbox = c(xmin = -26.89, ymin = 62.49, xmax = -11.15, ymax = 67.47),
    label_position = c(lon = -18.6737, lat = 64.7793)
  ),
  ISR = list(
    name = "Israel", iso2 = "IL",
    bbox = c(xmin = 33.16, ymin = 28.59, xmax = 36.98, ymax = 34.31),
    label_position = c(lon = 34.8479, lat = 30.9111)
  ),
  ITA = list(
    name = "Italy", iso2 = "IT",
    bbox = c(xmin = 5.26, ymin = 34.59, xmax = 19.86, ymax = 47.99),
    label_position = c(lon = 11.0769, lat = 44.7325)
  ),
  JAM = list(
    name = "Jamaica", iso2 = "JM",
    bbox = c(xmin = -79.33, ymin = 16.8, xmax = -75.23, ymax = 19.43),
    label_position = c(lon = -77.3188, lat = 18.1371)
  ),
  JEY = list(
    name = "Jersey", iso2 = "JE",
    bbox = c(xmin = -3.65, ymin = 48.27, xmax = -0.6, ymax = 50.17),
    label_position = c(lon = -2.0901, lat = 49.2208)
  ),
  JOR = list(
    name = "Jordan", iso2 = "JO",
    bbox = c(xmin = 33.86, ymin = 28.29, xmax = 40.38, ymax = 34.27),
    label_position = c(lon = 36.376, lat = 30.805)
  ),
  JPN = list(
    name = "Japan", iso2 = "JP",
    bbox = c(xmin = 121.63, ymin = 23.31, xmax = 147.13, ymax = 46.42),
    bbox_full = c(xmin = 121.63, ymin = 23.31, xmax = 155.29, ymax = 46.42),
    label_position = c(lon = 138.4422, lat = 36.1425)
  ),
  KAZ = list(
    name = "Kazakhstan", iso2 = "KZ",
    bbox = c(xmin = 44.85, ymin = 39.68, xmax = 88.95, ymax = 56.34),
    label_position = c(lon = 68.6855, lat = 49.0541)
  ),
  KEN = list(
    name = "Kenya", iso2 = "KE",
    bbox = c(xmin = 32.98, ymin = -5.58, xmax = 42.79, ymax = 5.93),
    label_position = c(lon = 37.9076, lat = 0.549)
  ),
  KGZ = list(
    name = "Kyrgyzstan", iso2 = "KG",
    bbox = c(xmin = 67.97, ymin = 38.29, xmax = 81.51, ymax = 44.17),
    label_position = c(lon = 74.5326, lat = 41.6685)
  ),
  KHM = list(
    name = "Cambodia", iso2 = "KH",
    bbox = c(xmin = 101.38, ymin = 9.51, xmax = 108.55, ymax = 15.61),
    label_position = c(lon = 104.5049, lat = 12.6476)
  ),
  KIR = list(
    name = "Kiribati", iso2 = "KI",
    bbox = c(xmin = -161.32, ymin = 0.81, xmax = -156.27, ymax = 5.63),
    bbox_full = c(xmin = 168.6, ymin = -12.36, xmax = 209.14, ymax = 5.63),
    label_position = c(lon = -157.3846, lat = 1.8204)
  ),
  KNA = list(
    name = "Saint Kitts and Nevis", iso2 = "KN",
    bbox = c(xmin = -63.81, ymin = 16.2, xmax = -61.59, ymax = 18.32),
    label_position = c(lon = -62.758, lat = 17.3366)
  ),
  KOR = list(
    name = "South Korea", iso2 = "KR",
    bbox = c(xmin = 123.44, ymin = 32.29, xmax = 133.03, ymax = 39.53),
    label_position = c(lon = 128.1295, lat = 36.3849)
  ),
  KWT = list(
    name = "Kuwait", iso2 = "KW",
    bbox = c(xmin = 45.48, ymin = 27.63, xmax = 49.49, ymax = 31),
    label_position = c(lon = 47.314, lat = 29.4136)
  ),
  LAO = list(
    name = "Laos", iso2 = "LA",
    bbox = c(xmin = 99.11, ymin = 13.01, xmax = 108.65, ymax = 23.4),
    label_position = c(lon = 102.5339, lat = 19.4318)
  ),
  LBN = list(
    name = "Lebanon", iso2 = "LB",
    bbox = c(xmin = 33.99, ymin = 32.15, xmax = 37.71, ymax = 35.59),
    label_position = c(lon = 35.9929, lat = 34.1334)
  ),
  LBR = list(
    name = "Liberia", iso2 = "LR",
    bbox = c(xmin = -12.39, ymin = 3.44, xmax = -6.47, ymax = 9.47),
    label_position = c(lon = -9.4604, lat = 6.4472)
  ),
  LBY = list(
    name = "Libya", iso2 = "LY",
    bbox = c(xmin = 8.2, ymin = 18.59, xmax = 26.25, ymax = 34.08),
    label_position = c(lon = 18.011, lat = 26.6389)
  ),
  LCA = list(
    name = "Saint Lucia", iso2 = "LC",
    bbox = c(xmin = -62.01, ymin = 12.81, xmax = -59.95, ymax = 15.02),
    label_position = c(lon = -60.9801, lat = 13.8924)
  ),
  LIE = list(
    name = "Liechtenstein", iso2 = "LI",
    bbox = c(xmin = 8.12, ymin = 46.15, xmax = 10.97, ymax = 48.17),
    label_position = c(lon = 9.5594, lat = 47.1114)
  ),
  LKA = list(
    name = "Sri Lanka", iso2 = "LK",
    bbox = c(xmin = 78.74, ymin = 5.02, xmax = 82.81, ymax = 10.73),
    label_position = c(lon = 80.7048, lat = 7.5811)
  ),
  LSO = list(
    name = "Lesotho", iso2 = "LS",
    bbox = c(xmin = 25.94, ymin = -31.56, xmax = 30.5, ymax = -27.67),
    label_position = c(lon = 28.2466, lat = -29.4802)
  ),
  LTU = list(
    name = "Lithuania", iso2 = "LT",
    bbox = c(xmin = 19.25, ymin = 52.98, xmax = 28.47, ymax = 57.35),
    label_position = c(lon = 24.0899, lat = 55.1037)
  ),
  LUX = list(
    name = "Luxembourg", iso2 = "LU",
    bbox = c(xmin = 4.28, ymin = 48.54, xmax = 7.94, ymax = 51.08),
    label_position = c(lon = 6.0776, lat = 49.7337)
  ),
  LVA = list(
    name = "Latvia", iso2 = "LV",
    bbox = c(xmin = 19.22, ymin = 54.76, xmax = 29.97, ymax = 58.98),
    label_position = c(lon = 25.4587, lat = 57.0669)
  ),
  MAC = list(
    name = "Macao S.A.R", iso2 = "MO",
    bbox = c(xmin = 112.54, ymin = 21.2, xmax = 114.57, ymax = 23.12),
    label_position = c(lon = 113.556, lat = 22.1297)
  ),
  MAF = list(
    name = "Saint Martin", iso2 = "MF",
    bbox = c(xmin = -64.1, ymin = 17.13, xmax = -62.06, ymax = 19.03),
    label_position = c(lon = -63.0494, lat = 18.0813)
  ),
  MAR = list(
    name = "Morocco", iso2 = "MA",
    bbox = c(xmin = -18.14, ymin = 20.52, xmax = 0.1, ymax = 36.83),
    label_position = c(lon = -7.1873, lat = 31.6507)
  ),
  MCO = list(
    name = "Monaco", iso2 = "MC",
    bbox = c(xmin = 6.1, ymin = 42.81, xmax = 8.71, ymax = 44.67),
    label_position = c(lon = 7.3983, lat = 43.7397)
  ),
  MDA = list(
    name = "Moldova", iso2 = "MD",
    bbox = c(xmin = 25.23, ymin = 44.56, xmax = 31.52, ymax = 49.39),
    label_position = c(lon = 28.4879, lat = 47.435)
  ),
  MDG = list(
    name = "Madagascar", iso2 = "MG",
    bbox = c(xmin = 42.21, ymin = -26.5, xmax = 51.51, ymax = -11.04),
    label_position = c(lon = 46.7042, lat = -18.6283)
  ),
  MDV = list(
    name = "Maldives", iso2 = "MV",
    bbox = c(xmin = 71.77, ymin = -1.59, xmax = 74.67, ymax = 8.01),
    label_position = c(lon = 73.5076, lat = 4.1744)
  ),
  MEX = list(
    name = "Mexico", iso2 = "MX",
    bbox = c(xmin = -119.45, ymin = 13.64, xmax = -85.62, ymax = 33.62),
    label_position = c(lon = -102.2894, lat = 23.92)
  ),
  MHL = list(
    name = "Marshall Islands", iso2 = "MH",
    bbox = c(xmin = 164.34, ymin = 3.67, xmax = 172.97, ymax = 15.51),
    label_position = c(lon = 171.1936, lat = 7.0826)
  ),
  MKD = list(
    name = "North Macedonia", iso2 = "MK",
    bbox = c(xmin = 19.21, ymin = 39.95, xmax = 24.25, ymax = 43.27),
    label_position = c(lon = 21.5558, lat = 41.5582)
  ),
  MLI = list(
    name = "Mali", iso2 = "ML",
    bbox = c(xmin = -13.27, ymin = 9.24, xmax = 5.24, ymax = 25.9),
    label_position = c(lon = -2.0385, lat = 18.6927)
  ),
  MLT = list(
    name = "Malta", iso2 = "MT",
    bbox = c(xmin = 13.05, ymin = 34.9, xmax = 15.7, ymax = 36.98),
    label_position = c(lon = 14.433, lat = 35.8929)
  ),
  MMR = list(
    name = "Myanmar", iso2 = "MM",
    bbox = c(xmin = 91.14, ymin = 8.89, xmax = 102.21, ymax = 29.44),
    label_position = c(lon = 95.8045, lat = 21.5739)
  ),
  MNE = list(
    name = "Montenegro", iso2 = "ME",
    bbox = c(xmin = 17.17, ymin = 40.95, xmax = 21.62, ymax = 44.45),
    label_position = c(lon = 19.1437, lat = 42.8031)
  ),
  MNG = list(
    name = "Mongolia", iso2 = "MN",
    bbox = c(xmin = 86.24, ymin = 40.68, xmax = 121.41, ymax = 53.03),
    label_position = c(lon = 104.1504, lat = 45.9975)
  ),
  MNP = list(
    name = "Northern Mariana Islands", iso2 = "MP",
    bbox = c(xmin = 143.93, ymin = 13.21, xmax = 146.84, ymax = 21.46),
    label_position = c(lon = 145.7344, lat = 15.1882)
  ),
  MOZ = list(
    name = "Mozambique", iso2 = "MZ",
    bbox = c(xmin = 29.19, ymin = -27.76, xmax = 41.87, ymax = -9.57),
    label_position = c(lon = 37.8379, lat = -13.9432)
  ),
  MRT = list(
    name = "Mauritania", iso2 = "MR",
    bbox = c(xmin = -18.11, ymin = 13.83, xmax = -3.8, ymax = 28.19),
    label_position = c(lon = -9.7403, lat = 19.5871)
  ),
  MSR = list(
    name = "Montserrat", iso2 = "MS",
    bbox = c(xmin = -63.18, ymin = 15.77, xmax = -61.19, ymax = 17.72),
    label_position = c(lon = -62.1883, lat = 16.7372)
  ),
  MTQ = list(
    name = "Martinique", iso2 = "MQ",
    bbox = c(xmin = -62.17, ymin = 13.5, xmax = -59.87, ymax = 15.78),
    label_position = c(lon = -61.0176, lat = 14.6564)
  ),
  MUS = list(
    name = "Mauritius", iso2 = "MU",
    bbox = c(xmin = 56.33, ymin = -21.42, xmax = 58.77, ymax = -19.08),
    bbox_full = c(xmin = 55.55, ymin = -21.42, xmax = 64.46, ymax = -9.42),
    label_position = c(lon = 57.5658, lat = -20.2995)
  ),
  MWI = list(
    name = "Malawi", iso2 = "MW",
    bbox = c(xmin = 31.71, ymin = -18.04, xmax = 36.85, ymax = -8.48),
    label_position = c(lon = 33.6081, lat = -13.3867)
  ),
  MYS = list(
    name = "Malaysia", iso2 = "MY",
    bbox = c(xmin = 98.73, ymin = -0.05, xmax = 120.19, ymax = 8.26),
    label_position = c(lon = 113.8371, lat = 2.5287)
  ),
  MYT = list(
    name = "Mayotte", iso2 = "YT",
    bbox = c(xmin = 44.11, ymin = -13.89, xmax = 46.22, ymax = -11.74),
    label_position = c(lon = 45.1375, lat = -12.8189)
  ),
  NAM = list(
    name = "Namibia", iso2 = "NA",
    bbox = c(xmin = 10.68, ymin = -29.86, xmax = 26.3, ymax = -16.05),
    label_position = c(lon = 17.1082, lat = -20.5753)
  ),
  NCL = list(
    name = "New Caledonia", iso2 = "NC",
    bbox = c(xmin = 162.63, ymin = -23.57, xmax = 172.33, ymax = -18.72),
    label_position = c(lon = 165.084, lat = -21.0647)
  ),
  NER = list(
    name = "Niger", iso2 = "NE",
    bbox = c(xmin = -0.84, ymin = 10.79, xmax = 16.96, ymax = 24.42),
    label_position = c(lon = 9.5044, lat = 17.4462)
  ),
  NFK = list(
    name = "Norfolk Island", iso2 = "NF",
    bbox = c(xmin = 166.87, ymin = -29.98, xmax = 169.04, ymax = -28.09),
    label_position = c(lon = 167.9545, lat = -29.033)
  ),
  NGA = list(
    name = "Nigeria", iso2 = "NG",
    bbox = c(xmin = 1.74, ymin = 3.37, xmax = 15.6, ymax = 14.78),
    label_position = c(lon = 7.5032, lat = 9.4398)
  ),
  NIC = list(
    name = "Nicaragua", iso2 = "NI",
    bbox = c(xmin = -88.63, ymin = 9.81, xmax = -81.79, ymax = 15.93),
    label_position = c(lon = -85.0693, lat = 12.6707)
  ),
  NIU = list(
    name = "Niue", iso2 = "NU",
    bbox = c(xmin = -170.91, ymin = -20.05, xmax = -168.82, ymax = -18.06),
    label_position = c(lon = -169.8626, lat = -19.046)
  ),
  NLD = list(
    name = "Netherlands", iso2 = "NL",
    bbox = c(xmin = 1.8, ymin = 49.84, xmax = 8.75, ymax = 54.46),
    label_position = c(lon = 5.6114, lat = 52.4222)
  ),
  NOR = list(
    name = "Norway", iso2 = NA_character_,
    bbox = c(xmin = 1.72, ymin = 57.09, xmax = 34, ymax = 72.07),
    label_position = c(lon = 9.68, lat = 61.3571)
  ),
  NPL = list(
    name = "Nepal", iso2 = "NP",
    bbox = c(xmin = 78.97, ymin = 25.44, xmax = 89.23, ymax = 31.32),
    label_position = c(lon = 83.6399, lat = 28.2979)
  ),
  NRU = list(
    name = "Nauru", iso2 = "NR",
    bbox = c(xmin = 166, ymin = -1.46, xmax = 167.86, ymax = 0.41),
    label_position = c(lon = 166.9326, lat = -0.5203)
  ),
  NZL = list(
    name = "New Zealand", iso2 = "NZ",
    bbox = c(xmin = 164.37, ymin = -53.5, xmax = 180.36, ymax = -33.24),
    bbox_full = c(xmin = 164.37, ymin = -53.5, xmax = 190.33, ymax = -7.64),
    label_position = c(lon = 172.787, lat = -39.759)
  ),
  OMN = list(
    name = "Oman", iso2 = "OM",
    bbox = c(xmin = 50.96, ymin = 15.74, xmax = 60.86, ymax = 27.29),
    label_position = c(lon = 57.3366, lat = 22.1204)
  ),
  PAK = list(
    name = "Pakistan", iso2 = "PK",
    bbox = c(xmin = 59.7, ymin = 22.79, xmax = 78.19, ymax = 37.96),
    label_position = c(lon = 68.5456, lat = 29.3284)
  ),
  PAN = list(
    name = "Panama", iso2 = "PA",
    bbox = c(xmin = -83.97, ymin = 6.3, xmax = -76.24, ymax = 10.53),
    label_position = c(lon = -80.3521, lat = 8.722)
  ),
  PCN = list(
    name = "Pitcairn Islands", iso2 = "PN",
    bbox = c(xmin = -131.76, ymin = -25.98, xmax = -123.77, ymax = -23.02),
    label_position = c(lon = -128.3175, lat = -24.3646)
  ),
  PER = list(
    name = "Peru", iso2 = "PE",
    bbox = c(xmin = -82.29, ymin = -19.24, xmax = -67.73, ymax = 0.87),
    label_position = c(lon = -72.9002, lat = -12.9767)
  ),
  PHL = list(
    name = "Philippines", iso2 = "PH",
    bbox = c(xmin = 115.98, ymin = 3.75, xmax = 127.59, ymax = 22.03),
    label_position = c(lon = 122.465, lat = 11.198)
  ),
  PLW = list(
    name = "Palau", iso2 = "PW",
    bbox = c(xmin = 130.22, ymin = 2.05, xmax = 135.64, ymax = 9),
    label_position = c(lon = 134.5802, lat = 7.5183)
  ),
  PNG = list(
    name = "Papua New Guinea", iso2 = "PG",
    bbox = c(xmin = 139.92, ymin = -12.54, xmax = 156.89, ymax = -0.44),
    label_position = c(lon = 143.9102, lat = -5.6953)
  ),
  POL = list(
    name = "Poland", iso2 = "PL",
    bbox = c(xmin = 12.52, ymin = 48.09, xmax = 25.74, ymax = 55.74),
    label_position = c(lon = 19.4905, lat = 51.9903)
  ),
  PRI = list(
    name = "Puerto Rico", iso2 = "PR",
    bbox = c(xmin = -68.9, ymin = 17.02, xmax = -64.29, ymax = 19.43),
    label_position = c(lon = -66.4811, lat = 18.2347)
  ),
  PRK = list(
    name = "North Korea", iso2 = "KP",
    bbox = c(xmin = 122.96, ymin = 36.77, xmax = 131.95, ymax = 43.91),
    label_position = c(lon = 126.4445, lat = 39.8853)
  ),
  PRT = list(
    name = "Portugal", iso2 = "PT",
    bbox = c(xmin = -10.73, ymin = 36.06, xmax = -4.97, ymax = 43.06),
    bbox_full = c(xmin = -32.52, ymin = 29.13, xmax = -4.97, ymax = 43.06),
    label_position = c(lon = -8.2718, lat = 39.6067)
  ),
  PRY = list(
    name = "Paraguay", iso2 = "PY",
    bbox = c(xmin = -63.68, ymin = -28.49, xmax = -53.22, ymax = -18.38),
    label_position = c(lon = -60.1464, lat = -21.6745)
  ),
  PSE = list(
    name = "Palestine", iso2 = "PS",
    bbox = c(xmin = 33.12, ymin = 30.31, xmax = 36.65, ymax = 33.45),
    label_position = c(lon = 35.2913, lat = 32.0474)
  ),
  PYF = list(
    name = "French Polynesia", iso2 = "PF",
    bbox = c(xmin = -155.56, ymin = -28.54, xmax = -133.92, ymax = -7.05),
    label_position = c(lon = -149.4616, lat = -17.6281)
  ),
  QAT = list(
    name = "Qatar", iso2 = "QA",
    bbox = c(xmin = 49.74, ymin = 23.66, xmax = 52.63, ymax = 27.06),
    label_position = c(lon = 51.1435, lat = 25.2374)
  ),
  REU = list(
    name = "Reunion", iso2 = "RE",
    bbox = c(xmin = 54.25, ymin = -22.27, xmax = 56.83, ymax = -19.96),
    label_position = c(lon = 55.5417, lat = -21.1195)
  ),
  ROU = list(
    name = "Romania", iso2 = "RO",
    bbox = c(xmin = 18.86, ymin = 42.75, xmax = 31.08, ymax = 49.18),
    label_position = c(lon = 24.9726, lat = 45.7332)
  ),
  RUS = list(
    name = "Russia", iso2 = "RU",
    bbox = c(xmin = 19.81, ymin = 40.29, xmax = 187.13, ymax = 82.76),
    bbox_full = c(xmin = 12.48, ymin = 40.29, xmax = 198.14, ymax = 82.76),
    label_position = c(lon = 44.6865, lat = 58.2494)
  ),
  RWA = list(
    name = "Rwanda", iso2 = "RW",
    bbox = c(xmin = 27.95, ymin = -3.73, xmax = 31.79, ymax = -0.16),
    label_position = c(lon = 30.1039, lat = -1.8972)
  ),
  SAU = list(
    name = "Saudi Arabia", iso2 = "SA",
    bbox = c(xmin = 33.5, ymin = 15.47, xmax = 56.71, ymax = 33.02),
    label_position = c(lon = 44.6996, lat = 23.8069)
  ),
  SDN = list(
    name = "Sudan", iso2 = "SD",
    bbox = c(xmin = 20.83, ymin = 7.78, xmax = 39.59, ymax = 23.13),
    label_position = c(lon = 29.2607, lat = 16.3307)
  ),
  SEN = list(
    name = "Senegal", iso2 = "SN",
    bbox = c(xmin = -18.48, ymin = 11.4, xmax = -10.43, ymax = 17.59),
    label_position = c(lon = -14.7786, lat = 15.1381)
  ),
  SGP = list(
    name = "Singapore", iso2 = "SG",
    bbox = c(xmin = 102.74, ymin = 0.36, xmax = 104.91, ymax = 2.35),
    label_position = c(lon = 103.8169, lat = 1.3666)
  ),
  SGS = list(
    name = "South Georgia and the Islands", iso2 = "GS",
    bbox = c(xmin = -39.69, ymin = -55.79, xmax = -34.18, ymax = -53.07),
    bbox_full = c(xmin = -39.91, ymin = -60.38, xmax = -24.42, ymax = -53.07),
    label_position = c(lon = -36.6837, lat = -54.3677)
  ),
  SHN = list(
    name = "Saint Helena", iso2 = "SH",
    bbox = c(xmin = -15.36, ymin = -16.92, xmax = -4.71, ymax = -6.97),
    bbox_full = c(xmin = -15.62, ymin = -41.3, xmax = -4.45, ymax = -6.97),
    label_position = c(lon = -5.7126, lat = -15.9505)
  ),
  SJM = list(
    name = "Svalbard and Jan Mayen", iso2 = "SJ",
    bbox = c(xmin = 4.28, ymin = 73.44, xmax = 39.84, ymax = 81.67),
    bbox_full = c(xmin = -15.32, ymin = 69.9, xmax = 39.84, ymax = 81.67),
    label_position = c(lon = 15.9011, lat = 78.61)
  ),
  SLB = list(
    name = "Solomon Islands", iso2 = "SB",
    bbox = c(xmin = 154.58, ymin = -13.19, xmax = 169.75, ymax = -5.7),
    label_position = c(lon = 159.1705, lat = -8.0295)
  ),
  SLE = list(
    name = "Sierra Leone", iso2 = "SL",
    bbox = c(xmin = -14.22, ymin = 6.02, xmax = -9.36, ymax = 10.9),
    label_position = c(lon = -11.7637, lat = 8.6174)
  ),
  SLV = list(
    name = "El Salvador", iso2 = "SV",
    bbox = c(xmin = -91.05, ymin = 12.26, xmax = -86.76, ymax = 15.35),
    label_position = c(lon = -88.8901, lat = 13.6854)
  ),
  SMR = list(
    name = "San Marino", iso2 = "SM",
    bbox = c(xmin = 11.11, ymin = 42.99, xmax = 13.77, ymax = 44.89),
    label_position = c(lon = 12.4412, lat = 43.9339)
  ),
  SOM = list(
    name = "Somalia", iso2 = "SO",
    bbox = c(xmin = 40.04, ymin = -2.6, xmax = 52.34, ymax = 12.89),
    label_position = c(lon = 45.1924, lat = 3.5689)
  ),
  SPM = list(
    name = "Saint Pierre and Miquelon", iso2 = "PM",
    bbox = c(xmin = -57.75, ymin = 45.85, xmax = -54.8, ymax = 48.04),
    label_position = c(lon = -56.3324, lat = 47.0403)
  ),
  SRB = list(
    name = "Republic of Serbia", iso2 = "RS",
    bbox = c(xmin = 17.52, ymin = 41.33, xmax = 24.31, ymax = 47.08),
    label_position = c(lon = 20.788, lat = 44.1899)
  ),
  SSD = list(
    name = "South Sudan", iso2 = "SS",
    bbox = c(xmin = 23.19, ymin = 2.59, xmax = 36.85, ymax = 13.12),
    label_position = c(lon = 30.3902, lat = 7.2305)
  ),
  STP = list(
    name = "São Tomé and Principe", iso2 = "ST",
    bbox = c(xmin = 5.56, ymin = -0.88, xmax = 8.37, ymax = 2.6),
    label_position = c(lon = 7.021, lat = 0.9709)
  ),
  SUR = list(
    name = "Suriname", iso2 = "SR",
    bbox = c(xmin = -58.98, ymin = 0.93, xmax = -53.08, ymax = 6.91),
    label_position = c(lon = -55.9109, lat = 4.144)
  ),
  SVK = list(
    name = "Slovakia", iso2 = "SK",
    bbox = c(xmin = 15.43, ymin = 46.85, xmax = 23.96, ymax = 50.51),
    label_position = c(lon = 19.0499, lat = 48.734)
  ),
  SVN = list(
    name = "Slovenia", iso2 = "SI",
    bbox = c(xmin = 12.02, ymin = 44.52, xmax = 17.86, ymax = 47.77),
    label_position = c(lon = 14.9153, lat = 46.0608)
  ),
  SWE = list(
    name = "Sweden", iso2 = "SE",
    bbox = c(xmin = 8.48, ymin = 54.44, xmax = 26.79, ymax = 69.94),
    label_position = c(lon = 19.0171, lat = 65.8592)
  ),
  SWZ = list(
    name = "eSwatini", iso2 = "SZ",
    bbox = c(xmin = 29.76, ymin = -28.22, xmax = 33.14, ymax = -24.83),
    label_position = c(lon = 31.4673, lat = -26.5337)
  ),
  SXM = list(
    name = "Sint Maarten", iso2 = "SX",
    bbox = c(xmin = -64.07, ymin = 17.12, xmax = -62.06, ymax = 18.97),
    label_position = c(lon = -63.0701, lat = 18.0409)
  ),
  SYC = list(
    name = "Seychelles", iso2 = "SC",
    bbox = c(xmin = 45.29, ymin = -10.66, xmax = 57.21, ymax = -2.89),
    label_position = c(lon = 55.4802, lat = -4.6767)
  ),
  SYR = list(
    name = "Syria", iso2 = "SY",
    bbox = c(xmin = 34.57, ymin = 31.41, xmax = 43.53, ymax = 38.23),
    label_position = c(lon = 38.2778, lat = 35.0066)
  ),
  TCA = list(
    name = "Turks and Caicos Islands", iso2 = "TC",
    bbox = c(xmin = -73.46, ymin = 20.39, xmax = -70.15, ymax = 22.86),
    label_position = c(lon = -71.7527, lat = 21.8166)
  ),
  TCD = list(
    name = "Chad", iso2 = "TD",
    bbox = c(xmin = 12.46, ymin = 6.55, xmax = 24.98, ymax = 24.35),
    label_position = c(lon = 18.645, lat = 15.143)
  ),
  TGO = list(
    name = "Togo", iso2 = "TG",
    bbox = c(xmin = -1.09, ymin = 5.2, xmax = 2.71, ymax = 12.04),
    label_position = c(lon = 1.0581, lat = 8.8072)
  ),
  THA = list(
    name = "Thailand", iso2 = "TH",
    bbox = c(xmin = 96.38, ymin = 4.73, xmax = 106.62, ymax = 21.35),
    label_position = c(lon = 101.0732, lat = 15.4597)
  ),
  TJK = list(
    name = "Tajikistan", iso2 = "TJ",
    bbox = c(xmin = 66.13, ymin = 35.78, xmax = 76.38, ymax = 41.94),
    label_position = c(lon = 72.5873, lat = 38.1998)
  ),
  TKM = list(
    name = "Turkmenistan", iso2 = "TM",
    bbox = c(xmin = 51.19, ymin = 34.24, xmax = 67.89, ymax = 43.69),
    label_position = c(lon = 58.6766, lat = 39.8552)
  ),
  TLS = list(
    name = "East Timor", iso2 = "TL",
    bbox = c(xmin = 123.11, ymin = -10.4, xmax = 128.23, ymax = -7.23),
    label_position = c(lon = 125.8547, lat = -8.8037)
  ),
  TON = list(
    name = "Tonga", iso2 = "TO",
    bbox = c(xmin = -177.2, ymin = -23.24, xmax = -172.93, ymax = -14.66),
    label_position = c(lon = -175.163, lat = -21.21)
  ),
  TTO = list(
    name = "Trinidad and Tobago", iso2 = "TT",
    bbox = c(xmin = -62.85, ymin = 9.14, xmax = -59.6, ymax = 12.25),
    label_position = c(lon = -60.9184, lat = 10.9989)
  ),
  TUN = list(
    name = "Tunisia", iso2 = "TN",
    bbox = c(xmin = 6.33, ymin = 29.33, xmax = 12.71, ymax = 38.25),
    label_position = c(lon = 9.0079, lat = 33.6873)
  ),
  TUR = list(
    name = "Turkey", iso2 = "TR",
    bbox = c(xmin = 24.43, ymin = 34.92, xmax = 46.04, ymax = 43),
    label_position = c(lon = 34.5083, lat = 39.3454)
  ),
  TUV = list(
    name = "Tuvalu", iso2 = "TV",
    bbox = c(xmin = 175.21, ymin = -10.32, xmax = 180.82, ymax = -4.77),
    label_position = c(lon = 179.2096, lat = -8.5137)
  ),
  TWN = list(
    name = "Taiwan", iso2 = "CN-TW",
    bbox = c(xmin = 117.27, ymin = 21, xmax = 123.01, ymax = 26.19),
    label_position = c(lon = 120.8682, lat = 23.6524)
  ),
  TZA = list(
    name = "United Republic of Tanzania", iso2 = "TZ",
    bbox = c(xmin = 28.4, ymin = -12.63, xmax = 41.37, ymax = -0.08),
    label_position = c(lon = 34.9592, lat = -6.0519)
  ),
  UGA = list(
    name = "Uganda", iso2 = "UG",
    bbox = c(xmin = 28.64, ymin = -2.38, xmax = 35.91, ymax = 5.12),
    label_position = c(lon = 32.9486, lat = 1.9726)
  ),
  UKR = list(
    name = "Ukraine", iso2 = "UA",
    bbox = c(xmin = 20.63, ymin = 44.31, xmax = 41.67, ymax = 53.27),
    label_position = c(lon = 32.1409, lat = 49.7247)
  ),
  UMI = list(
    name = "United States Minor Outlying Islands", iso2 = "UM",
    bbox = c(xmin = -163.31, ymin = -1.29, xmax = -159.1, ymax = 7.35),
    bbox_full = c(xmin = 165.59, ymin = -1.29, xmax = 286.03, ymax = 29.12),
    label_position = c(lon = -162.0765, lat = 5.8808)
  ),
  URY = list(
    name = "Uruguay", iso2 = "UY",
    bbox = c(xmin = -59.55, ymin = -35.88, xmax = -52, ymax = -29.19),
    label_position = c(lon = -55.9669, lat = -32.9611)
  ),
  USA = list(
    name = "United States of America", iso2 = "US",
    bbox = c(xmin = -126.14, ymin = 23.64, xmax = -65.57, ymax = 50.27),
    bbox_full = c(xmin = 169.51, ymin = 18, xmax = 295.98, ymax = 72.32),
    label_position = c(lon = -97.4826, lat = 39.5385)
  ),
  UZB = list(
    name = "Uzbekistan", iso2 = "UZ",
    bbox = c(xmin = 54.67, ymin = 36.28, xmax = 74.46, ymax = 46.46),
    label_position = c(lon = 64.0054, lat = 41.6936)
  ),
  VAT = list(
    name = "Vatican", iso2 = "VA",
    bbox = c(xmin = 11.22, ymin = 41, xmax = 13.68, ymax = 42.81),
    label_position = c(lon = 12.4534, lat = 41.9033)
  ),
  VCT = list(
    name = "Saint Vincent and the Grenadines", iso2 = "VC",
    bbox = c(xmin = -62.39, ymin = 11.68, xmax = -60.19, ymax = 14.28),
    label_position = c(lon = -61.3359, lat = 13.0879)
  ),
  VEN = list(
    name = "Venezuela", iso2 = "VE",
    bbox = c(xmin = -74.33, ymin = -0.25, xmax = -58.87, ymax = 16.61),
    label_position = c(lon = -64.5994, lat = 7.1825)
  ),
  VGB = list(
    name = "British Virgin Islands", iso2 = "VG",
    bbox = c(xmin = -65.73, ymin = 17.43, xmax = -63.31, ymax = 19.65),
    label_position = c(lon = -64.6366, lat = 18.4266)
  ),
  VIR = list(
    name = "United States Virgin Islands", iso2 = "VI",
    bbox = c(xmin = -66, ymin = 16.78, xmax = -63.6, ymax = 19.29),
    label_position = c(lon = -64.7792, lat = 17.7467)
  ),
  VNM = list(
    name = "Vietnam", iso2 = "VN",
    bbox = c(xmin = 101.13, ymin = 7.66, xmax = 110.46, ymax = 24.27),
    label_position = c(lon = 105.3873, lat = 21.7154)
  ),
  VUT = list(
    name = "Vanuatu", iso2 = "VU",
    bbox = c(xmin = 165.55, ymin = -21.16, xmax = 170.87, ymax = -12.16),
    label_position = c(lon = 166.9088, lat = -15.3715)
  ),
  WLF = list(
    name = "Wallis and Futuna", iso2 = "WF",
    bbox = c(xmin = -179.12, ymin = -15.22, xmax = -175.19, ymax = -12.31),
    label_position = c(lon = -178.1374, lat = -14.2864)
  ),
  WSM = list(
    name = "Samoa", iso2 = "WS",
    bbox = c(xmin = -173.72, ymin = -14.96, xmax = -170.5, ymax = -12.56),
    label_position = c(lon = -172.4382, lat = -13.6391)
  ),
  XKX = list(
    name = "Kosovo", iso2 = NA_character_,
    bbox = c(xmin = 18.77, ymin = 40.94, xmax = 23.03, ymax = 44.17),
    label_position = c(lon = 20.8607, lat = 42.5936)
  ),
  YEM = list(
    name = "Yemen", iso2 = "YE",
    bbox = c(xmin = 41.59, ymin = 11.21, xmax = 55.5, ymax = 19.9),
    label_position = c(lon = 45.8744, lat = 15.3282)
  ),
  ZAF = list(
    name = "South Africa", iso2 = "ZA",
    bbox = c(xmin = 15.36, ymin = -35.73, xmax = 34, ymax = -21.22),
    bbox_full = c(xmin = 15.13, ymin = -47.87, xmax = 39.32, ymax = -21.22),
    label_position = c(lon = 23.6657, lat = -29.7088)
  ),
  ZMB = list(
    name = "Zambia", iso2 = "ZM",
    bbox = c(xmin = 21.02, ymin = -18.97, xmax = 34.63, ymax = -7.29),
    label_position = c(lon = 26.3953, lat = -14.6608)
  ),
  ZWE = list(
    name = "Zimbabwe", iso2 = "ZW",
    bbox = c(xmin = 24.24, ymin = -23.3, xmax = 34.03, ymax = -14.71),
    label_position = c(lon = 29.9254, lat = -18.9116)
  )
)
