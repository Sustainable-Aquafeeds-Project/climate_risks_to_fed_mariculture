#' Spatially combine comparable rasters or stacks into one
#'
#' Mosaics any number of `SpatRaster` objects (single- or multi-layer) that
#' cover different parts of the same area — e.g. the 1-3 FAO fishing regions
#' belonging to one country — into a single `SpatRaster` covering the union of
#' their extents. Layer structure is preserved, so the function is agnostic to
#' how you eventually organise species/year (one layer per stack, or many).
#'
#' @param ... `SpatRaster` objects, file paths, `raster`-package objects,
#'   `SpatRasterCollection`s, or (nested) lists of any of these. Anything empty
#'   or `NULL` is dropped.
#' @param fun How to resolve cells covered by more than one input:
#'   `"mean"`, `"median"`, `"modal"`, `"min"`, `"max"`, `"sum"` (passed to
#'   `terra::mosaic()`), or `"first"` / `"last"` (uses `terra::merge()`, which
#'   is cheaper when inputs are known not to overlap).
#' @param match_layers If `TRUE`, multi-layer inputs whose layer names are the
#'   same set but in a different order are reordered to match the first input.
#'   If `FALSE`, layers are combined in positional order.
#' @param align `"snap"` (default) resamples inputs whose cells don't line up
#'   onto the first input's grid. `"strict"` errors instead, which is safer if
#'   you don't expect misalignment.
#' @param method Resampling method used when `align = "snap"`. Defaults to
#'   `"near"` for categorical rasters and `"bilinear"` otherwise.
#' @param filename,overwrite Optional; passed to the writing step.
#' @param quiet Suppress informational messages.
#'
#' @return A single `SpatRaster` with the same number of layers and the same
#'   layer names as the first input.
combine_rasters <- function(...,
                            fun          = "mean",
                            match_layers = TRUE,
                            align        = c("snap", "strict"),
                            method       = NULL,
                            filename     = "",
                            overwrite    = FALSE,
                            quiet        = FALSE) {

  align <- match.arg(align)
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("`terra` is required.", call. = FALSE)
  }

  # ---- 1. gather and normalise inputs ------------------------------------
  flatten <- function(x) {
    if (inherits(x, "SpatRasterCollection")) {
      lapply(seq_len(length(x)), function(i) x[i])
    } else if (is.list(x)) {
      do.call(c, c(list(list()), lapply(x, flatten)))
    } else {
      list(x)
    }
  }

  rl <- do.call(c, c(list(list()), lapply(list(...), flatten)))
  rl <- Filter(function(x) !is.null(x), rl)
  if (!length(rl)) stop("No rasters supplied.", call. = FALSE)

  rl <- lapply(rl, function(x) {
    if (inherits(x, "SpatRaster")) x else terra::rast(x)  # paths, raster::* objs
  })

  # ---- 2. check comparability --------------------------------------------
  ref   <- rl[[1]]
  n_lyr <- terra::nlyr(ref)
  ref_n <- names(ref)

  for (i in seq_along(rl)) {
    x <- rl[[i]]

    if (!terra::same.crs(x, ref)) {
      stop("Input ", i, " has a different CRS to the first input. ",
           "Reproject deliberately (terra::project()) before combining, so ",
           "you control the method and target grid.", call. = FALSE)
    }
    if (terra::nlyr(x) != n_lyr) {
      stop("Input ", i, " has ", terra::nlyr(x), " layer(s) but the first has ",
           n_lyr, ". All inputs must have the same layer structure.",
           call. = FALSE)
    }
    if (!isTRUE(all.equal(terra::res(x), terra::res(ref), tolerance = 1e-6))) {
      stop("Input ", i, " has resolution ",
           paste(signif(terra::res(x), 6), collapse = " x "),
           " but the first has ",
           paste(signif(terra::res(ref), 6), collapse = " x "),
           ". Resample to a common resolution first.", call. = FALSE)
    }

    if (match_layers && n_lyr > 1L && !identical(names(x), ref_n)) {
      if (setequal(names(x), ref_n) && !anyDuplicated(names(x))) {
        if (!quiet) message("Reordering layers of input ", i,
                            " to match the first input.")
        rl[[i]] <- x[[ref_n]]
      } else {
        stop("Layer names of input ", i, " don't match the first input: ",
             paste(setdiff(names(x), ref_n), collapse = ", "),
             ". Fix the names, or set match_layers = FALSE to combine ",
             "positionally.", call. = FALSE)
      }
    }
  }

  # ---- 3. snap grids that don't line up -----------------------------------
  tol    <- terra::res(ref) * 1e-4
  o_ref  <- terra::origin(ref)

  for (i in seq_along(rl)) {
    if (any(abs(terra::origin(rl[[i]]) - o_ref) > tol)) {
      if (align == "strict") {
        stop("Input ", i, " is not aligned to the first input's grid ",
             "(different cell origin). Use align = 'snap' to resample it.",
             call. = FALSE)
      }
      m <- if (!is.null(method)) {
        method
      } else if (any(terra::is.factor(rl[[i]]))) "near" else "bilinear"

      if (!quiet) message("Snapping input ", i,
                          " onto the reference grid (method = '", m, "').")

      tmpl <- terra::rast(terra::align(terra::ext(rl[[i]]), ref, snap = "out"),
                          resolution = terra::res(ref),
                          crs        = terra::crs(ref))
      rl[[i]] <- terra::resample(rl[[i]], tmpl, method = m)
    }
  }

  # ---- 4. combine ----------------------------------------------------------
  if (length(rl) == 1L) {
    out <- rl[[1]]
    if (nzchar(filename)) {
      out <- terra::writeRaster(out, filename, overwrite = overwrite)
    }
    return(out)
  }

  sc <- terra::sprc(rl)

  out <- if (fun %in% c("first", "last")) {
    terra::merge(sc, first = identical(fun, "first"), na.rm = TRUE,
                 filename = filename, overwrite = overwrite)
  } else {
    terra::mosaic(sc, fun = fun, filename = filename, overwrite = overwrite)
  }

  names(out) <- ref_n
  if (terra::has.time(ref)) terra::time(out) <- terra::time(ref)
  out
}


# ---- usage ----------------------------------------------------------------
# Single-layer rasters, three FAO regions for one country:
#   aus <- combine_rasters(r_5701, r_5702, r_5703)
#
# Stacks (e.g. one layer per year), passed as a list:
#   aus <- combine_rasters(list(stk_5701, stk_5702, stk_5703))
#
# Straight from file paths:
#   aus <- combine_rasters(list.files("data/AUS", "\\.tif$", full.names = TRUE))
#
# Whole dataset, grouped by country. `combine_rasters()` accepts the list that
# lapply hands it, so this needs no wrapper:
#   by_country <- split(rasters, meta$country) |> lapply(combine_rasters)
#
# Regions that are known not to overlap — cheaper, and errs loudly if they do
# overlap in a way you didn't expect (you'd see a hard edge):
#   aus <- combine_rasters(r_5701, r_5702, fun = "first")
