# check_production_files.R
# Validates the structure and integrity of production result files.
#
# Each file should be a list of 75 elements named "2025"-"2099".
# Each element should contain:
#   - $days: a numeric vector of length N
#   - $stats: a list of 15 response matrices, each with N rows
#
# The 15 expected response matrices:
#   weight, dw, T_response, rel_feeding, ing_act, anab, catab,
#   total_excr, metab, weight_scaled, dw_scaled, ing_act_scaled, total_excr_scaled,
#   food_prov, food_prov_scaled
#
# Additional data quality checks (Check 6):
#   - Weight must increase at some point (not monotonically decreasing)
#   - dw (weight change) must not be negative for the entire growth period
#   - rel_feeding and ing_act must be mostly numeric (not predominantly NaN/NA)
#   - anab, total_excr, food_prov, and food_prov_scaled must not be all zeros
#   If feeding issues are detected, the issue is reported as FEEDING_FAILED;
#   otherwise, growth-only issues are reported as NO_GROWTH.
#
# Additional data quality check (Check 7):
#   - Harvest biomass (weight_scaled, last finite value of the year) must not
#     increase by more than 250% from 2025 to 2099
#     (flagged as BIOMASS_INCREASE_EXCESSIVE), nor decrease by more than 90%
#     over the same period (flagged as BIOMASS_DECREASE_EXCESSIVE).

library(here)
library(qs2)

source(here("src", "dirs.R"))

# ---------- Configuration ----------
expected_years <- as.character(2025:2099)
expected_n_years <- length(expected_years) # 75

expected_response_names <- c("weight", "dw", "T_response", "rel_feeding", "ing_act", "anab", "catab", "total_excr", "metab", "weight_scaled", "dw_scaled", "ing_act_scaled", "total_excr_scaled", "food_prov", "food_prov_scaled")
expected_n_responses <- length(expected_response_names) # 15

# ---------- Find production files ----------
if (length(production_files) == 0) {
  stop("No production files found matching 'MC5000.qs' in: ", prod_path)
}

cat(sprintf("Found %d production files to check.\n\n", length(production_files)))

# ---------- Storage for issues ----------
issues <- list()

add_issue <- function(file, year = NA, issue_type, detail) {
  issues[[length(issues) + 1]] <<- data.frame(
    file = basename(file),
    year = as.character(year),
    issue_type = issue_type,
    detail = detail,
    stringsAsFactors = FALSE
  )
}

# ---------- Check each file ----------
for (i in seq_along(production_files)) {
  fpath <- production_files[i]
  fname <- basename(fpath)
  
  if (i %% 50 == 0 || i == 1) {
    cat(sprintf("  Checking file %d / %d ...\n", i, length(production_files)))
  }
  
  # Try to read the file
  dat <- tryCatch(
    qs_read(fpath),
    error = function(e) {
      add_issue(fpath, NA, "READ_ERROR", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(dat)) next
  
  # --- Check 1: Number and naming of top-level elements ---
  n_elements <- length(dat)
  element_names <- names(dat)
  
  if (n_elements != expected_n_years) {
    add_issue(fpath, NA, "WRONG_N_ELEMENTS",
              sprintf("Expected %d elements, got %d", expected_n_years, n_elements))
  }
  
  missing_years <- setdiff(expected_years, element_names)
  extra_years <- setdiff(element_names, expected_years)
  
  if (length(missing_years) > 0) {
    add_issue(fpath, NA, "MISSING_YEARS",
              paste("Missing:", paste(missing_years, collapse = ", ")))
  }
  if (length(extra_years) > 0) {
    add_issue(fpath, NA, "EXTRA_YEARS",
              paste("Extra:", paste(extra_years, collapse = ", ")))
  }
  
  # --- Check each year element ---
  # Collect all days vectors to compare consistency across years
  days_vectors <- list()
  
  for (yr in element_names) {
    yr_dat <- dat[[yr]]
    
    # Check 2: Does "days" exist?
    if (!"days" %in% names(yr_dat)) {
      add_issue(fpath, yr, "MISSING_DAYS", "No 'days' element found")
      next # can't check length matching without days
    }
    
    days_vec <- yr_dat[["days"]]
    n_days <- length(days_vec)
    days_vectors[[yr]] <- days_vec
    
    # Check 3: Does "stats" exist and have all 13 responses?
    if (!"stats" %in% names(yr_dat)) {
      add_issue(fpath, yr, "MISSING_STATS", "No 'stats' element found")
      next
    }
    
    stats <- yr_dat[["stats"]]
    stat_names <- names(stats)
    
    missing_responses <- setdiff(expected_response_names, stat_names)
    if (length(missing_responses) > 0) {
      add_issue(fpath, yr, "MISSING_RESPONSES",
                paste("Missing:", paste(missing_responses, collapse = ", ")))
    }
    
    # Check 4: Do all present response matrices have nrow == length(days) and ncol == 2?
    for (resp in intersect(expected_response_names, stat_names)) {
      resp_mat <- stats[[resp]]
      
      if (!is.matrix(resp_mat)) {
        add_issue(fpath, yr, "NOT_A_MATRIX",
                  sprintf("'%s' is not a matrix (class: %s)",
                          resp, paste(class(resp_mat), collapse = ", ")))
        next
      }
      
      n_rows <- nrow(resp_mat)
      n_cols <- ncol(resp_mat)
      
      if (n_rows != n_days) {
        add_issue(fpath, yr, "ROW_MISMATCH",
                  sprintf("'%s' has %d rows but 'days' has length %d",
                          resp, n_rows, n_days))
      }
      
      if (n_cols != 2) {
        add_issue(fpath, yr, "WRONG_NCOL",
                  sprintf("'%s' has %d columns, expected 2",
                          resp, n_cols))
      }
    }
    
    # Check 6: Data quality — detect feeding failure and no-growth conditions
    # Only run if the required matrices are present and valid
    has_weight    <- "weight"      %in% stat_names && is.matrix(stats[["weight"]])
    has_dw        <- "dw"          %in% stat_names && is.matrix(stats[["dw"]])
    has_relfeed   <- "rel_feeding" %in% stat_names && is.matrix(stats[["rel_feeding"]])
    has_ingact    <- "ing_act"     %in% stat_names && is.matrix(stats[["ing_act"]])
    has_anab      <- "anab"        %in% stat_names && is.matrix(stats[["anab"]])
    has_totalexcr <- "total_excr"  %in% stat_names && is.matrix(stats[["total_excr"]])
    has_foodprov  <- "food_prov"   %in% stat_names && is.matrix(stats[["food_prov"]])
    has_foodprovs <- "food_prov_scaled" %in% stat_names && is.matrix(stats[["food_prov_scaled"]])
    
    # Feeding checks (use column 1 = mean)
    feeding_failed <- FALSE
    feeding_details <- character(0)
    
    if (has_relfeed) {
      rf_vals <- stats[["rel_feeding"]][, 1]
      rf_na_frac <- sum(!is.finite(rf_vals))
      if (rf_na_frac > 100) {
        feeding_failed <- TRUE
        feeding_details <- c(feeding_details,
          sprintf("there are %.0f NA rel_feeding values", rf_na_frac))
      }
    }
    
    if (has_ingact) {
      ia_vals <- stats[["ing_act"]][, 1]
      ia_na_frac <- sum(!is.finite(ia_vals))
      if (ia_na_frac > 1000) {
        feeding_failed <- TRUE
        feeding_details <- c(feeding_details,
          sprintf("there are %.0f NA ing_act values", ia_na_frac))
      }
    }
    
    if (has_anab) {
      anab_vals <- stats[["anab"]][, 1]
      if (all(anab_vals == 0, na.rm = TRUE)) {
        feeding_failed <- TRUE
        feeding_details <- c(feeding_details, "anab is all zeros")
      }
    }
    
    if (has_totalexcr) {
      te_vals <- stats[["total_excr"]][, 1]
      if (all(te_vals == 0, na.rm = TRUE)) {
        feeding_failed <- TRUE
        feeding_details <- c(feeding_details, "total_excr is all zeros")
      }
    }
    
    if (has_foodprov) {
      fp_vals <- stats[["food_prov"]][, 1]
      if (all(fp_vals == 0, na.rm = TRUE)) {
        feeding_failed <- TRUE
        feeding_details <- c(feeding_details, "food_prov is all zeros")
      }
    }
    
    if (has_foodprovs) {
      fps_vals <- stats[["food_prov_scaled"]][, 1]
      if (all(fps_vals == 0, na.rm = TRUE)) {
        feeding_failed <- TRUE
        feeding_details <- c(feeding_details, "food_prov_scaled is all zeros")
      }
    }
    
    # Growth checks (use column 1 = mean)
    no_growth <- FALSE
    growth_details <- character(0)
    
    if (has_weight) {
      wt_vals <- stats[["weight"]][, 1]
      wt_finite <- wt_vals[is.finite(wt_vals)]
      if (length(wt_finite) > 1 && all(diff(wt_finite) <= 0)) {
        no_growth <- TRUE
        growth_details <- c(growth_details, "weight never increases")
      }
    }
    
    if (has_dw) {
      dw_vals <- stats[["dw"]][, 1]
      dw_finite <- dw_vals[is.finite(dw_vals)]
      if (length(dw_finite) > 0 && all(dw_finite < 0)) {
        no_growth <- TRUE
        growth_details <- c(growth_details, "dw is negative throughout")
      }
    }
    
    # Report: feeding failure supersedes generic no-growth
    if (feeding_failed) {
      all_details <- c(feeding_details, growth_details)
      add_issue(fpath, yr, "FEEDING_FAILED",
                paste("Feeding failed, fish not growing:",
                      paste(all_details, collapse = "; ")))
    } else if (no_growth) {
      add_issue(fpath, yr, "NO_GROWTH",
                paste("Fish are not growing:",
                      paste(growth_details, collapse = "; ")))
    }
  }
  
  # --- Check 7: Harvest biomass (weight_scaled) must not change excessively from 2025 to 2099 ---
  if (all(c("2025", "2099") %in% element_names)) {
    stats_2025 <- dat[["2025"]][["stats"]]
    stats_2099 <- dat[["2099"]][["stats"]]
    
    has_ws_2025 <- !is.null(stats_2025) && "weight_scaled" %in% names(stats_2025) && is.matrix(stats_2025[["weight_scaled"]])
    has_ws_2099 <- !is.null(stats_2099) && "weight_scaled" %in% names(stats_2099) && is.matrix(stats_2099[["weight_scaled"]])
    
    if (has_ws_2025 && has_ws_2099) {
      ws_2025_finite <- stats_2025[["weight_scaled"]][, 1]
      ws_2025_finite <- ws_2025_finite[is.finite(ws_2025_finite)]
      ws_2099_finite <- stats_2099[["weight_scaled"]][, 1]
      ws_2099_finite <- ws_2099_finite[is.finite(ws_2099_finite)]
      
      if (length(ws_2025_finite) > 0 && length(ws_2099_finite) > 0) {
        # Harvest biomass: the final (end-of-year) biomass value for each year
        harvest_biomass_2025 <- ws_2025_finite[length(ws_2025_finite)]
        harvest_biomass_2099 <- ws_2099_finite[length(ws_2099_finite)]
        
        if (is.finite(harvest_biomass_2025) && harvest_biomass_2025 > 0) {
          pct_change <- (harvest_biomass_2099 - harvest_biomass_2025) / harvest_biomass_2025 * 100
          if (pct_change > 250) {
            add_issue(fpath, NA, "BIOMASS_INCREASE_EXCESSIVE",
                      sprintf("Harvest biomass (weight_scaled) increased by %.1f%% from 2025 (%.4f) to 2099 (%.4f)",
                              pct_change, harvest_biomass_2025, harvest_biomass_2099))
          } else if (pct_change < -90) {
            add_issue(fpath, NA, "BIOMASS_DECREASE_EXCESSIVE",
                      sprintf("Harvest biomass (weight_scaled) decreased by %.1f%% from 2025 (%.4f) to 2099 (%.4f)",
                              pct_change, harvest_biomass_2025, harvest_biomass_2099))
          }
        }
      }
    }
  }
  
  # --- Check 5: days vectors must be identical across all years ---
  if (length(days_vectors) > 1) {
    ref_days <- days_vectors[[1]]
    ref_yr <- names(days_vectors)[1]
    for (k in seq_along(days_vectors)[-1]) {
      comp_yr <- names(days_vectors)[k]
      comp_days <- days_vectors[[k]]
      if (!identical(ref_days, comp_days)) {
        # Determine the nature of the difference
        if (length(ref_days) != length(comp_days)) {
          add_issue(fpath, comp_yr, "DAYS_LENGTH_DIFFERS",
                    sprintf("'days' has length %d but year %s has length %d",
                            length(comp_days), ref_yr, length(ref_days)))
        } else {
          add_issue(fpath, comp_yr, "DAYS_VALUES_DIFFER",
                    sprintf("'days' values differ from year %s (same length %d)",
                            ref_yr, length(ref_days)))
        }
      }
    }
  }
}

# ---------- Compile and report results ----------
if (length(issues) == 0) {
  cat("\n========================================\n")
  cat("ALL CHECKS PASSED\n")
  cat(sprintf("All %d production files have correct structure.\n", length(production_files)))
  cat("========================================\n")
} else {
  issues_df <- do.call(rbind, issues)
  
  cat("\n========================================\n")
  cat("ISSUES FOUND\n")
  cat("========================================\n\n")
  
  # Summary by issue type
  issue_summary <- table(issues_df$issue_type)
  cat("Summary of issues:\n")
  for (it in names(issue_summary)) {
    cat(sprintf("  %-20s : %d\n", it, issue_summary[it]))
  }
  cat(sprintf("\nTotal issues: %d across %d files\n\n",
              nrow(issues_df),
              length(unique(issues_df$file))))
  
  # Detailed report
  cat("--- Detailed report ---\n\n")
  for (f in unique(issues_df$file)) {
    cat(sprintf("File: %s\n", f))
    f_issues <- issues_df[issues_df$file == f, ]
    for (j in seq_len(nrow(f_issues))) {
      row <- f_issues[j, ]
      yr_label <- if (is.na(row$year)) "file-level" else paste0("year ", row$year)
      cat(sprintf("  [%s] %s - %s\n", yr_label, row$issue_type, row$detail))
    }
    cat("\n")
  }
  
  # Save to CSV for further inspection
  issues_csv <- file.path(outs_path, "production_check_issues.csv")
  write.csv(issues_df, issues_csv, row.names = FALSE)
  cat(sprintf("Issues saved to: %s\n", issues_csv))
  
  # ---------- High-error-count report ----------
  # File-level issues (year is NA) are more severe than year-level issues.
  # A file qualifies as "high error" if it has ANY file-level issue,
  # or more than 10 year-level issues.
  year_level_threshold <- 10
  
  file_level_counts <- aggregate(
    is.na(issues_df$year) ~ issues_df$file,
    FUN = sum
  )
  names(file_level_counts) <- c("file", "n_file_level_issues")
  
  year_level_counts <- aggregate(
    !is.na(issues_df$year) ~ issues_df$file,
    FUN = sum
  )
  names(year_level_counts) <- c("file", "n_year_level_issues")
  
  file_issue_counts <- merge(file_level_counts, year_level_counts, by = "file", all = TRUE)
  file_issue_counts[is.na(file_issue_counts)] <- 0
  
  high_error_files <- file_issue_counts[
    file_issue_counts$n_file_level_issues >= 1 |
      file_issue_counts$n_year_level_issues > year_level_threshold,
  ]
  high_error_files <- high_error_files[
    order(-high_error_files$n_file_level_issues, -high_error_files$n_year_level_issues),
  ]
  
  high_error_csv <- file.path(outs_path, "production_check_high_error_files.csv")
  write.csv(high_error_files, high_error_csv, row.names = FALSE)
  
  cat(sprintf("\n%d files flagged as high-error (any file-level issue, or > %d year-level issues) saved to: %s\n",
              nrow(high_error_files), year_level_threshold, high_error_csv))
}
