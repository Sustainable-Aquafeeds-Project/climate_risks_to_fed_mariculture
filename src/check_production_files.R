# check_production_files.R
# Validates the structure and integrity of production result files.
#
# Each file should be a list of 75 elements named "2025"-"2099".
# Each element should contain:
#   - $days: a numeric vector of length N
#   - $stats: a list of 13 response matrices, each with N rows
#
# The 13 expected response matrices:
#   weight, dw, T_response, rel_feeding, ing_act, anab, catab,
#   total_excr, metab, weight_scaled, dw_scaled, ing_act_scaled, total_excr_scaled

library(here)
library(qs2)

source(here("src", "dirs.R"))

# ---------- Configuration ----------
expected_years <- as.character(2025:2099)
expected_n_years <- length(expected_years) # 75

expected_response_names <- c("weight", "dw", "T_response", "rel_feeding", "ing_act", "anab", "catab", "total_excr", "metab", "weight_scaled", "dw_scaled", "ing_act_scaled", "total_excr_scaled")
expected_n_responses <- length(expected_response_names) # 13

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
}
