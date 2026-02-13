# ============================================================================
# Helper Functions for Twin Study Analysis
# ============================================================================

#' Extract and Format Goodness-of-Fit Statistics
#'
#' @param fitted_model A fitted mxModel object
#' @return Invisible (prints results to console)
fitGofs <- function(fitted_model) {
  summary_obj <- summary(fitted_model)
  
  cat("\n=== Goodness-of-Fit Statistics ===\n")
  cat(sprintf("  -2 Log-Likelihood: %.2f\n", summary_obj$Minus2LogLikelihood))
  cat(sprintf("  Degrees of Freedom: %d\n", summary_obj$degreesOfreedom))
  cat(sprintf("  AIC: %.2f\n", summary_obj$AIC))
  cat(sprintf("  BIC: %.2f\n", summary_obj$BIC))
  
  if (!is.null(summary_obj$chi)) {
    cat(sprintf("  Chi-square: %.2f\n", summary_obj$chi))
    cat(sprintf("  p-value: %.4f\n", summary_obj$p))
  }
  
  if (!is.null(summary_obj$CFI)) {
    cat(sprintf("  CFI: %.3f\n", summary_obj$CFI))
  }
  if (!is.null(summary_obj$TLI)) {
    cat(sprintf("  TLI: %.3f\n", summary_obj$TLI))
  }
  if (!is.null(summary_obj$RMSEA)) {
    cat(sprintf("  RMSEA: %.3f\n", summary_obj$RMSEA))
  }
  
  cat("\n")
  invisible()
}


#' Extract and Format Parameter Estimates
#'
#' @param fitted_model A fitted mxModel object
#' @return Invisible (prints results to console)
fitEsts <- function(fitted_model) {
  params <- fitted_model$output$estimate
  
  cat("\n=== Parameter Estimates ===\n")
  if (!is.null(params) && length(params) > 0) {
    param_df <- data.frame(
      Parameter = names(params),
      Estimate = params,
      row.names = NULL
    )
    print(param_df, row.names = FALSE)
  } else {
    cat("  No parameter estimates available.\n")
  }
  
  cat("\n")
  invisible()
}


#' Reshape Wide Twin Data to Long Format
#'
#' @param data Data frame with twin data in wide format
#' @param id Column name for twin pair identifier
#' @return Data frame in long format
fast.reshape <- function(data, id = "tvparnr") {
  # This is a placeholder for the mets::fast.reshape function
  # If mets package is not available, implement basic reshape here
  
  if (requireNamespace("mets", quietly = TRUE)) {
    return(mets::fast.reshape(data, id = id))
  } else {
    stop("Package 'mets' is required for fast.reshape.\n",
         "Install with: install.packages('mets')")
  }
}


#' Scale Variables Using Grand Mean and SD
#'
#' @param data Data frame containing twin data
#' @param var1 Column name for twin 1 variable
#' @param var2 Column name for twin 2 variable
#' @return Data frame with scaled variables
scale_grand <- function(data, var1, var2) {
  # Combine both twins for grand scaling
  all_vals <- c(data[[var1]], data[[var2]])
  grand_mean <- mean(all_vals, na.rm = TRUE)
  grand_sd <- sd(all_vals, na.rm = TRUE)
  
  # Scale using grand stats
  data[[var1]] <- (data[[var1]] - grand_mean) / grand_sd
  data[[var2]] <- (data[[var2]] - grand_mean) / grand_sd
  
  return(data)
}


#' Check Data Quality
#'
#' @param data Data frame to check
#' @param vars Variables to check
#' @param threshold Missing data threshold (percent)
#' @return Invisible (prints warnings)
check_data_quality <- function(data, vars, threshold = 20) {
  cat("\n=== Data Quality Checks ===\n")
  
  for (var in vars) {
    if (!var %in% colnames(data)) {
      warning(sprintf("Variable '%s' not found in data.", var))
      next
    }
    
    n_missing <- sum(is.na(data[[var]]))
    pct_missing <- (n_missing / nrow(data)) * 100
    
    if (pct_missing > threshold) {
      warning(sprintf("Variable '%s': %.1f%% missing data (> %.1f%% threshold)",
                      var, pct_missing, threshold))
    } else {
      cat(sprintf("  %s: %.1f%% missing\n", var, pct_missing))
    }
  }
  
  cat("\n")
  invisible()
}


#' Print Model Comparison Table
#'
#' @param comparison mxCompare output object
#' @return Invisible (prints formatted table)
print_model_comparison <- function(comparison) {
  cat("\n=== Model Comparison ===\n")
  print(comparison, digits = 3)
  cat("\n")
  invisible()
}


#' Save Results to CSV
#'
#' @param result_df Data frame with results
#' @param filename Output filename
#' @param output_dir Output directory
#' @return Path to saved file
save_results <- function(result_df, filename, output_dir = "output/results") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  filepath <- file.path(output_dir, filename)
  write.csv(result_df, filepath, row.names = FALSE)
  
  cat(sprintf("Results saved: %s\n", filepath))
  return(filepath)
}


#' Run Model with Error Handling
#'
#' @param model mxModel object
#' @param intervals Compute confidence intervals?
#' @param use_tryhard Use mxTryHard instead of mxRun?
#' @param max_attempts Number of optimization attempts (for mxTryHard)
#' @return Fitted model or NULL if failed
run_model_safe <- function(model, intervals = FALSE, use_tryhard = TRUE, max_attempts = 5) {
  tryCatch({
    if (use_tryhard) {
      fitted <- mxTryHard(model, intervals = intervals, extraTries = max_attempts - 1)
    } else {
      fitted <- mxRun(model, intervals = intervals)
    }
    
    # Check convergence
    if (fitted$output$status$code != 0) {
      warning(sprintf("Model did not converge properly. Status: %d - %s",
                      fitted$output$status$code,
                      fitted$output$status$statusMsg))
    }
    
    return(fitted)
  }, error = function(e) {
    warning(sprintf("Model fitting failed: %s", e$message))
    return(NULL)
  })
}


#' Extract Correlation from Covariance
#'
#' @param cov Covariance value
#' @param var1 Variance 1
#' @param var2 Variance 2
#' @return Correlation coefficient
cov2cor_value <- function(cov, var1, var2) {
  if (var1 <= 0 || var2 <= 0) {
    warning("Non-positive variance detected")
    return(NA)
  }
  return(cov / sqrt(var1 * var2))
}


#' Print Session Info for Reproducibility
#'
#' @return Invisible (prints session info)
print_session_info <- function() {
  cat("\n=== Session Information ===\n")
  print(sessionInfo())
  cat("\n")
  invisible()
}
