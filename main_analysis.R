# ============================================================================
# Twin Study Heritability Analysis - Main Pipeline
# ============================================================================
# This script runs the complete heritability analysis for twin data using
# saturated models with covariate adjustment in OpenMx.
# ============================================================================

# ---- Setup ----

# Clear workspace
rm(list = ls())

# Load required packages
required_packages <- c("OpenMx", "psych", "polycor", "mets")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.\n", pkg),
         sprintf("Install with: install.packages('%s')", pkg))
  }
}

# Load configuration and helper functions
source("config.R")
source("helper_functions.R")

cat("\n============================================================\n")
cat("Twin Study Heritability Analysis with OpenMx\n")
cat("============================================================\n\n")

# Print session info for reproducibility
print_session_info()

# ---- Load Data ----

cat("Loading data...\n")
twinData_raw <- read.csv(file = DATA_PATH, header = TRUE, sep = ",")
cat(sprintf("  Data loaded: %d rows, %d columns\n", nrow(twinData_raw), ncol(twinData_raw)))

# Descriptive statistics
if (PRINT_DESCRIPTIVES) {
  cat("\nDataset Overview:\n")
  print(str(twinData_raw))
}

# ---- Analysis Loop ----

cat("\n============================================================\n")
cat(sprintf("Starting analysis for %d phenotypes...\n", length(PHENOTYPES)))
cat("============================================================\n\n")

for (phenotype in PHENOTYPES) {
  
  cat(sprintf("\n--- Analyzing: %s ---\n", phenotype))
  
  # Define variable names
  selVars <- c(paste0(phenotype, "1"), paste0(phenotype, "2"))
  covVars <- c(AGE_COLS, SEX_COLS, ETIV_COLS)
  
  # Reload and prepare data for this phenotype
  twinData <- twinData_raw
  twinData <- fast.reshape(twinData, id = TWIN_ID_COL)
  
  # Remove cases with missing values
  complete_cases <- complete.cases(twinData[, selVars[1]]) & 
                    complete.cases(twinData[, selVars[2]])
  rows_with_nan <- which(!complete_cases)
  
  if (length(rows_with_nan) > 0) {
    tvparnrs_with_nan <- twinData[[TWIN_ID_COL]][rows_with_nan]
    twinData <- twinData[!(twinData[[TWIN_ID_COL]] %in% tvparnrs_with_nan), ]
    cat(sprintf("  Removed %d pairs with missing data\n", length(rows_with_nan)))
  }
  
  # Scale variables using grand mean/SD
  if (GRAND_SCALE) {
    # Scale phenotype
    twinData <- scale_grand(twinData, selVars[1], selVars[2])
    
    # Scale covariates
    twinData <- scale_grand(twinData, AGE_COLS[1], AGE_COLS[2])
    twinData <- scale_grand(twinData, ETIV_COLS[1], ETIV_COLS[2])
    
    cat("  Variables scaled by grand mean/SD\n")
  }
  
  # Split by zygosity
  mzData <- subset(twinData, get(ZYGOSITY_COL) == 1, c(selVars, covVars))
  dzData <- subset(twinData, get(ZYGOSITY_COL) == 2, c(selVars, covVars))
  mzData <- as.data.frame(mzData)
  dzData <- as.data.frame(dzData)
  
  cat(sprintf("  MZ pairs: %d\n", nrow(mzData)))
  cat(sprintf("  DZ pairs: %d\n", nrow(dzData)))
  
  # Check minimum sample size
  if (nrow(mzData) < MIN_SAMPLE_SIZE || nrow(dzData) < MIN_SAMPLE_SIZE) {
    warning(sprintf("Insufficient sample size for %s. Skipping.", phenotype))
    next
  }
  
  # Define output files
  output_txt <- file.path(LOG_DIR, paste0(phenotype, "_saturated.txt"))
  output_tests_csv <- file.path(RESULTS_DIR, paste0(phenotype, "_model_tests.csv"))
  output_cov_csv <- file.path(RESULTS_DIR, paste0(phenotype, "_covariate_tests.csv"))
  
  # Start logging
  sink(file = output_txt)
  
  cat(sprintf("\n=== Analysis: %s ===\n", phenotype))
  cat(sprintf("Date: %s\n", Sys.time()))
  cat(sprintf("MZ sample: %d pairs\n", nrow(mzData)))
  cat(sprintf("DZ sample: %d pairs\n", nrow(dzData)))
  
  # Print descriptive statistics
  cat("\n--- Descriptive Statistics ---\n")
  cat("\nMZ Means:\n")
  print(colMeans(mzData, na.rm = TRUE))
  cat("\nDZ Means:\n")
  print(colMeans(dzData, na.rm = TRUE))
  cat("\nMZ Covariance Matrix:\n")
  print(cov(mzData, use = "complete"))
  cat("\nDZ Covariance Matrix:\n")
  print(cov(dzData, use = "complete"))
  
  # ---- Model Specification ----
  
  nv <- 1  # Number of variables per twin
  nt <- 2  # Number of twins
  ntv <- nv * nt
  
  # Covariate definition matrices
  defSex <- mxMatrix(type = "Full", nrow = 1, ncol = nt, free = FALSE, 
                     labels = paste0("data.", SEX_COLS), name = "Sex")
  defAge <- mxMatrix(type = "Full", nrow = 1, ncol = nt, free = FALSE, 
                     labels = paste0("data.", AGE_COLS), name = "Age")
  defEtiv <- mxMatrix(type = "Full", nrow = 1, ncol = nt, free = FALSE, 
                      labels = paste0("data.", ETIV_COLS), name = "eTIV")
  
  # Regression coefficient matrices
  betaSex <- mxMatrix(type = "Full", nrow = 1, ncol = nv, free = TRUE, 
                      values = SV_BETA_SEX, labels = "betaS", name = "bS")
  betaAge <- mxMatrix(type = "Full", nrow = 1, ncol = nv, free = TRUE, 
                      values = SV_BETA_AGE, labels = "betaA", name = "bA")
  betaEtiv <- mxMatrix(type = "Full", nrow = 1, ncol = nv, free = TRUE, 
                       values = SV_BETA_ETIV, labels = "betaE", name = "bE")
  
  # Mean matrices
  meanMZ <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                     values = SV_MEANS, labels = c("mMZ1", "mMZ2"), name = "meanMZ")
  meanDZ <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, 
                     values = SV_MEANS, labels = c("mDZ1", "mDZ2"), name = "meanDZ")
  
  # Expected means with covariates
  expMeanMZ <- mxAlgebra(expression = meanMZ + Sex %x% bS + Age %x% bA + eTIV %x% bE,
                         name = "expMeanMZ")
  expMeanDZ <- mxAlgebra(expression = meanDZ + Sex %x% bS + Age %x% bA + eTIV %x% bE,
                         name = "expMeanDZ")
  
  # Covariance matrices
  covMZ <- mxMatrix(type = "Symm", nrow = ntv, ncol = ntv, free = TRUE,
                    values = c(SV_VARIANCE, SV_VARIANCE * 0.5, SV_VARIANCE),
                    lbound = c(LB_VARIANCE, LB_CORRELATION, LB_VARIANCE),
                    labels = c("vMZ1", "cMZ21", "vMZ2"),
                    name = "covMZ")
  
  covDZ <- mxMatrix(type = "Symm", nrow = ntv, ncol = ntv, free = TRUE,
                    values = c(SV_VARIANCE, SV_VARIANCE * 0.5, SV_VARIANCE),
                    lbound = c(LB_VARIANCE, LB_CORRELATION, LB_VARIANCE),
                    labels = c("vDZ1", "cDZ21", "vDZ2"),
                    name = "covDZ")
  
  # Data objects
  dataMZ <- mxData(observed = mzData, type = "raw")
  dataDZ <- mxData(observed = dzData, type = "raw")
  
  # Expectation objects
  expMZ <- mxExpectationNormal(covariance = "covMZ", means = "expMeanMZ", dimnames = selVars)
  expDZ <- mxExpectationNormal(covariance = "covDZ", means = "expMeanDZ", dimnames = selVars)
  
  # Fit function
  funML <- mxFitFunctionML()
  
  # Model objects
  pars <- list(betaAge, betaSex, betaEtiv)
  defs <- list(defAge, defSex, defEtiv)
  
  modelMZ <- mxModel(pars, defs, meanMZ, expMeanMZ, covMZ, dataMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, defs, meanDZ, expMeanDZ, covDZ, dataDZ, expDZ, funML, name = "DZ")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ"))
  
  # Confidence intervals
  ciCov <- mxCI(c("MZ.covMZ", "DZ.covDZ"))
  ciMean <- mxCI(c("MZ.meanMZ", "DZ.meanDZ"))
  
  # Build saturated model
  modelSAT <- mxModel("Saturated_Model", pars, modelMZ, modelDZ, multi, ciCov, ciMean)
  
  # ---- Fit Saturated Model ----
  
  cat("\n--- Fitting Saturated Model ---\n")
  fitSAT <- run_model_safe(modelSAT, intervals = FALSE, use_tryhard = USE_MXTRYHARD)
  
  if (is.null(fitSAT)) {
    warning(sprintf("Saturated model failed for %s. Skipping.", phenotype))
    sink()
    closeAllConnections()
    next
  }
  
  sumSAT <- summary(fitSAT)
  print(sumSAT)
  fitGofs(fitSAT)
  fitEsts(fitSAT)
  
  # ---- Fit Nested Models ----
  
  cat("\n--- Fitting Nested Models ---\n")
  
  # Equal means within zygosity
  modelEMO <- mxModel(fitSAT, name = "Equal_Means_Within_Zyg")
  modelEMO <- omxSetParameters(modelEMO, label = c("mMZ1", "mMZ2"), 
                               free = TRUE, values = SV_MEANS, newlabels = "mMZ")
  modelEMO <- omxSetParameters(modelEMO, label = c("mDZ1", "mDZ2"), 
                               free = TRUE, values = SV_MEANS, newlabels = "mDZ")
  fitEMO <- run_model_safe(modelEMO, intervals = FALSE, use_tryhard = USE_MXTRYHARD)
  
  # Equal variances within zygosity
  modelEMVO <- mxModel(fitEMO, name = "Equal_Var_Within_Zyg")
  modelEMVO <- omxSetParameters(modelEMVO, label = c("vMZ1", "vMZ2"), 
                                free = TRUE, values = SV_VARIANCE, newlabels = "vMZ")
  modelEMVO <- omxSetParameters(modelEMVO, label = c("vDZ1", "vDZ2"), 
                                free = TRUE, values = SV_VARIANCE, newlabels = "vDZ")
  fitEMVO <- run_model_safe(modelEMVO, intervals = FALSE, use_tryhard = USE_MXTRYHARD)
  
  # Equal means and variances across zygosity
  modelEMVZ <- mxModel(fitEMVO, name = "Equal_Means_Var_Across_Zyg")
  modelEMVZ <- omxSetParameters(modelEMVZ, label = c("mMZ", "mDZ"), 
                                free = TRUE, values = SV_MEANS, newlabels = "mZ")
  modelEMVZ <- omxSetParameters(modelEMVZ, label = c("vMZ", "vDZ"), 
                                free = TRUE, values = SV_VARIANCE, newlabels = "vZ")
  modelEMVZ <- mxModel(modelEMVZ,
                       mxAlgebra(MZ.covMZ[2, 1] / MZ.covMZ[1, 1], name = "rMZ"),
                       mxAlgebra(DZ.covDZ[2, 1] / DZ.covDZ[1, 1], name = "rDZ"),
                       mxCI(c("rMZ", "rDZ")))
  fitEMVZ <- run_model_safe(modelEMVZ, intervals = COMPUTE_CI, 
                           use_tryhard = USE_MXTRYHARD, max_attempts = MAX_ATTEMPTS)

  if (!is.null(fitEMVZ)) {
    rMZ_val <- mxEval(rMZ, fitEMVZ)
    rDZ_val <- mxEval(rDZ, fitEMVZ)
    cat(sprintf("\nEstimated twin correlations (rMZ, rDZ): %.4f, %.4f\n", rMZ_val, rDZ_val))

    if (COMPUTE_CI) {
      ci_table <- summary(fitEMVZ)$CI
      if (!is.null(ci_table)) {
        cat("\nCorrelation CIs (rMZ, rDZ):\n")
        print(ci_table[rownames(ci_table) %in% c("rMZ", "rDZ"), , drop = FALSE])
      }
    }
  }
  
  # ---- Model Comparison ----
  
  cat("\n--- Model Comparison: Nested Models ---\n")
  if (!is.null(fitEMO) && !is.null(fitEMVO) && !is.null(fitEMVZ)) {
    model_comparison <- mxCompare(fitSAT, list(fitEMO, fitEMVO, fitEMVZ))
    print(model_comparison)
    
    # Save model comparison
    save_results(as.data.frame(model_comparison), 
                basename(output_tests_csv), 
                dirname(output_tests_csv))
  }
  
  # ---- Covariate Tests ----
  
  cat("\n--- Testing Covariate Effects ---\n")
  
  modelNoAge <- omxSetParameters(fitEMVZ, labels = "betaA", values = 0, 
                                 free = FALSE, name = "No_Age")
  modelNoSex <- omxSetParameters(fitEMVZ, labels = "betaS", values = 0, 
                                 free = FALSE, name = "No_Sex")
  modelNoETIV <- omxSetParameters(fitEMVZ, labels = "betaE", values = 0, 
                                  free = FALSE, name = "No_eTIV")
  
  fitNoAge <- run_model_safe(modelNoAge, intervals = FALSE, use_tryhard = USE_MXTRYHARD)
  fitNoSex <- run_model_safe(modelNoSex, intervals = FALSE, use_tryhard = USE_MXTRYHARD)
  fitNoETIV <- run_model_safe(modelNoETIV, intervals = FALSE, use_tryhard = USE_MXTRYHARD)
  
  if (!is.null(fitNoAge) && !is.null(fitNoSex) && !is.null(fitNoETIV)) {
    cov_comparison <- mxCompare(fitEMVZ, list(fitNoAge, fitNoSex, fitNoETIV))
    print(cov_comparison)
    
    # Save covariate comparison
    save_results(as.data.frame(cov_comparison), 
                basename(output_cov_csv), 
                dirname(output_cov_csv))
  }
  
  # Stop logging
  sink()
  closeAllConnections()
  
  cat(sprintf("  Analysis complete for %s\n", phenotype))
  cat(sprintf("  Results saved to %s\n", output_txt))
  
}

# ---- Final Summary ----

cat("\n============================================================\n")
cat("Analysis Pipeline Complete\n")
cat("============================================================\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
cat("\nCheck the following directories:\n")
cat(sprintf("  - Logs: %s\n", LOG_DIR))
cat(sprintf("  - Results: %s\n", RESULTS_DIR))
cat("\n")
