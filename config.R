# ============================================================================
# Twin Study Heritability Analysis - Configuration File
# ============================================================================
# Edit this file to customize the analysis for your dataset
# ============================================================================

# ---- Data Configuration ----

# Path to twin data CSV file
DATA_PATH <- "data/twin_data.csv"

# Twin pair identifier column
TWIN_ID_COL <- "tvparnr"

# Zygosity column (1 = MZ, 2 = DZ)
ZYGOSITY_COL <- "zyg1"

# Covariate column names (must be in format: var1, var2 for each twin)
AGE_COLS <- c("age1", "age2")
SEX_COLS <- c("sex1", "sex2")
ETIV_COLS <- c("eTIV1", "eTIV2")

# Phenotypes to analyze (suffix 1,2 will be added automatically)
PHENOTYPES <- c(
  "Total_Cerebel_Vol",
  "Left_VIIB",
  "Corpus_Medullare",
  "Left_Crus_I",
  "Left_Crus_II",
  "Left_I_III",
  "Left_IV",
  "Left_IX",
  "Left_V",
  "Left_VI",
  "Left_VIIIA",
  "Left_VIIIB",
  "Left_X",
  "Right_Crus_I",
  "Right_Crus_II",
  "Right_I_III",
  "Right_IX",
  "Right_V",
  "Right_VI",
  "Right_VIIB",
  "Right_VIIIA",
  "Right_VIIIB",
  "Right_X",
  "Right_IV",
  "Vermis_IX",
  "Vermis_VI",
  "Vermis_VII",
  "Vermis_VIII",
  "Vermis_X"
)

# ---- Output Configuration ----

# Output directory for results
OUTPUT_DIR <- "output"

# Create subdirectories
LOG_DIR <- file.path(OUTPUT_DIR, "logs")
RESULTS_DIR <- file.path(OUTPUT_DIR, "results")

# ---- Model Configuration ----

# Starting values for model optimization
SV_MEANS <- 0          # Mean starting value
SV_VARIANCE <- 1.0     # Variance starting value

# Lower bounds for parameters
LB_VARIANCE <- 0.0001  # Variances must stay positive
LB_CORRELATION <- -0.99 # Correlations can be negative but bounded

# Starting values for regression coefficients
SV_BETA_AGE <- -0.1    # Age effect starting value
SV_BETA_SEX <- 0       # Sex effect starting value
SV_BETA_ETIV <- 0.2    # eTIV effect starting value

# ---- Optimization Configuration ----

# Use mxTryHard for robust optimization? (TRUE recommended)
USE_MXTRYHARD <- TRUE

# Number of optimization attempts (for mxTryHard)
MAX_ATTEMPTS <- 5

# Compute confidence intervals?
COMPUTE_CI <- TRUE

# ---- Data Quality Checks ----

# Print descriptive statistics?
PRINT_DESCRIPTIVES <- TRUE

# Minimum sample size per group (MZ/DZ) to proceed
MIN_SAMPLE_SIZE <- 20

# Print warnings for variables with > X% missing data
MISSING_DATA_THRESHOLD <- 20  # percent

# ---- Scaling Options ----

# Scale variables by grand mean/SD across both twins?
# If TRUE: uses combined stats from T1 and T2
# If FALSE: would require user to pre-scale data
GRAND_SCALE <- TRUE

# ============================================================================
# Do not edit below this line unless you understand the implications
# ============================================================================

# Create output directories if they don't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
if (!dir.exists(LOG_DIR)) {
  dir.create(LOG_DIR, recursive = TRUE)
}
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# Verify data file exists
if (!file.exists(DATA_PATH)) {
  stop(sprintf("Data file not found: %s\nPlease check DATA_PATH in config.R", DATA_PATH))
}

cat("Configuration loaded successfully.\n")
cat(sprintf("Data source: %s\n", DATA_PATH))
cat(sprintf("Number of phenotypes: %d\n", length(PHENOTYPES)))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
