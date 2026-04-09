# ============================================================================
# Test Configuration - Single Phenotype Run (No overwrite of main outputs)
# ============================================================================

# ---- Data Configuration ----
data_candidates <- c(
	"data/twin_data.csv"
)

existing_candidates <- data_candidates[file.exists(data_candidates)]
if (length(existing_candidates) == 0) {
	stop("No valid DATA_PATH found for test run. Update data_candidates in config_test.R")
}
DATA_PATH <- existing_candidates[1]

TWIN_ID_COL <- "tvparnr"
ZYGOSITY_COL <- "zyg1"

AGE_COLS <- c("age1", "age2")
SEX_COLS <- c("sex1", "sex2")
ETIV_COLS <- c("eTIV1", "eTIV2")

# Restrict to one phenotype for quick testing
PHENOTYPES <- c("Total_Cerebel_Vol")

# ---- Output Configuration ----
OUTPUT_DIR <- file.path(
	"output/test_runs",
	format(Sys.time(), "%Y%m%d_%H%M%S")
)
LOG_DIR <- file.path(OUTPUT_DIR, "logs")
RESULTS_DIR <- file.path(OUTPUT_DIR, "results")

# ---- Model Configuration ----
SV_MEANS <- 0
SV_VARIANCE <- 1.0
LB_VARIANCE <- 0.0001
LB_CORRELATION <- -0.99

SV_BETA_AGE <- -0.1
SV_BETA_SEX <- 0
SV_BETA_ETIV <- 0.2

# ---- Optimization Configuration ----
USE_MXTRYHARD <- TRUE
MAX_ATTEMPTS <- 5
COMPUTE_CI <- TRUE

# ---- Data Quality Checks ----
PRINT_DESCRIPTIVES <- TRUE
MIN_SAMPLE_SIZE <- 20
MISSING_DATA_THRESHOLD <- 20

# ---- Scaling ----
GRAND_SCALE <- TRUE

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

cat("Test configuration loaded.\n")
cat(sprintf("Data source: %s\n", DATA_PATH))
cat(sprintf("Phenotypes: %s\n", paste(PHENOTYPES, collapse = ", ")))
cat(sprintf("Test output directory: %s\n", OUTPUT_DIR))
