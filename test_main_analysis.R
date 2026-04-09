# ============================================================================
# Quick Test Runner - Single phenotype without overwriting main outputs
# ============================================================================

# Tell main_analysis.R to use test config
Sys.setenv(HERIT_CONFIG = "config_test.R")

# Run main pipeline
source("main_analysis.R")
