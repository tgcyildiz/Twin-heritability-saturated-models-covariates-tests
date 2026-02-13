# Quick Start Guide

## Installation

1. **Clone this repository**:
   ```bash
   git clone https://github.com/yourusername/twin-heritability-analysis.git
   cd twin-heritability-analysis
   ```

2. **Install R dependencies**:
   ```r
   install.packages(c("OpenMx", "psych", "polycor", "mets"))
   ```

## Setup Your Data

1. Place your twin data CSV file in the `data/` directory
2. Ensure your data matches the expected format (see `data/README_DATA.txt`)

## Configure the Analysis

Edit `config.R` to specify:
- Path to your data file
- Phenotypes to analyze
- Output directory
- Model parameters (optional)

Example:
```r
DATA_PATH <- "data/your_twin_data.csv"
PHENOTYPES <- c("Total_Cerebel_Vol", "Left_VIIB")
OUTPUT_DIR <- "output"
```

## Run the Analysis

From R or RStudio:
```r
source("main_analysis.R")
```

Or from command line:
```bash
Rscript main_analysis.R
```

## Check Results

Results will be saved in the `output/` directory:
- `output/logs/` - Detailed text output for each phenotype
- `output/results/` - CSV tables with model comparisons

## Interpreting Results

### Model Comparison Table
- Compare `-2LL` values between models
- Lower AIC/BIC indicates better fit
- p-value < 0.05: reject the constrained model

### Covariate Tests
- Significant p-value: covariate improves model fit
- Beta estimates: direction and magnitude of effects

### Parameter Estimates
- `vMZ1, vMZ2`: Variances for MZ twins
- `vDZ1, vDZ2`: Variances for DZ twins
- `cMZ21, cDZ21`: Covariances (twin correlations)
- Higher MZ than DZ correlation suggests genetic influence

## Troubleshooting

**Model won't converge:**
- Check data for outliers
- Increase `MAX_ATTEMPTS` in config.R
- Set `USE_MXTRYHARD <- TRUE`

**Missing data errors:**
- Ensure variable names match between data and config
- Check for sufficient complete pairs

**Memory issues:**
- Analyze phenotypes in batches
- Reduce number of phenotypes in `PHENOTYPES` vector

## Getting Help

- See full documentation in README.md
- Check OpenMx documentation: https://openmx.ssri.psu.edu/
- Report issues on GitHub

## Citation

If you use this code, please cite:
- This repository
- OpenMx: Boker et al. (2011) Psychometrika
- Relevant methodology papers for your domain
