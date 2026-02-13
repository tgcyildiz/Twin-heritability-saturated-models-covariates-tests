# Twin Study Heritability Analysis with OpenMx

A comprehensive R pipeline for analyzing heritability of quantitative traits in twin studies using saturated and nested models in OpenMx. This analysis includes covariate adjustments (age, sex, total intracranial volume) and provides confidence intervals for genetic and environmental correlations.

## Overview

This pipeline implements:
- **Saturated Model**: Baseline model with all parameters free
- **Nested Models**: Progressive constraints on means, variances, and correlations
- **Covariate Adjustment**: Regression effects of age, sex, and eTIV
- **Model Comparison**: Likelihood ratio tests and effect sizes
- **Confidence Intervals**: Bootstrap CI estimation for parameters

## Project Structure

```
├── README.md                 # This file
├── config.R                  # User configuration (paths, variables)
├── main_analysis.R           # Main analysis pipeline
├── helper_functions.R        # Custom functions (estimate extraction, model fitting)
├── data/                     # Data directory (user-provided)
│   └── twin_data.csv         # Twin study data (ID, age, sex, eTIV, phenotypes, zygosity)
└── output/                   # Results directory (auto-created)
    ├── logs/                 # Text output files
    ├── results/              # CSV tables
    └── figures/              # Optional plots
```

## Data Format

Expected CSV file with columns:
- `tvparnr`: Twin pair ID
- `zyg1`: Zygosity (1=MZ, 2=DZ)
- `age1`, `age2`: Age at time point 1 and 2
- `sex1`, `sex2`: Sex for each twin
- `eTIV1`, `eTIV2`: Total intracranial volume
- Phenotype columns: Volume measurements (matched pairs with suffixes 1, 2)
  - e.g., `Total_Cerebel_Vol1`, `Total_Cerebel_Vol2`

## Usage

### Quick Start

1. **Set configuration** in `config.R`:
   ```r
   data_path <- "path/to/your/data.csv"
   phenotypes <- c("Total_Cerebel_Vol", "Left_VIIB", ...)
   output_dir <- "path/to/output"
   ```

2. **Run analysis**:
   ```r
   source("main_analysis.R")
   ```

3. **Outputs**:
   - `{phenotype}std.txt` - Detailed model results and goodness-of-fit
   - `{phenotype}_testsstd.csv` - Model comparison table
   - `{phenotype}_covstd.csv` - Covariate effects model comparison

### Custom Configuration

Edit `config.R` to:
- Specify data file path
- Select phenotypes to analyze
- Change output directory
- Adjust model starting values
- Modify convergence parameters

## References

### Key Functions (via OpenMx)

- `mxModel()`: Model specification
- `mxRun()`: Model optimization
- `mxTryHard()`: Robust optimization with multiple starting values
- `mxCompare()`: Nested model comparison
- `mxCI()`: Confidence interval estimation

### Methods

Twin models assume:
- MZ correlation = 1.0 (genetic effects)
- DZ correlation = 0.5 (genetic effects)
- Equal environments assumption

Saturated models estimate all means, variances, and covariances without constraints. Nested models test specific hypotheses by progressively constraining parameters.

## Dependencies

```r
library(OpenMx)      # Structural equation modeling
library(psych)       # Descriptive statistics
library(polycor)     # Polychoric correlations (optional)
library(mets)        # Twin family methods
```

Install with:
```r
install.packages(c("OpenMx", "psych", "polycor", "mets"))
```

## Output Interpretation

### Model Comparison Table
- **-2LL**: Minus 2 log-likelihood (lower = better fit)
- **df**: Degrees of freedom
- **AIC/BIC**: Information criteria
- **p**: P-value from likelihood ratio test

### Covariate Table
- **Model**: Tests removing each covariate from full model
- **p < 0.05**: Covariate significantly improves fit

### Parameter Estimates
- **vMZ1, vMZ2**: MZ variances
- **vDZ1, vDZ2**: DZ variances  
- **cMZ21, cDZ21**: Covariances (genetic correlation)
- **bA, bS, bE**: Age, sex, eTIV regression coefficients

## Troubleshooting

**Model fails to converge:**
- Check data for outliers or missing patterns
- Adjust starting values in `config.R`
- Increase tolerance/iteration settings

**Negative variances:**
- Scale variables appropriately
- Check data quality

**Singular covariance matrix:**
- Ensure sufficient sample size
- Check for multicollinearity in covariates

## License

[Specify your license]

## Contact

[Your information]

## Citation

Please cite this repository and the following references:

- Neale, M. C., & Cardon, L. R. (1992). *Methodology for genetic studies of twins and families*. Kluwer Academic.
- Boker, S., Neale, M., Maes, H., et al. (2011). OpenMx: An open source extended structural equation modeling framework. *Psychometrika*, 76(2), 306.
