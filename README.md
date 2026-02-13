# PopPK Analysis

Population pharmacokinetic analysis application built with R Shiny. Supports non-compartmental analysis (NCA) and compartmental modeling (1-, 2-, 3-compartment) using non-linear mixed effects (NLME) estimation.

## Features

- **Data import**: CSV and Excel (XLS/XLSX) files with flexible column name matching
- **Non-Compartmental Analysis**: AUC (linear-log trapezoidal), Cmax, Tmax, half-life, clearance, volume of distribution, MRT
- **Compartmental Modeling**: 1-, 2-, and 3-compartment models via `nlme` with between-subject random effects
- **Routes**: IV bolus, IV infusion, oral, SC, IM, IP, and other extravascular routes
- **Species**: Human, dog, cat, horse, cattle, rat, mouse, rabbit, pig, sheep, goat, NHP, guinea pig, ferret
- **Visualization**: Mean±SEM concentration-time plots (linear and semi-log), individual profiles, model fit overlays, goodness-of-fit diagnostics
- **Export**: CSV parameter tables and PNG plot downloads

## Requirements

- R >= 4.0
- Packages: `shiny`, `DT`, `ggplot2`, `nlme`, `readxl`, `scales`

## Setup

```r
# Install dependencies
source("install_packages.R")

# Run locally
shiny::runApp()
```

## Deploy to shinyapps.io

```r
# Install rsconnect if needed
install.packages("rsconnect")

# Configure your account (one-time)
rsconnect::setAccountInfo(
  name = "your-account-name",
  token = "your-token",
  secret = "your-secret"
)

# Deploy
rsconnect::deployApp()
```

## Input Data Format

Upload a CSV or Excel file in **wide format**: the first column is Time, and each subsequent column contains the concentrations for one animal/subject. The column headers become the subject IDs.

Example:

| Time | Animal1 | Animal2 | Animal3 |
|------|---------|---------|---------|
| 0    | 0       | 0       | 0       |
| 0.5  | 14.5    | 16.1    | 13.2    |
| 1    | 22.2    | 24.7    | 19.8    |
| 2    | 18.7    | 20.1    | 16.9    |
| ...  | ...     | ...     | ...     |

The Time column is detected automatically (accepts headers like `Time`, `Time_h`, `Timepoint`, etc.). If no match is found, the first column is used as Time.

Sample datasets are included in the `data/` folder.

## Analysis Methods

### Non-Compartmental Analysis (NCA)
- AUC calculated using the linear-log trapezoidal rule
- Terminal half-life estimated via log-linear regression on the best-fit terminal phase
- Clearance, volume of distribution, and MRT derived from AUC and lambda_z
- Per-subject parameters with population summary statistics (Mean, SD, SEM, Median, Min, Max)

### Compartmental Modeling
- Non-linear mixed effects models fitted using `nlme::nlme()`
- Parameters estimated on the log scale to ensure positivity
- Between-subject variability on CL (or k) and V (or V1) as random effects
- Automated initial parameter estimation from the data
- Model diagnostics: observed vs. predicted, residual plots, Q-Q plot
- AIC/BIC for model comparison
