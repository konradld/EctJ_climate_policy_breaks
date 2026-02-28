# Replication Package for 'Bayesian Indicator-Saturated Regression for Climate Policy Evaluation'

This repository contains all necessary files to reproduce the paper 'Bayesian Indicator-Saturated Regression for Climate Policy Evaluation' by Konrad, Vashold and Crespo Cuaresma

## Repository Structure

### Overview
Important: **First run the `inst/install_bespoke_mombf_3.5.4.R` file** to install the bespoke version of the `mombf` package! 
- All files to replicate the results of the paper are in `scripts`.
- The `output` folder contain `simulation` and `application` where results and all figures can be found.

### Core Functions

- **`R/estimate_bisam_fun.R`**  
  Main estimation function for Bayesian structural break analysis. This is the key function of the repository.

- **`R/pip_window_fun.R`**  
  Calculates Joint Posterior Inclusion Probabilities (PIPs) and selects break points according to a specified threshold.

- **`R/contr_sim_breaks_fun.R`**  
  Data simulation function that generates panel data with structural breaks and outliers of specified sizes and location.

### Installation

1. Clone this repository:
```r
# In your terminal
git clone https://github.com/konradld/EctJ_climate_policy_breaks.git
cd EctJ_climate_policy_breaks
```

2. Install the bespoke `mombf` package:
```r
source("inst/install_bespoke_mombf_3.5.4.R")
```

## Contact

<lucas.konrad@wu.ac.at>

## Lizense
