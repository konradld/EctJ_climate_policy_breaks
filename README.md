# Overview

This repository contains all necessary files to reproduce results in the the paper 'Bayesian Indicator-Saturated Regression' by Konrad, Vashold and Crespo Cuaresma

It is organized into two main components: 

- **Simulation Study**: replicates results in Section 3 of the manuscript as well as Sections S1 and S4 in the Online Appendix. 

- **Transport Emissions and Climate Policies**: replicaes results in Section 4 of the manuscript as well as Section S5 in the Online Appendix


# Data availability

All data sources are publicly available.

The transport emissions and covariate data used for the analysis in Section 4 in are provided by Koch et al (2022a) and are hosted on Zenodo ([doi:10.5281/zenodo.6768562](https://doi.org/10.5281/zenodo.6768563)).

## Statement about Rights

We certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

We certify that the author(s) have documented permission to redistribute and publish the data contained in this replication package. The empirical data are publicly available from Zenodo, and no restricted or confidential data are included.


# Hardware and Software

## Hardware 

For the creation of outputs as they are included in the manuscript and Online Appendix the following setups have been used:

- Lenovo Thinkpad X13 Yoga Gen 2 (used for replication estimation and creation of paper figures):
  - Intel(R) Core(TM) i5-1134G7 @ 2.40Ghz
  - 16GB RAM
  
- High-preformance computing cluster (used for simulation estimation): 
  - Lenovo Thinksystem SR650 V2
  - 2x Intel Xeon Platinum 8358 with 32 cores @2.6Ghz each per computing node; 26 nodes 
  - 8GB RAM used for each computing instance; up to 15TB available

## Software:

- Local environment (used for creation of paper figures):
  - Windows 11 Enterprise 25H2 (x64)
  - R (4.6.0)
    - required packages: `mombf`, `mvtnorm`, `Matrix`, `matrixStats`, `dplyr`, `glmnet`, `RcppEigen`, `extraDistr`, `stringr`, `gets`, `getspanel`, `ggplot2`, `invgamma`
    - All packages except `mombf` can be installed using `install.packages("pkg-name")`.
    - Use bespoke `mombf` version 3.5.4 for exact reproduction of results (see installation instructions below). 

- High-performance computing cluster (used for simulation estimation):
  - Linux 5.15.0.173.generic Kernel
  - Ubuntu 22.04.5 LTS distribution
  - SLURM job scheduler


# Repository Installation

1. Clone this repository:
```r
# In your terminal
git clone https://github.com/konradld/EctJ_climate_policy_breaks.git
cd EctJ_climate_policy_breaks
```

2. Install the modified `mombf` package:
```r
source("inst/install_bespoke_mombf_3.5.4.R")
```


# Directory Structure

## inst

Contains installation files for the modified `mombf` package (see installation instructions above). 


## R 

Contains the core function used to simulate data, estimate the Bayesian ISR model, and summarize model output. These files are sourced when executing code of files included in the `scripts` folder.

- **`contr_sim_breaks_fun.R`**: Data simulation function that generates panel data with structural breaks and outliers of specified sizes and location.

- **`estimate_bisam_fun.R`**: Main estimation function for Bayesian structural break analysis. This is the key function of the repository.

- **`pip_window_fun.R`**: Calculates Joint Posterior Inclusion Probabilities (PIPs) and selects break points according to a specified threshold.


## data

Contains the dataset on European transport emissions and relevant covariates as provided by Koch et al (2022a).


## scripts

Contains scripts to estimate models for both the simulation study and the replication of the study by Koch et al (2022b). The subfolder `adapted_for_cluster` contains scripts adapted for the simulation study, together with an example `slurm` file for the setup and execution of code at a high-intensity computing cluster.

See the section on Replication Instructions for a more detailed description.


## output

Contains all outputs necessary for the replication of the paper's results. This includes binary files containing the results for the simulation study and the replication of the results by Koch et al (2022) as well as produced figures as presented in the main manuscript and the Online Appendix.

The directory contains subfolders pertaining to and collecting all results and visualizations regarding the simulation and the empirical study in respective subfolders.


# Replication Instructions

The individual files in the repository for the replication of results use relative paths to source other files. The working directory should be set to the root folder of the provided replication package (e.g., with `setwd()`).


## Descriptive Graphs

- **00_imom_pdf_fig_1.R**: Creates Figure 1 of the main manuscript illustrating the probability density functions of the iMom prior for different scale parameters tau.

- **05_prior_rates_fig_S1.R**: Creates Figure S1 of the Online Appendix that illustrates evidence-accumulation rate for different specifications of the slab component for the step-shift design. It creates an intermediary binary `prior_rate_comparison.RDS` in the `output` folder, holding the numerical results.


## Simulation Study

Most estimations for the simulation study are carried out at a high-performance computing cluster that allows for the distribution of parallel tasks across many instances. The folder `scripts/adapted_for_cluster` contains the specific files that were used to create the results stored as binary files as provided in the subfolders of `output/simulation`. This includes `R` files to be executed for estimation of specifications as well as an example `slurm` file that provides the structure to send jobs to the high-computing cluster used. The latter needs to be adapted to the specific cluster setup used together with lines 10-18 in the `R` files.


- **10_simulation_estimation_example.R**: Provides an example file for an individual simulation run. Options that are adjustable (also in other simulation estimation files) include: 
  - Break Environment: Break size (`rel_effect`), break environment (`setup`)
  - Sample Dimensions: Number of cross-sectional units (`Ni`) and time periods (`Nt`)
  - Covariates and fixed effects: Number of external regressors (`Nx`) and which fixed effects to use (`ife` for unit fixed effects, `tfe` for time fixed effect, both for two-way fixed effects)
  - Outlier detection: Whether to include and check for outliers (`do_check_outlier`)
  - BISAM-specific settings: Prior setup for prior inclusion probabilities (`bern` for fixing it to 0.5, `beta-bern` to impose hierarchical Beta-Bernoulli setup)
  - GETS-specific settings:

For all simulation results, the ISR model is estimated with the Bayesian approach proposed in the paper (BISAM), the block-search algorithm GETS, and an adaptive verion of the LASSO (ALASSO). Detection results (exact-date and within a window of \pm 1) are compared to the true breaks as simulated and collected in a summary table. This summary table is saved as binary file for further analyses and visualizations.

To obtain the results necessary for the replication of figures 2, 3, S4, S5, and S6, the following files in the `scripts/adapted_for_cluster` folder need to be run first (potentially adapted to the user-specific setup): 

- **01_simulation_estimation_SD.R**: Creates simulation results using different break sizes (\{0.5, 1, 1.5, 2, 3, 5, 10\}), expressed as multiples of the standard deviation of the variance of the error term, in a sparse and a dense break environment (four units with exactly one break or eight units with breaks where four feature two). Set `date` in line 34 to the current date, followed by suffix `_SD`.

- **02_simulation_estimation_BN.R**: Creates simulation results using a break size of 3, expressed in terms of the standard deviation of the variance of the error term, varying the number of breaks present in the data on the interval \{1,...,20\}. Set `date` in line 34 to the current date, followed by suffix `_BN`.

- **03_simulation_estimation_variations.R**: Creates simulation results for break sizes (\{2, 3, 5\}) in a sparse and a dense break environment, varying different settings in the setup and estimation methods. See the description of the example simulation file above for an overview of the varied settings, results are depicted in Figues S5 and S6 in the Online Appendix. Note that for each specification variation, the corresponding setting in lines 33-45 should be adapted before running the scripts for each setting individually. Set `date` in line 49 to the current date, followed by suffix `_variations`.

- **simulation_cluster_setup.slurm**: Example file for the setup of running the other files in the directory on a computing cluster distributing the tasks on individual nodes.

Set `DATE` to the date that the simulations were carried out and the results are stored in the `output/simulation` folder, and `FIGURE` to the according figure.


The following files then create the figures pertaining to the simulation study as presented in the manuscript and the Online Appendix:

- **11_simulation_figs_2-3-S4.R**: Creates figures 2 and 3 of the manuscript as well as figure S4 of the Online Appendix using results stored in the directory `output/simulation`. For figures 2 and 3 the code uses results stored for a certain date in the subdirectory with suffix `_SD`, for figure S4 it uses the subdirectory with suffix `_BN`. For this code to be executable, simulation results have to be obtained first, see descriptions above.

- **12_simulation_figs_S2-S3.R**: Creates figures S2 and S3 of the Online Appendix, showcasing the sparse and the dense break environments for simulations.

- **13_simulation_figs_S5-S6.R**: Creates figures S5 and S6 of the Online Appendix using results stored in the directory `output/simulation`, showing the performance of BISAM relative to GETS for various additional settings (see description of adjustable options above). Reads in results stored for a certain date in the subdirectory with suffix `_variations`.


## Empirical Study

For creating either Figure 4 in the main manuscript or Figure S7 in the Online Appendix, the following scripts have to be run sequentially, with `FIGURE` set respectively in the code:

- **21_emissions_estimation.R**: Estimates the ISR model using both the Bayesian approach proposed in the paper and the GETS procedure as described in Pretis and Schwarz (2026). With `FIGURE = 4`, a fixed prior on the prior inclusion probabilities of 0.5 is imposed, with `FIGURE = S7` a hierarchical Beta-Bernoulli prior is used.

- **22_emissions_figs_4-S7.R**: Creates the plots presented as Figure 4 in the main manuscript or Figure S7 in the Online Appendix based on estimation results from the previous file.


## Expected Running Times

The running times below were recorded using the setup for hard- and software specified above.

For files run on the local setup:

| Program                             | Expected Running Time (seconds) |
|-------------------------------------|---------------------------------|
| 00_imom_pdf_fig_1.R                 | <1                              |
| 05_prior_rate_fig_S1.R              | 15                              |
| 10_simulation_estimation_example.R  | 168                             |
| 11_simulation_figs_2-3_S4.R         | 18 for Figs 2/3, 38 for Fig S4  |
| 12_simulation_figs_S2-S3.R          | 408 for each figure             |
| 13_simulation_figs_S5-S6.R          | 80 for each figure              |
| 20_emissions_estimation.R           | 2500 for Fig 4, 3000 for Fig S7 |
| 21_emissions_figs_4-S7.R            | 5 for each figure               |


For files run on a high-performance computing cluster (approximate times denoted per replication run): 

| Program                                | Expected Running Time (seconds) per Run |
|----------------------------------------|-----------------------------------------|
| 01_simulation_estimation_SD.R          | 200                                     |
| 02_simulation_estimation_BN.R          | 200                                     |
| 03_simulation_estimation_variations.R  | 200-3000 depending on specification     |

Total expected running times are highly dependent on the possibilities for distributing tasks among individual node instances.


# Citations

Koch, N., Naumann, L., Pretis, F., Ritter, N., & Schwarz, M. (2022a). Attributing agnostically-detected large reductions in road CO2 emissions to policy mixes (code and data) (Version v1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6768563.

Koch, N., Naumann, L., Pretis, F., Ritter, N., & Schwarz, M. (2022b). Attributing agnostically-detected large reductions in road CO2 emissions to policy mixes. Nature Energy, 7, 844–853. http://dx.doi.org/10.1038/s41560-022-01095-6

Pretis, F., and Schwarz, M. (2026). Discovering What Mattered: Detecting Unknown Treatment as Breaks in Panel Models. SSRN. http://dx.doi.org/10.2139/ssrn.4022745 


# News & Contact

Updates about and the newest version of the methodology introduced can be found in the [bisam](https://github.com/konradld/bisam) (Bayesian Indicator Saturated Models) repository.

For questions or comments, please email: [lucas.konrad@wu.ac.at](mailto:lucas.konrad@wu.ac.at) or [lukas.vashold@wu.ac.at](mailto:lukas.vashold@wu.ac.at).

