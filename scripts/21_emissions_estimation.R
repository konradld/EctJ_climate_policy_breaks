################################################################################
# Plotting of replication results for Koch et al (2022)
################################################################################

rm(list = ls())

#===============================================================================
#           Which Figure to replicate (see paper for numbering)
              FIGURE <- "4" # 4 for figure in main test, S7 for appendix
              DATE <- "2026-07-21"
#===============================================================================


# ==============================================================================
# SETUP
# ==============================================================================

library(dplyr)
library(stringr)
library(gets)
library(getspanel)
library(Matrix)
library(mombf)

beta_prior <- "f"
beta_scale <- 10
check_outl <- TRUE
outl_scale <- 10
sis_prior <- "imom"
tau <- mombf::priorp2g(0.05, q = 1, prior = "iMom")
incl_prior <- ifelse(FIGURE == "4", "bern", "beta_bern")

# ==============================================================================
# CREATE SAVING DIRECTORY
# ==============================================================================

dir_res <- sprintf("./output/emissions/%s/", DATE)
dir.create(dir_res, showWarnings = FALSE)


# Load and prepare data ---------------------------------------------------

dir_data <- "./data/CO2DriversEU_dataset_CLEAN.csv"
data <- read.csv(dir_data)[-1]

# Group specification
EU15   <- c("Austria", "Belgium", "Germany", "Denmark", "Spain", "Finland",
            "France", "United Kingdom", "Ireland", "Italy", "Luxembourg", 
            "Netherlands", "Greece", "Portugal", "Sweden")
data_    <- data[, c('country','year','ltransport.emissions','lgdp','lgdp_sq','lpop')]
dat      <- filter(data_, country %in% EU15, year>=1995)

# i_names <- unique(dat$country)
# t_names <- unique(dat$year)
# 
# n <- length(i_names)
# t <- length(t_names)

# ==============================================================================
# ESTIMATE BISAM
# ==============================================================================

source("./R/estimate_bisam_fun.R")

# Data processing
I_INDEX <- 1
T_INDEX <- 2
Y_INDEX <- 3
DO_CONST <- FALSE
DO_INDIV_FE <- TRUE
DO_TIME_FE <- TRUE 
DO_CENTER_Y <- FALSE
DO_SCALE_Y <- FALSE
DO_CENTER_X <- FALSE
DO_SCALE_X <- FALSE

# MCMC settings
NDRAW <- 50000L
NBURN <- 25000L

# Prior settings
BETA_VARIANCE_SCALE <- beta_scale
BETA_PRIOR <- beta_prior

STEP_SIZE_PRIOR <- sis_prior
TAU <- tau

STEP_INCL_PRIOR <- incl_prior
STEP_INCL_PROB <- 0.5
STEP_INCL_ALPHA <- 1
STEP_INCL_BETA <- 1

DO_CHECK_OUTLIER <- check_outl
OUTLIER_SCALE <- outl_scale
OUTLIER_INCL_ALPHA <- 1
OUTLIER_INCL_BETA <- 10

DO_CLUSTER_S2 <- TRUE
SIGMA2_SHAPE <- NULL
SIGMA2_RATE <- NULL
SIGMA2_HYPER_P <- 0.9
DO_SV <- FALSE

ssvs_i <- estimate_bisam(
  data = dat,
  do_constant = DO_CONST,
  do_individual_fe = DO_INDIV_FE,
  do_time_fe = DO_TIME_FE,
  y_index = Y_INDEX,
  i_index = I_INDEX,
  t_index = T_INDEX,
  do_center_y = DO_CENTER_Y,
  do_scale_y = DO_SCALE_Y,
  do_center_x = DO_CENTER_X,
  do_scale_x = DO_SCALE_X,
  Ndraw = NDRAW,
  Nburn = NBURN,
  
  beta_prior = BETA_PRIOR,
  step_size_prior = STEP_SIZE_PRIOR,
  step_incl_prior = STEP_INCL_PRIOR,
  beta_variance_scale = BETA_VARIANCE_SCALE,
  sigma2_shape = SIGMA2_SHAPE,
  sigma2_rate = SIGMA2_RATE,
  sigma2_hyper_p = SIGMA2_HYPER_P,
  step_incl_prob = STEP_INCL_PROB,
  step_incl_alpha = STEP_INCL_ALPHA,
  step_incl_beta = STEP_INCL_BETA,
  step_size_scale = TAU,
  do_cluster_s2 = DO_CLUSTER_S2,
  do_check_outlier = DO_CHECK_OUTLIER,
  outlier_incl_alpha = OUTLIER_INCL_ALPHA,
  outlier_incl_beta = OUTLIER_INCL_BETA,
  outlier_scale = OUTLIER_SCALE,
  do_sv = DO_SV
)

dir_save <- sprintf(paste0(dir_res, 
                           "ssvs_outlier-%s-%s_beta-%s-%s_sisprior-%s-%s_inclprior-%s.RDS"), 
                    check_outl,
                    outl_scale,
                    beta_prior, 
                    beta_scale,
                    sis_prior,
                    tau,
                    incl_prior)  
saveRDS(ssvs_i, dir_save)

# ==============================================================================
# ESTIMATE GETS
# ==============================================================================

formula <- "ltransport.emissions ~ lgdp + lgdp_sq + lpop"
index   <- c("country", "year")

p_value <- 0.05

gets_i <- isatpanel(
  data    = dat,
  formula = as.formula(formula),
  index   = index,
  effect  = "twoways",
  iis     = TRUE,
  jsis    = FALSE,
  fesis   = TRUE, 
  t.pval  = p_value,
  print.searchinfo = FALSE)

dir_save <- sprintf(paste0(dir_res, "gets_%s.RDS"), 
                    p_value)
saveRDS(gets_i, dir_save)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================