################################################################################
# Simulation Setup Plots
#
# Description: Creates plots of representative simulation setups
################################################################################

rm(list = ls())


# ==============================================================================
#                         USER CONFIGURATION
# ==============================================================================

FIGURE <- "S2"          # S2 (sparse), S3 (dense)


# ==============================================================================
# SETUP
# ==============================================================================

library(mombf)
set.seed(12345)
setup_type <- if(FIGURE == "S2") "sparse" else "dense"

# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================

# Data dimensions
Ni <- 10          # number of sim. observations
Nt <- 30          # number of sim. time periods
NX <- 0           # number of regressors

# Model structure
DO_CONST <- TRUE    # inclusion of a constant
DO_INDIV_FE <- FALSE     # inclusion of indiv. fixed effects
DO_TIME_FE <- FALSE     # inclusion of time fixed effects
DO_OUTLIERS <- FALSE      # inclusion of indicator saturation
DO_STEP_SATURATION <- TRUE      # inclusion of stepshift saturation

# Error distribution
ERROR_SD <- 1   # standard deviation of the error

# Stepshift characteristics
STEP_MEAN_REL <- 3      # relative mean of size of stepshift in error.sd

# Outlier positions
POS_OUTL <- 0
OUTL_MEAN <- 0

# Break positions
if(setup_type == "sparse") {
  POS_STEP <- c(43, 108, 169, 221)
} else if (setup_type == "dense") {
  POS_STEP <- c(7, 22, 43, 80, 108, 115, 127, 144, 169, 190, 200, 221)
} else {
  stop("Simulation setup not available.")
}
POS_STEP_IN_Z <- POS_STEP - 2 * (POS_STEP %/% Nt + 1) - (POS_STEP %/% Nt)
STEP_MEAN_ABS <- STEP_MEAN_REL * ERROR_SD
S2_TRUE <- ERROR_SD^2

# ==============================================================================
# DATA SIMULATION
# ==============================================================================

source("./R/contr_sim_breaks_fun.R")

sim <- contr_sim_breaks(
  n = Ni, 
  t = Nt, 
  nx = NX, 
  iis = DO_OUTLIERS, 
  sis = DO_STEP_SATURATION,
  const = DO_CONST, 
  ife = DO_INDIV_FE, 
  tfe = DO_TIME_FE,
  pos.outl = POS_OUTL, 
  pos.step = POS_STEP,
  outl.mean = OUTL_MEAN, 
  step.mean = STEP_MEAN_ABS,
  error.sd = ERROR_SD
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

make_Z <- function(n, t, i_index = 1, t_index = 2) {
  z <- lower.tri(matrix(1, nrow = t, ncol = t), diag = TRUE)[, -c(1:2, t)]
  mode(z) <- "integer"
  Z <- kronecker(diag(n), z)
  mode(Z) <- "integer"
  
  sis_grid <- expand.grid(unique(data[, t_index])[-c(1:2, t)], 
                          unique(data[, i_index]))
  iis_grid <- expand.grid(unique(data[, t_index]), 
                          unique(data[, i_index]))
  colnames(Z) <- paste0("sis.", sis_grid[, "Var2"], ".", sis_grid[, "Var1"])
  
  return(Z)
}

# ==============================================================================
# TRUE MODEL FIT
# ==============================================================================

data <- sim$data
y <- data[, 3]
Z <- make_Z(Ni, Nt)
true_fit <- lm(y ~ Z[, POS_STEP_IN_Z])

df <- cbind(as.data.frame(data), as.data.frame(Z[, POS_STEP_IN_Z]))
names(df) <- c("n","t","y", paste0("x",1:length(POS_STEP)))
if (DO_INDIV_FE & DO_TIME_FE) {
  true_fit <- fixest::feols(as.formula(paste0("y~", paste0("x", 1:length(POS_STEP), collapse = "+"), "| n + t")) , data = df)
} else if (DO_INDIV_FE) {
  true_fit <- fixest::feols(as.formula(paste0("y~", paste0("x", 1:length(POS_STEP), collapse = "+"), "| n")) , data = df)  
} else if (DO_TIME_FE) {
  true_fit <- fixest::feols(as.formula(paste0("y~", paste0("x", 1:length(POS_STEP), collapse = "+"), "| t")), data = df)
} else {
  true_fit <- fixest::feols(as.formula(paste0("y~", paste0("x", 1:length(POS_STEP), collapse = "+"))), data = df)
}

# ==============================================================================
# MODEL ESTIMATION PARAMETERS
# ==============================================================================

# Data processing
I_INDEX <- 1
T_INDEX <- 2
Y_INDEX <- 3
DO_CENTER_Y <- FALSE
DO_SCALE_Y <- FALSE
DO_CENTER_X <- FALSE
DO_SCALE_X <- FALSE

# MCMC settings
NDRAW <- 10000L
NBURN <- 2000L

# Prior specification
BETA_PRIOR <- "f"
STEP_SIZE_PRIOR <- "imom"
STEP_INCL_PRIOR <- "bern"

# Prior settings
BETA_VARIANCE_SCALE <- 100

SIGMA2_SHAPE <- NULL
SIGMA2_RATE <- NULL
SIGMA2_HYPER_P <- 0.90

if (STEP_SIZE_PRIOR == "imom") {
  STEP_SIZE_SCALE <- priorp2g(0.01, 1, nu = 1, prior = "iMom")
} else if (STEP_SIZE_PRIOR == "mom") {
  STEP_SIZE_SCALE <- priorp2g(0.01, 1, nu = 1, prior = "normalMom")
} else {
  stop("selected prior not implemented")
}

STEP_INCL_PROB <- 0.5

# Variance settings
DO_CLUSTER_S2 <- FALSE
DO_CHECK_OUTLIER <- FALSE
DO_SV <- FALSE

# ==============================================================================
# RUN MODEL
# ==============================================================================

source("./R/estimate_bisam_fun.R")

mod <- estimate_bisam(
  data = data,
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
  step_size_scale = STEP_SIZE_SCALE,
  do_cluster_s2 = DO_CLUSTER_S2,
  do_check_outlier = DO_CHECK_OUTLIER,
  do_sv = DO_SV
)


# ==============================================================================
# PLOTTING
# ==============================================================================

# Extract sample dimensions
t_mod <- mod$meta$sample['t']
n_mod <- mod$meta$sample['n']
N_mod <- mod$meta$sample['N']

# Color palette
COL_MAIN <- "#2166AC"    # Blue
COL_STEP <- "#B2182B"    # Red
COL_FIT  <- "#1B9E77"    # Teal
COL_IIS  <- "#762A83"    # Purple
COL_GRID <- "gray70"

pdf(sprintf("./output/simulation/figure_%s.pdf", FIGURE),
    width = 16, height = 9)

par(
  mfrow = c(2, 1),
  mar  = c(1, 4.5, 1.5, 1),
  oma  = c(1, 0, 2, 0),
  cex.axis = 1,
  cex.lab  = 1.25,
  las = 1
)

# ---- Plot 1: Omega coefficients ----

omega_with_breaks <- mod$coefs$omega

omega_plot <- rep(NA, length(data[, 3]))
current_pos <- 1
for(i in 1:n_mod) {
  start_idx <- (i-1) * (t_mod) + 3
  end_idx   <- i * (t_mod) - 1
  segment_length <- end_idx - start_idx + 1
  
  omega_plot[start_idx:end_idx] <- omega_with_breaks[current_pos:(current_pos + segment_length - 1)]
  current_pos <- current_pos + segment_length
}

x_coords_a <- 1:length(omega_plot)

plot(
  x_coords_a, omega_plot,
  type = "l",
  col  = COL_MAIN,
  lwd  = 2,
  ylim = c(0, 1),
  xlab = "",
  ylab = "PIP",
  main = "",
  xaxt = "n",
  frame.plot = FALSE
)

abline(h = c(0, 0.25, 0.5, 0.75, 1), lty = 3, col = COL_GRID, lwd = 0.8)
abline(v = 0:Ni * Nt,  lty = 2, col = COL_GRID, lwd = 0.8)
abline(v = POS_STEP,   lty = 2, col = COL_STEP,  lwd = 1)

legend(
  "topright",
  legend = c("PIP", "True Breaks"),
  col    = c(COL_MAIN, COL_STEP),
  lty    = c(1, 2),
  lwd    = c(2, 1.5),
  bty    = "n",
  cex    = 1.25
)

mtext(side = 3, line = 0.5, text = "Panel A: PIP estimates",
      adj = 0, font = 2, cex = 1.5)

midpoints_a <- sapply(1:Ni, function(i) {
  mean(c((i-1)*Nt + 1, i*Nt)) - .5
})
text(x = midpoints_a, y = 0,
     labels = paste0("Unit ", seq_len(Ni)),
     pos = 1, xpd = TRUE, cex = 1.25)

# ---- Plot 2: Response variable ----
par(mar = c(3, 4.5, 0.5, 1))

plot(
  data[, 3], 
  cex = 0.6,
  pch = 19, 
  col = adjustcolor(COL_MAIN, alpha.f = 0.4),
  xlab = "",
  ylab = "Response (y)",
  main = "",
  xaxt = "n",
  frame.plot = FALSE
)

# Reference lines
abline(v = 0:Ni * Nt, lty = 2, col = COL_GRID, lwd = 0.8)
abline(v = POS_STEP, col = COL_STEP, lty = 2, lwd = 1.5)
# abline(h = sim$true.const, lty = 3, col = COL_GRID, lwd = 0.8)
# abline(h = sim$true.const + STEP_MEAN_ABS, lty = 3, col = COL_GRID, lwd = 0.8)

# Fitted lines with breaks every 30th value
n_fit <- length(true_fit$fitted.values)

for (start in seq(1, n_fit, by = 30)) {
  end <- min(start + 29, n_fit)
  lines(start:end, true_fit$fitted.values[start:end], col = scales::alpha(COL_STEP, 0.75), lwd = 2)
  lines(start:end, mod$fitted[start:end],            col = scales::alpha(COL_FIT, 0.75), lwd = 2)
}

# Legend
legend(
  "topright",
  legend = c("Observed", "True Mean", "BISAM-fit", "True Breaks"),
  col = c(
    adjustcolor(COL_MAIN, alpha.f = 0.9),
    COL_STEP,
    COL_FIT, 
    COL_STEP
  ),
  pch = c(19, NA, NA, NA),
  lty = c(NA, 1, 1, 2),
  lwd = c(NA, 2, 2, 1.5),
  bty = "n",
  cex = 1.25,
  bg = "white"
)

mtext(
  side = 3, 
  line = -0.3, 
  text = "Panel B: Fitted values of y", 
  adj = 0, 
  font = 2, 
  cex = 1.5
)

# Add unit labels
midpoints_a <- sapply(1:Ni, function(i) {
  start_idx <- (i-1) * Nt + 1
  end_idx <- i * (Nt)
  mean(c(start_idx, end_idx)) - .5
})
text(
  x = midpoints_a,
  y = min(data[, 3], na.rm = TRUE),
  labels = paste0("Unit ", seq_len(Ni)),
  pos = 1,
  xpd = TRUE,
  cex = 1.25
)

# Overall title
title(
  main = bquote(
    "Step size: " ~  
      .(round(STEP_MEAN_REL, 2)) ~ sigma ~ 
      # " / " ~ tau ~ "=" ~ .(round(TAU, 2)) ~ 
      " / " ~ hat(sigma) ~ "=" ~ .(round(sqrt(mod$coefs$sigma2), 2))
  ),
  outer = TRUE, 
  cex.main = 2, 
  font.main = 2
)

dev.off()