################################################################################
# Comparison of 
#
# Description: 
#
# Figure S1: vary sparse vs dense setting (SD breakpoints), 
#
# Output: Comparison file, plot
################################################################################

rm(list = ls())

# ==============================================================================
#                         USER CONFIGURATION
# ==============================================================================

DATE   <- "2026-07-24"

# ==============================================================================
# 1. SETUP
# ==============================================================================

library(mombf)
set.seed(20260716)

TS       <- c(20, 40, 80, 160, 320, 640)
REPS     <- 200
STEP_SD  <- 3      # true break size in error SD
NI       <- 10     # units, for the multiplicity burden log K = log(Ni*(T-3))

#   iMOM : priorp2g(0.01, 1, nu=1, prior="iMom")      = 3.3174
#   pMOM : priorp2g(0.01, 1, nu=1, prior="normalMom") = 8.7084
#   normal slab N(0, tau*sigma^2): P(|g|<=sigma) = 2*pnorm(1/sqrt(tau)) - 1 = 0.01
#          => tau = (1/qnorm(0.505))^2.
TAU_IMOM <- priorp2g(0.01, 1, nu = 1, prior = "iMom")
TAU_MOM  <- priorp2g(0.01, 1, nu = 1, prior = "normalMom")
TAU_NORM <- (1 / qnorm(0.505))^2

prior_arms <- list(
  list(key = "zellner_g_n",   fam = "local", lab = "Local: Zellner",
       f = function(n) zellnerprior(tau = n), 
       col = "#E69F00", lty = 1),
  list(key = "normalid_cal",  fam = "local", lab = "Local: normal slab",
       f = function(n) normalidprior(tau = TAU_NORM), 
       col = "#E69F00", lty = 3),
  list(key = "pmom",          fam = "mom",   lab = "pMOM",
       f = function(n) momprior(tau = TAU_MOM), 
       col = "#CC79A7", lty = 2),
  list(key = "pimom",         fam = "imom",  lab = "piMOM",
       f = function(n) imomprior(tau = TAU_IMOM), 
       col = "#0072B2", lty = 1)
)


# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

# create SIS block
sis_block <- function(Ti) {
  1 * lower.tri(matrix(1, Ti, Ti), diag = TRUE)[, 3:(Ti - 1), drop = FALSE]
}

# log BF for including candidate j against the null model, on data y.
log_bf <- function(y, x, pr) {
  m0 <- nlpMarginal(sel = integer(0), y = y, x = x, priorCoef = pr,
                    logscale = TRUE, method = "Laplace")
  m1 <- nlpMarginal(sel = 1, y = y, x = x, priorCoef = pr,
                    logscale = TRUE, method = "Laplace")
  m1 - m0
}

# simulation function
#   H0 arm: pure-noise unit, evidence AGAINST a randomly chosen spurious candidate.
#   H1 arm: a true STEP_SD break at a random date, evidence FOR the true candidate.
run_all <- function() {
  res <- expand.grid(T = TS, arm = c("H0", "H1"), key = sapply(prior_arms, `[[`, "key"),
                     stringsAsFactors = FALSE)
  res$mean <- NA_real_; res$se <- NA_real_
  t0 <- Sys.time()
  for (Ti in TS) {
    L <- sis_block(Ti)
    K <- ncol(L)
    for (a in prior_arms) {
      pr <- a$f(Ti)
      v0 <- v1 <- numeric(REPS)
      for (r in seq_len(REPS)) {
        j <- sample(K, 1)
        # --- H0: null is TRUE at candidate j
        y0 <- as.numeric(scale(rnorm(Ti), scale = FALSE))
        v0[r] <- log_bf(y0, L[, j, drop = FALSE], pr)
        # --- H1: a real break of STEP_SD sigma at candidate j
        y1 <- as.numeric(scale(STEP_SD * L[, j] + rnorm(Ti), scale = FALSE))
        v1[r] <- log_bf(y1, L[, j, drop = FALSE], pr)
      }
      i0 <- res$T == Ti & res$arm == "H0" & res$key == a$key
      i1 <- res$T == Ti & res$arm == "H1" & res$key == a$key
      res$mean[i0] <- mean(v0); res$se[i0] <- sd(v0) / sqrt(REPS)
      res$mean[i1] <- mean(v1); res$se[i1] <- sd(v1) / sqrt(REPS)
    }
    cat(sprintf("  T = %4d done (%.1f min)\n", Ti,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  res
}

# summary table function
rate_table <- function(res) {
  out <- do.call(rbind, lapply(prior_arms, function(a) {
    d <- res[res$key == a$key & res$arm == "H0", ]
    d <- d[order(d$T), ]
    f_log <- lm(d$mean ~ log(d$T))
    f_sqrt <- lm(d$mean ~ sqrt(d$T))
    h <- res[res$key == a$key & res$arm == "H1", ]
    h <- h[order(h$T), ]
    f_lin <- lm(h$mean ~ h$T)
    data.frame(
      prior        = a$lab,
      slope_logT   = coef(f_log)[2],  r2_logT  = summary(f_log)$r.squared,
      slope_sqrtT  = coef(f_sqrt)[2], r2_sqrtT = summary(f_sqrt)$r.squared,
      favours      = ifelse(AIC(f_log) < AIC(f_sqrt), "log T (polynomial)",
                            "sqrt(T) (root-exponential)"),
      H1_slope_T   = coef(f_lin)[2],  H1_r2 = summary(f_lin)$r.squared,
      stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}


# ==============================================================================
# 3. RUN SIMULATIONS
# ==============================================================================

cat(sprintf("tau: iMOM %.4f | pMOM %.4f | normal slab %.1f\n",
            TAU_IMOM, TAU_MOM, TAU_NORM))
cat(sprintf("simulating %d reps x %d T-values x %d priors x 2 arms ...\n",
            REPS, length(TS), length(prior_arms)))
res <- run_all()
tab <- rate_table(res)
cat("\n--- Evidence AGAINST a spurious break (H0 true) ---\n")
print(format(tab[, c("prior", "slope_logT", "r2_logT", "slope_sqrtT",
                     "r2_sqrtT", "favours")], digits = 3), row.names = FALSE)
cat("\n--- Evidence FOR a true break (H1 true): slope in T ---\n")
print(format(tab[, c("prior", "H1_slope_T", "H1_r2")], digits = 3),
      row.names = FALSE)

# save for later use
saveRDS(list(res = res, rates = tab, reps = REPS, step_sd = STEP_SD, Ni = NI,
             tau = c(imom = TAU_IMOM, mom = TAU_MOM, normal = TAU_NORM)),
        "./output/simulation/prior_rate_comparison.RDS")

cat("\nwrote", "./output/simulation/prior_rate_comparison.RDS", "\n")

# ==============================================================================
# 4. PLOT SETTINGS
# ==============================================================================

colors <- list(
  ssvs   = "#0072B2",
  gets   = "#D55E00", 
  alasso = "#009E73",
  
  local  = "#E69F00", 
  mom    = "#CC79A7", 
  imom   = "#0072B2", 

  ink    = "#1A1A1A",
  mute   = "#7A7A7A",
  grid   = "#E5E5E5"
)

settings <- list(
  width      = 6.5,
  height     = 2.55,
  ps         = 9,
  family     = "serif",
  lwd_line   = 1.4,
  lwd_axis   = 0.7,
  lwd_grid   = 0.5,
  lwd_ref    = 0.9,
  cex_pt     = 0.62,
  cex_lab    = 1.0,
  cex_axis   = 0.9,
  cex_leg    = 0.85
)

setup_plot <- function(mfrow = c(1, 1), mar = c(2.9, 3.1, 1.2, 0.7),
                       mgp = c(1.85, 0.55, 0)) {
  par(mfrow = mfrow, mar = mar, mgp = mgp, las = 1,
      cex = 1, cex.lab = settings$cex_lab, cex.axis = settings$cex_axis,
      col.lab = colors$ink, col.axis = colors$ink, col.main = colors$ink,
      family = settings$family, tcl = -0.2, bty = "n", xpd = FALSE)
}

draw_grid <- function(h = NULL, v = NULL) {
  if (!is.null(h)) abline(h = h, col = colors$grid, lwd = settings$lwd_grid, lty = 1)
  if (!is.null(v)) abline(v = v, col = colors$grid, lwd = settings$lwd_grid, lty = 1)
}

add_clean_axes <- function(at_x = NULL, lab_x = TRUE, at_y = NULL, lab_y = TRUE) {
  axis(1, at = at_x, labels = lab_x, col = colors$mute, col.axis = colors$ink,
       lwd = settings$lwd_axis, 
       cex.lab = settings$cex_axis, cex.axis = settings$cex_axis)
  axis(2, at = at_y, labels = lab_y, col = colors$mute, col.axis = colors$ink,
       lwd = settings$lwd_axis, 
       cex.lab = settings$cex_axis, cex.axis = settings$cex_axis)
}

add_panel_label <- function(lab) {
  mtext(sprintf("(%s)", lab), side = 3, line = 0.25, adj = 0,
        cex = settings$cex_lab, font = 1, col = colors$ink)
}

draw_ref <- function(h = NULL, v = NULL) {
  if (!is.null(h)) abline(h = h, col = colors$ink, lty = 2, lwd = settings$lwd_ref)
  if (!is.null(v)) abline(v = v, col = colors$ink, lty = 2, lwd = settings$lwd_ref)
}

draw_series <- function(x, y, col, pch = 16, lty = 1, lwd = settings$lwd_line,
                        pts = TRUE) {
  lines(x, y, col = col, lwd = lwd, lty = lty)
  if (pts) points(x, y, col = col, pch = pch, cex = settings$cex_pt)
}

draw_band <- function(x, lo, hi, col, alpha = 0.15) {
  polygon(c(x, rev(x)), c(lo, rev(hi)), col = adjustcolor(col, alpha),
          border = NA)
}

add_label_end <- function(x, y, txt, col, pos = 4, off = 0.3, cex = settings$cex_leg) {
  n <- length(x)
  text(x[n], y[n], txt, col = col, pos = pos, offset = off, cex = cex, xpd = NA)
}

add_legend <- function(where, legend, col, lty = 1, pch = 16, ...) {
  legend(where, legend = legend, col = col, lty = lty, pch = pch,
         lwd = settings$lwd_line, bty = "n", cex = settings$cex_leg,
         pt.cex = settings$cex_pt, seg.len = 1.8, ...)
}

add_legend_outer <- function(legend, col, lty = 1, pch = 16, ncol = length(legend)) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("top", legend = legend, col = col, lty = lty, pch = pch,
         lwd = settings$lwd_line, bty = "n", cex = settings$cex_leg,
         pt.cex = settings$cex_pt, seg.len = 1.8, horiz = (ncol == length(legend)),
         ncol = ncol, xpd = NA)
}

get <- function(key, arm) {
  d <- res[res$key == key & res$arm == arm, ]
  d[order(d$T), ]
}

# ==============================================================================
# 4. PLOT RESULTS
# ==============================================================================

cairo_pdf("./output/simulation/figure_S1.pdf", 
          width = settings$width, height =  settings$height,
          pointsize = settings$ps, family = settings$family)

par(mfrow = c(1, 3), mar = c(2.9, 3.0, 1.2, 1.6), mgp = c(1.85, 0.55, 0), las = 1,
    cex = 1, cex.lab = settings$cex_lab, cex.axis = settings$cex_axis,
    col.lab = colors$ink, col.axis = colors$ink, col.main = colors$ink,
    family = settings$family, tcl = -0.2, bty = "n", xpd = FALSE)

xt <- c(20, 80, 320, 640)
h0 <- lapply(prior_arms, function(a) get(a$key, "H0"))
lo_a <- min(unlist(lapply(h0, function(d) d$mean)))
yt_a <- seq(-50, 0, 10)
ylim_a <- c(min(lo_a * 1.06, -10), 2)

# ---- (a) evidence AGAINST a spurious break ----------------------------------
plot(NA, xlim = range(TS), ylim = ylim_a, xlab = "Time periods (T)",
     ylab = "log Bayes factor", main = "", axes = FALSE, log = "x")
draw_grid(h = yt_a)
abline(h = 0, col = colors$mute, lwd = colors$lwd_axis)
for (i in seq_along(prior_arms))
  draw_series(h0[[i]]$T, h0[[i]]$mean, prior_arms[[i]]$col, lty = prior_arms[[i]]$lty)
add_clean_axes(at_x = xt, at_y = yt_a)
add_label_end(h0[[1]]$T, h0[[1]]$mean, "local", colors$local, pos = 3, off = 0.15)
add_label_end(h0[[3]]$T, h0[[3]]$mean, "pMOM",  colors$mom,   pos = 3, off = 0.15)
add_label_end(h0[[4]]$T, h0[[4]]$mean, "piMOM", colors$imom,  pos = 3, off = 0.15)
add_panel_label("a")

# ---- (b) evidence FOR a true break ------------------------------------------
h1 <- lapply(prior_arms, function(a) get(a$key, "H1"))
ylim_b <- c(0, max(unlist(lapply(h1, function(d) d$mean))) * 1.08)

plot(NA, xlim = range(TS), ylim = ylim_b, xlab = "Time periods (T)",
     ylab = "log Bayes factor", main = "", axes = FALSE, log = "x")
draw_grid(h = seq(0, 120, 20))
for (i in seq_along(prior_arms))
  draw_series(h1[[i]]$T, h1[[i]]$mean, prior_arms[[i]]$col, lty = prior_arms[[i]]$lty)
add_clean_axes(at_x = xt, at_y = seq(0, 120, 20))
add_legend("topleft",
           legend = sapply(prior_arms, `[[`, "lab"),
           col    = sapply(prior_arms, `[[`, "col"),
           lty    = sapply(prior_arms, `[[`, "lty"),
           pch    = 16)
add_panel_label("b")

# ---- (c) evidence vs the multiplicity burden --------------------------------
logK <- log(NI * (TS - 3))
ev   <- lapply(h0, function(d) -d$mean)
ylim_c <- c(0, max(unlist(ev)) * 1.08)

plot(NA, xlim = range(TS), ylim = ylim_c, xlab = "Time periods (T)",
     ylab = "Evidence vs. burden (log-odds)", main = "", axes = FALSE, log = "x")
draw_grid(h = seq(0, 50, 10))
for (i in seq_along(prior_arms))
  draw_series(TS, ev[[i]], prior_arms[[i]]$col, lty = prior_arms[[i]]$lty)
lines(TS, logK, col = colors$ink, lty = 2, lwd = settings$lwd_ref)
add_clean_axes(at_x = xt, at_y = seq(0, 50, 10))
add_label_end(TS, logK,    "log K", colors$ink,   pos = 1, off = 0.25)
add_label_end(TS, ev[[3]], "pMOM",  colors$mom,   pos = 3, off = 0.15)
add_label_end(TS, ev[[4]], "piMOM", colors$imom,  pos = 1, off = 0.25)
add_label_end(TS, ev[[1]], "local", colors$local, pos = 1, off = 0.25)
add_panel_label("c")

dev.off()


# ==============================================================================
# 4. PRINT NUMBERS FOR TEXT
# ==============================================================================

cat("\nRates:\n")
print(format(tab, digits = 3), row.names = FALSE)
cat(sprintf("\ntau: iMOM %.4f | pMOM %.4f | normal slab %.1f\n",
            TAU_IMOM, TAU_MOM, TAU_NORM))
bf <- sapply(prior_arms, function(a) exp(tail(get(a$key, "H0")$mean, 1)))
cat(sprintf("BF against one spurious break at T = %d:\n", max(TS)))
for (i in seq_along(prior_arms)) cat(sprintf("  %-20s %.3e\n", prior_arms[[i]]$lab, bf[i]))
cat(sprintf("\npiMOM vs Zellner at T = %d: %.3g x less evidence for a spurious break\n",
            max(TS), bf[1] / bf[4]))
