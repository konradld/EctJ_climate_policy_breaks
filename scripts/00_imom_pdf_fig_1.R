################################################################################
# iMom prior PDF
#
# Figure 1: PDF of iMom prior with different scale paramters tau
#
# Output: Plot of iMom PDF
################################################################################


library(mombf)

# Generated data
tau_5pct <- priorp2g(0.05, 1, prior = "iMom")
tau_1pct <- priorp2g(0.01, 1, prior = "iMom")

x <- seq(-7, 7, length.out = 4001)

y5 <- dimom(x, tau = tau_5pct)
y1 <- dimom(x, tau = tau_1pct)

# COLORS
orange <- "#cf5e00"
blue <- "#0074b7"
orange_fill <- grDevices::adjustcolor(orange, alpha.f = 0.28)
blue_fill <- grDevices::adjustcolor(blue, alpha.f = 0.28)


pdf("./output/figure_1.pdf", width = 5.833333, height = 2.083333,
    family = "serif", pointsize = 11, useDingbats = FALSE)

par(mar = c(1.5, 2, 0.5, 2), mgp = c(1.45, 0.26, 0),
    tcl = -0.25, xaxs = "i", yaxs = "i", family = "serif")
plot(NA, xlim = c(-7.3, 7.3), ylim = c(0, 0.20), axes = FALSE,
     xlab = "", ylab = "")

# The overlaid shaded probability masses have bounds -1 and 1
polygon(c(-1, x[x > -1 & x < 1], 1),
        c(0, y5[x > -1 & x < 1], 0), col = orange_fill, border = NA)
polygon(c(-1, x[x > -1 & x < 1], 1),
        c(0, y1[x > -1 & x < 1], 0), col = blue_fill, border = NA)
lines(x, y5, col = orange, lwd = 2.0)
lines(x, y1, col = blue, lwd = 2.0, lty = 2)

# Axes
axis(1, at = c(-6, -4, -2, -1, 0, 1, 2, 4, 6),
     labels = c("-6", "-4", "-2", "-1", "0", "1", "2", "4", "6"),
     cex.axis = 0.87, col.axis = "black", col.ticks = "black", 
     lwd = 0, lwd.ticks = 0.5,
     tcl = 0.25)
axis(2, at = c(0.05, 0.10, 0.15, 0.20), pos = 0,
     labels = c("0.05", "0.10", "0.15", "0.20"), las = 1,
     cex.axis = 0.87, col.axis = "black", col.ticks = "black", 
     lwd = 0, lwd.ticks = 0.5, 
     tcl = 0.25)
axis(2, at = c(0, 
               0.01, 0.02, 0.03, 0.04, 
               0.06, 0.07, 0.08, 0.09, 
               0.11, 0.12, 0.13, 0.14, 
               0.16, 0.17, 0.18, 0.19), pos = 0,
     labels = rep(NA, 17), las = 1,
     cex.axis = 0.01, col.axis = "black", col.ticks = "black", 
     lwd = 0, lwd.ticks = 0.5,
     tcl = 0.15)
segments(-7.3, 0, 7.3, 0, col = "black", lwd = 1.35)
segments(0, 0, 0, 0.2, col = "black", lwd = 1.35)

# Legend
segments(-5.78, 0.169, -5.20, 0.169, col = orange, lwd = 2.0)
text(-5.16, 0.169, expression(iMOM~(tau == 1.92)), adj = c(0, 0.5), cex = 0.83)
segments(-5.78, 0.144, -5.20, 0.144, col = blue, lwd = 2.0, lty = 2)
text(-5.16, 0.144, expression(iMOM~(tau == 3.32)), adj = c(0, 0.5), cex = 0.83)

# Inset callouts and arrows
rect(-2.46, 0.042, -1.55, 0.077, col = "#f1d6c0", border = "black", lwd = 0.9)
rect(-2.46, 0.002, -1.55, 0.039, col = "#b9d8e8", border = "black", lwd = 0.9)
text(-2.00, 0.0595, "5%", cex = 0.85, font = 3)
text(-2.00, 0.0205, "1%", cex = 0.85, font = 3)
arrows(-1.54, 0.0595, -1.03, 0.0555, length = 0.09, angle = 25, lwd = 0.95)
arrows(-1.54, 0.0205, -1.03, 0.0175, length = 0.09, angle = 25, lwd = 0.95)

text(7.5, -0.0055, expression(gamma/sigma), adj = c(0, 0),
     cex = 1.0, xpd = NA, col = )

dev.off()

