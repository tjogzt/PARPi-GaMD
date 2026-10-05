#!/usr/bin/env Rscript
# Fig_SI_Sampling_Toolkit: sampling-resolution synthesis (Figure S17)
#
# Purpose:  Render the sampling-resolution synthesis figure (Fig_SI_Sampling_Toolkit): truncation scan + timescale scaling.
# Author:   Tao Zhu (tjogzt@gmail.com)
# Created:  2026-10-04 (header standardised 2026-10-05)
# Inputs:   data/analysis/p0_2_trunc_scan.csv ; data/analysis/s2_window_timescale.csv ; data/analysis/timescale_scaling.csv
# Run:      Rscript code/40_sampling_toolkit_figure.R   (from the repository root)
set.seed(49)
s2 <- read.csv("data/analysis/s2_window_timescale.csv", stringsAsFactors = FALSE)
ts <- read.csv("data/analysis/timescale_scaling.csv", stringsAsFactors = FALSE)
tr <- read.csv("data/analysis/p0_2_trunc_scan.csv", stringsAsFactors = FALSE)

CINNABAR <- "#C23531"
NIRA <- "#4DAF4A"; OLA <- "#E41A1C"; RUCA <- "#984EA3"
TALA <- "#FF7F00"; VELI <- "#377EB8"

pa <- s2[s2$window_ns == 2, ]
ord <- c("sys2_niraparib", "sys2_olaparib", "sys2_rucaparib", "sys2_nira_prot", "sys2_ruca_prot")
pa <- pa[match(ord, pa$tag), ]
lab <- c("niraparib", "olaparib", "rucaparib", "niraparib\n(prot.)", "rucaparib\n(prot.)")
col <- c(NIRA, OLA, RUCA, adjustcolor(NIRA, 0.45), adjustcolor(RUCA, 0.45))

draw_panels <- function() {
  par(mfrow = c(1, 3), family = "Arial", mar = c(4.2, 4.6, 2.5, 0.7),
      mgp = c(1.9, 0.5, 0), tcl = -0.25, cex.axis = 0.75, cex.lab = 0.8, las = 1)
  par(cex = 1)   # base R reduces base cex to 0.66 for mfrow layouts (>=3 panels); pin it

  # ---------- Panel A ----------
  plot(NA, xlim = c(0, 21.5), ylim = c(0.4, 6.0), xlab = "", ylab = "", axes = FALSE)
  rect(3.5, 0.4, 8.9, 6.0, col = adjustcolor("#888888", 0.16), border = NA)
  abline(v = 3.5, col = "#888888", lty = 3, lwd = 0.7); abline(v = 8.9, col = "#888888", lty = 3, lwd = 0.7)
  for (i in 1:5) {
    y <- 5.5 - i
    rect(0, y - 0.32, pa$c3_span_sd[i], y + 0.32, col = col[i], border = col[i], lwd = 0.7)
    text(pa$c3_span_sd[i] + 0.35, y, sprintf("%.1f", pa$c3_span_sd[i]), cex = 0.73, adj = 0)
  }
  axis(1, at = c(0, 5, 10, 15, 20))
  # axis title split to two lines: keeps the full term inside the narrow panel slot
  mtext("window-to-window SD", side = 1, line = 1.8, cex = 0.8)
  mtext("(kcal/mol)", side = 1, line = 2.7, cex = 0.8)
  xpd_bak <- par(xpd = TRUE)
  text(-0.5, 5.5 - (1:5), lab, adj = 1, cex = 0.75, xpd = NA)
  par(xpd = FALSE)
  mtext("signal range", side = 3, at = 6.2, line = -1.25, cex = 0.73, col = "#666666")
  mtext("A", side = 3, at = 0.1, line = 0.6, adj = 0, font = 2, cex = 0.92)
  mtext("noise vs signal", side = 3, at = 2.6, line = 0.6, adj = 0, cex = 0.73)

  # ---------- Panel B ----------
  plot(NA, xlim = c(1.2, 24.8), ylim = c(1.6, 4.5), xlab = "window length (ns)", ylab = "span SD (kcal/mol)", axes = FALSE)
  for (lg in c("talazoparib", "veliparib")) {
    seeds <- unique(ts$seed[ts$ligand == lg]); ltys <- c(1, 2, 3)
    for (k in seq_along(seeds)) {
      d <- ts[ts$ligand == lg & ts$seed == seeds[k], ]
      lines(d$window_ns, d$span_sd, col = ifelse(lg == "talazoparib", TALA, VELI), lty = ltys[k], lwd = 1.1)
    }
  }
  axis(1, at = c(2, 6, 12, 18, 24)); axis(2, at = c(2, 3, 4))
  legend("topright", legend = c("talazoparib (3 seeds)", "veliparib (3 seeds)"), col = c(TALA, VELI),
         lty = 1, lwd = 1.2, bty = "n", cex = 0.75, seg.len = 1.5, y.intersp = 1.05)
  mtext("B", side = 3, at = 1.3, line = 0.6, adj = 0, font = 2, cex = 0.92)
  mtext("noise vs sampling length", side = 3, at = 3.9, line = 0.6, adj = 0, cex = 0.73)

  # ---------- Panel C ----------
  plot(NA, xlim = c(5.5, 27.5), ylim = c(-1.05, 0.72), xlab = "production length (ns)", ylab = "Spearman rho", axes = FALSE)
  abline(h = 0, col = "#CCCCCC", lwd = 0.8)
  abline(h = -0.975, col = "#888888", lty = 2, lwd = 0.9)
  lines(tr$T_ns, tr$rho, col = CINNABAR, lwd = 1.2)
  points(tr$T_ns, tr$rho, pch = 16, col = CINNABAR, cex = 0.85)
  pcol <- ifelse(tr$p < 0.2, CINNABAR, "#666666")
  for (i in seq_len(nrow(tr))) {
    dy <- ifelse(tr$rho[i] >= 0, 0.15, ifelse(tr$rho[i] < -0.6, -0.19, -0.20))
    dx <- ifelse(tr$T_ns[i] > 25, -1.8, 0)  # keep the rightmost label inside the panel region
    text(tr$T_ns[i] + dx, tr$rho[i] + dy, sprintf("p=%.3f", tr$p[i]), cex = 0.73, col = pcol[i])
  }
  axis(1, at = c(8, 12, 16, 20, 24)); axis(2, at = c(-1.0, -0.5, 0, 0.5))
  text(6.1, -0.58, "p<0.05 needs", adj = 0, cex = 0.73, col = "#666666")
  text(6.1, -0.84, "|rho| >= 0.975", adj = 0, cex = 0.73, col = "#666666")
  mtext("C", side = 3, at = 5.7, line = 0.6, adj = 0, font = 2, cex = 0.92)
  mtext("association vs length", side = 3, at = 9.0, line = 0.6, adj = 0, cex = 0.73)
}

W <- 6.7; H <- 2.75
cairo_pdf("figures/pdf/Fig_SI_Sampling_Toolkit.pdf", width = W, height = H, pointsize = 11)
draw_panels()
dev.off()
png("figures/png/Fig_SI_Sampling_Toolkit.png", width = W, height = H, units = "in", res = 300, family = "Arial", pointsize = 11)
draw_panels()
dev.off()
cat("done. Panel A SDs:\n"); print(data.frame(system = pa$tag, sd = round(pa$c3_span_sd, 2)))
