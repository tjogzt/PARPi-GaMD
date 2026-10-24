#!/usr/bin/env Rscript
# Fig_SI_Sampling_Toolkit: sampling-resolution synthesis (Figure S17)
#
# Purpose:  Render the sampling-resolution synthesis figure (Fig_SI_Sampling_Toolkit,
#           SI Figure S17): truncation scan, window timescale, timescale scaling.
# Created:  2026-10-04 (header standardised 2026-10-05)
# Inputs:   data/analysis/p0_2_trunc_scan.csv ; data/analysis/s2_window_timescale.csv ;
#           data/analysis/timescale_scaling.csv
# Outputs:  figures/pdf/Fig_SI_Sampling_Toolkit.pdf
# Depends:  base R graphics only (cairo_pdf)
# Run:      Rscript code/40_sampling_toolkit_figure.R   (from the repository root)
set.seed(49)
s2 <- read.csv("data/analysis/s2_window_timescale.csv", stringsAsFactors = FALSE)
ts <- read.csv("data/analysis/timescale_scaling.csv", stringsAsFactors = FALSE)
tr <- read.csv("data/analysis/p0_2_trunc_scan.csv", stringsAsFactors = FALSE)

CINNABAR <- "#C23531"
NIRA <- "#4DAF4A"; OLA <- "#E41A1C"; RUCA <- "#984EA3"
TALA <- "#FF7F00"; VELI <- "#377EB8"

pa <- s2[s2$window_ns == 2, ]
# (2026-10-15) bars sorted by span SD, descending (largest at top), rather than
# the previous semantic grouping (three standard systems, then the two controls).
pa <- pa[order(-pa$c3_span_sd), ]
lab_map <- c(sys2_niraparib = "niraparib", sys2_olaparib = "olaparib",
             sys2_rucaparib = "rucaparib", sys2_nira_prot = "niraparib\n(prot.)",
             sys2_ruca_prot = "rucaparib\n(prot.)")
col_map <- c(sys2_niraparib = NIRA, sys2_olaparib = OLA, sys2_rucaparib = RUCA,
             sys2_nira_prot = adjustcolor(NIRA, 0.45), sys2_ruca_prot = adjustcolor(RUCA, 0.45))
lab <- lab_map[pa$tag]
col <- col_map[pa$tag]

draw_panels <- function() {
  par(mfrow = c(1, 3), family = "Arial", mar = c(3.6, 4.4, 2.2, 0.4),
      mgp = c(1.6, 0.5, 0), tcl = -0.25, cex.axis = 1.23, cex.lab = 1.25, las = 1)

  # ---------- Panel A ----------
  plot(NA, xlim = c(0, 21.5), ylim = c(0.52, 6.55), xlab = "", ylab = "", axes = FALSE)
  # (2026-10-15) signal-range band: extends from the raised baseline up to just
  # above the topmost bar (top bar top edge 5.97; band top 6.15).
  rect(2.6, 0.52, 5.3, 6.15, col = adjustcolor("#888888", 0.16), border = NA)
  segments(2.6, 0.52, 2.6, 6.15, col = "#888888", lty = 3, lwd = 0.7)
  segments(5.3, 0.52, 5.3, 6.15, col = "#888888", lty = 3, lwd = 0.7)
  for (i in 1:5) {
    y <- 6.85 - 1.2 * i
    rect(0, y - 0.32, pa$c3_span_sd[i], y + 0.32, col = col[i], border = col[i], lwd = 0.7)
    text(pa$c3_span_sd[i] + 0.35, y, sprintf("%.1f", pa$c3_span_sd[i]), cex = 1.23, adj = 0)
  }
  axis(1, at = c(0, 5, 10, 15, 20))
  # (2026-10-15) manual placement: the auto-centred label ran past the panel
  # boundary and its closing parenthesis was clipped.
  mtext("window-to-window SD (kcal/mol)", side = 1, line = 1.6, at = 9.5, adj = 0.5, cex = 0.84)
  xpd_bak <- par(xpd = TRUE)
  text(-0.5, 6.85 - 1.2 * (1:5), lab, adj = 1, cex = 1.23, xpd = NA)
  par(xpd = FALSE)
  mtext("signal range", side = 3, at = 3.95, line = -1.25, cex = 0.82, col = "#666666")
  mtext("A", side = 3, at = 0.1, line = 0.6, adj = 0, font = 2, cex = 0.92)
  mtext("noise vs signal", side = 3, at = 2.6, line = 0.6, adj = 0, cex = 0.82)

  # ---------- Panel B ----------
  plot(NA, xlim = c(1.2, 24.8), ylim = c(1.6, 4.5), xlab = "window length (ns)", ylab = "span SD (kcal/mol)", axes = FALSE)
  for (lg in c("talazoparib", "veliparib")) {
    seeds <- unique(ts$seed[ts$ligand == lg]); ltys <- c(1, 2, 3)
    for (k in seq_along(seeds)) {
      d <- ts[ts$ligand == lg & ts$seed == seeds[k], ]
      lines(d$window_ns, d$span_sd, col = ifelse(lg == "talazoparib", TALA, VELI), lty = ltys[k], lwd = 1.1)
    }
  }
  axis(1, at = c(2, 4, 6, 12, 24), labels = FALSE)
  text(c(2, 4, 6, 12, 24), par("usr")[3], labels = c("2", "4", "6", "12", "24"),
       pos = 1, offset = 0.45, cex = 1.15, xpd = NA); axis(2, at = c(2, 3, 4))
  legend("topright", legend = c("talazoparib (3 seeds)", "veliparib (3 seeds)"), col = c(TALA, VELI),
         lty = 1, lwd = 1.2, bty = "n", cex = 1.23, seg.len = 1.5, y.intersp = 1.05)
  mtext("B", side = 3, at = 1.3, line = 0.6, adj = 0, font = 2, cex = 0.92)
  mtext("noise vs sampling length", side = 3, at = 3.9, line = 0.6, adj = 0, cex = 0.82)

  # ---------- Panel C ----------
  # (2026-10-15) ylab moved ~0.85 margin lines further out: at mgp[1]=1.6 the
  # "Spearman rho" label touched the tick numbers (bbox overlap 1.8-5.5 pt).
  par(mgp = c(2.45, 0.5, 0))
  plot(NA, xlim = c(5.5, 27.5), ylim = c(-1.05, 0.72), xlab = "", ylab = "Spearman rho", axes = FALSE)
  abline(h = 0, col = "#CCCCCC", lwd = 0.8)
  abline(h = -0.975, col = "#888888", lty = 2, lwd = 0.9)
  lines(tr$T_ns, tr$rho, col = CINNABAR, lwd = 1.2)
  points(tr$T_ns, tr$rho, pch = 16, col = CINNABAR, cex = 1.0)
  pcol <- ifelse(tr$p < 0.2, CINNABAR, "#666666")
  for (i in seq_len(nrow(tr))) {
    dy <- ifelse(tr$rho[i] >= 0, 0.12, ifelse(tr$rho[i] < -0.6, -0.10, -0.155))
    text(tr$T_ns[i], tr$rho[i] + dy, sprintf("p=%.2f", tr$p[i]), cex = 1.23, col = pcol[i])
  }
  axis(1, at = c(8, 10, 14, 19.3, 23, 26), labels = FALSE)
  text(c(8, 10, 14, 19.3, 23, 26), par("usr")[3], labels = c("8", "10", "14", "19.3", "23", "26"),
       pos = 1, offset = 0.45, cex = 1.15, xpd = NA); axis(2, at = c(-1.0, -0.5, 0, 0.5))
  mtext("production length (ns)", side = 1, line = 1.6, cex = 0.84)
  text(8.2, -0.87, "p<0.05 needs |rho| >= 0.975", adj = 0, cex = 1.23, col = "#666666")
  mtext("C", side = 3, at = 5.7, line = 0.6, adj = 0, font = 2, cex = 0.92)
  mtext("association vs length", side = 3, at = 9.0, line = 0.6, adj = 0, cex = 0.82)
}

W <- 7.0; H <- 3.0
cairo_pdf("figures/pdf/Fig_SI_Sampling_Toolkit.pdf", width = W, height = H, pointsize = 11)
draw_panels()
dev.off()
png("figures/pdf/Fig_SI_Sampling_Toolkit.png", width = W, height = H, units = "in", res = 300, family = "Arial", pointsize = 11)
draw_panels()
dev.off()
cat("done. Panel A SDs:\n"); print(data.frame(system = pa$tag, sd = round(pa$c3_span_sd, 2)))
