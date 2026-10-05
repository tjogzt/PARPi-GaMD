#!/usr/bin/env Rscript
# 34-helix_rmsf_heatmap.R — Per-helix ΔRMSF (vs APO) heatmap (SI Figure S16)
# =============================================================================
# Purpose: diverging heatmap of the per-helix mean RMSF relative to the APO
#          control (ΔRMSF = system − APO) for SI section
#          "Per-helix RMSF decomposition of the HD".
# Inputs:  values transcribed from SI Table S9 (rmsf_recomp pipeline; seven S2
#          systems x five HD helices alphaB--alphaF). The Δ matrix is computed
#          in-script from the APO row; the transcription is self-checked
#          against the table's value range (1.8--7.4).
# Outputs: results/figures/Fig_SI_Helix_RMSF_Heatmap.pdf
# Style:   diverging scale (石青 blue = compressed vs APO; 朱砂 red = expanded),
#          symmetric ±2.5 Å; signed cell values annotated; Arial >= 8 pt.
#          (Δ version adopted by the user 2026-10-05, replacing the earlier
#          absolute-means rendering of the same data.)
# =============================================================================
set.seed(49)

sys_all <- c("APO", "AZD5305", "Niraparib", "Olaparib", "Rucaparib", "Talazoparib", "Veliparib")
hel <- c("\u03B1B", "\u03B1C", "\u03B1D", "\u03B1E", "\u03B1F")
m <- rbind(
  c(6.21, 2.25, 3.86, 4.44, 2.91),   # APO (baseline row; not plotted)
  c(7.36, 3.30, 4.67, 5.14, 2.64),
  c(3.83, 2.25, 3.25, 2.78, 2.22),
  c(4.77, 2.74, 4.29, 3.34, 2.62),
  c(4.18, 2.06, 3.01, 3.00, 2.08),
  c(4.33, 2.34, 3.47, 3.59, 3.13),
  c(5.19, 1.86, 3.38, 4.03, 2.54))
stopifnot(identical(dim(m), c(7L, 5L)), all(m >= 1.8 & m <= 7.4))  # Table S9 range
d <- sweep(m[2:7, , drop = FALSE], 2, m[1, ], "-")                 # 6 x 5 ΔRMSF
sys <- sys_all[2:7]
stopifnot(max(abs(d)) <= 2.4)                                      # sane magnitudes

zlim <- c(-2.5, 2.5)
dr <- colorRampPalette(c("#24466E", "#3D6BA8", "#A9BFD6", "#F8F5ED",
                         "#E9B49E", "#C23531", "#8A1F1A"))(100)
colf <- function(v) dr[pmax(1, pmin(100, round((v - zlim[1]) / diff(zlim) * 99) + 1))]

draw_delta <- function() {
  par(family = "Arial", ps = 8, cex = 8/12, col.axis = "grey20", xpd = NA,
      oma = c(0.4, 0.4, 0.4, 0.4))
  layout(matrix(1:2, nrow = 1), widths = c(1, 0.19))
  # ---- main heatmap (6 rows) ----
  par(mar = c(2.0, 6.8, 0.6, 0.2))
  plot(0, 0, type = "n", xlim = c(0, 5), ylim = c(0, 6), axes = FALSE,
       xlab = "", ylab = "")
  for (i in 1:6) for (j in 1:5) {
    v <- d[i, j]
    rect(j - 1, 6 - i, j, 7 - i, col = colf(v), border = "white", lwd = 0.8)
    lab <- if (abs(v) < 0.005) "0.00" else sprintf("%+.2f", v)
    text(j - 0.5, 6.5 - i, lab, cex = 8/12,
         col = ifelse(abs(v) > 1.65, "white", "grey15"))
  }
  for (j in 1:5) mtext(hel[j], 1, at = j - 0.5, line = 0.35, cex = 8.5/12)
  for (i in 1:6) mtext(sys[i], 2, at = 6.5 - i, line = 0.35, las = 1, adj = 1, cex = 8.5/12)
  # ---- diverging colour bar (aligned top/bottom with the main panel) ----
  par(mar = c(2.0, 0.1, 0.6, 1.8))
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = zlim, axes = FALSE,
       xlab = "", ylab = "")
  yy <- seq(zlim[1], zlim[2], length.out = 100)
  for (k in 1:99) rect(0.0, yy[k], 0.34, yy[k + 1], col = colf((yy[k] + yy[k + 1]) / 2), border = NA)
  rect(0.0, zlim[1], 0.34, zlim[2], border = "grey30", lwd = 0.5)
  for (v in -2:2) {
    lw <- if (v == 0) 1.4 else 0.9
    segments(0.34, v, if (v == 0) 0.52 else 0.475, v, lwd = lw, col = "grey20")
    text(0.52, v, v, adj = c(0, 0.5), cex = 8/12, col = "grey20")
  }
  mtext("\u0394RMSF vs APO (\u00C5)", 4, line = 0.55, cex = 8/12)
}

cairo_pdf("results/figures/Fig_SI_Helix_RMSF_Heatmap.pdf",
          width = 118/25.4, height = 66/25.4, pointsize = 8)
draw_delta(); invisible(dev.off())
cat("Saved: results/figures/Fig_SI_Helix_RMSF_Heatmap.pdf (ΔRMSF vs APO)\n")
