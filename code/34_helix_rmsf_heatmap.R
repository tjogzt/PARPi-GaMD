#!/usr/bin/env Rscript
# 34_helix_rmsf_heatmap.R — Per-helix ΔRMSF (vs APO) heatmap (SI Figure S16)
# =============================================================================
# Purpose: diverging heatmap of the per-helix mean RMSF relative to the APO
#          control (ΔRMSF = system − APO) for SI section
#          "Per-helix RMSF decomposition of the HD".
# Inputs:  values transcribed from SI Table S9 (rmsf_recomp pipeline; seven S2
#          systems x five HD helices alphaB--alphaF). The Δ matrix is computed
#          in-script from the APO row; the transcription is self-checked
#          against the table's value range (1.8--7.4).
# Outputs: figures/pdf/Fig_SI_Helix_RMSF_Heatmap.pdf
# Style:   diverging scale (stone-blue = compressed vs APO; cinnabar red = expanded),
#          symmetric ±2.5 Å; signed cell values annotated; Arial >= 8 pt.
#          (Δ version adopted by the user 2026-10-05, replacing the earlier
#          absolute-means rendering of the same data.)
# =============================================================================
#
# Purpose:  Per-helix dRMSF (vs APO) heatmap for the HD subdomain (SI Figure S16).
# Created:  2026-10-03 (values updated 2026-10-06; header standardised 2026-10-05)
# Inputs:   data/01_curated/helix_rmsf_table_s9.csv (SI Table S9 values)
# Outputs:  figures/pdf/Fig_SI_Helix_RMSF_Heatmap.pdf
# Depends:  base R graphics only (cairo_pdf)
# Run:      Rscript code/34_helix_rmsf_heatmap.R   (from the repository root)
set.seed(49)

sys_all <- c("APO", "AZD5305", "Niraparib", "Olaparib", "Rucaparib", "Talazoparib", "Veliparib")
hel <- c("\u03B1B", "\u03B1C", "\u03B1D", "\u03B1E", "\u03B1F")
# Values: data/01_curated/helix_rmsf_table_s9.csv (SI Table S9; rerun_202609 recomputation).
s9 <- read.csv("data/01_curated/helix_rmsf_table_s9.csv", stringsAsFactors = FALSE)
s9 <- s9[match(sys_all, s9$system), ]
stopifnot(!any(is.na(s9$system)))
m <- as.matrix(s9[, c("alphaB", "alphaC", "alphaD", "alphaE", "alphaF")])
rownames(m) <- s9$system
stopifnot(identical(dim(m), c(7L, 5L)), all(m >= 1.8 & m <= 7.4))  # Table S9 range
d <- sweep(m[2:7, , drop = FALSE], 2, m[1, ], "-")                 # 6 x 5 ΔRMSF
sys <- sys_all[2:7]
stopifnot(max(abs(d)) <= 2.4)                                      # sane magnitudes

zlim <- c(-2.5, 2.5)
dr <- colorRampPalette(c("#24466E", "#3D6BA8", "#A9BFD6", "#F8F5ED",
                         "#E9B49E", "#C23531", "#8A1F1A"))(100)
colf <- function(v) dr[pmax(1, pmin(100, round((v - zlim[1]) / diff(zlim) * 99) + 1))]

draw_delta <- function() {
  par(family = "Arial", ps = 12, cex = 9/12, col.axis = "grey20", xpd = NA,
      oma = c(0.4, 0.4, 0.4, 0.4))
  layout(matrix(1:2, nrow = 1), widths = c(1, 0.19))
  # ---- main heatmap (6 rows) ----
  par(mar = c(2.0, 6.8, 1.5, 0.2))
  plot(0, 0, type = "n", xlim = c(0, 5), ylim = c(0, 6), axes = FALSE,
       xlab = "", ylab = "")
  for (i in 1:6) for (j in 1:5) {
    v <- d[i, j]
    rect(j - 1, 6 - i, j, 7 - i, col = colf(v), border = "white", lwd = 0.8)
    lab <- if (abs(v) < 0.005) "0.00" else sprintf("%+.2f", v)
    text(j - 0.5, 6.5 - i, lab, cex = 9/12,
         col = ifelse(abs(v) > 1.65, "white", "grey15"))
  }
  for (j in 1:5) mtext(hel[j], 1, at = j - 0.5, line = 0.35, cex = 9/12)
  for (i in 1:6) mtext(sys[i], 2, at = 6.5 - i, line = 0.35, las = 1, adj = 1, cex = 9/12)
  mtext("A", 3, at = 0, adj = 0, line = 0.15, font = 2, cex = 10/12)
  # ---- diverging colour bar (aligned top/bottom with the main panel) ----
  par(mar = c(2.0, 0.1, 1.5, 1.8))
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = zlim, axes = FALSE,
       xlab = "", ylab = "")
  yy <- seq(zlim[1], zlim[2], length.out = 100)
  for (k in 1:99) rect(0.0, yy[k], 0.34, yy[k + 1], col = colf((yy[k] + yy[k + 1]) / 2), border = NA)
  rect(0.0, zlim[1], 0.34, zlim[2], border = "grey30", lwd = 0.5)
  for (v in -2:2) {
    lw <- if (v == 0) 1.4 else 0.9
    segments(0.34, v, if (v == 0) 0.52 else 0.475, v, lwd = lw, col = "grey20")
    text(0.52, v, v, adj = c(0, 0.5), cex = 9/12, col = "grey20")
  }
  mtext("\u0394RMSF vs APO (\u00C5)", 4, line = 0.55, cex = 9/12)
}

cairo_pdf("figures/pdf/Fig_SI_Helix_RMSF_Heatmap.pdf",
          width = 118/25.4, height = 60/25.4, pointsize = 8)
draw_delta(); invisible(dev.off())
cat("Saved: figures/pdf/Fig_SI_Helix_RMSF_Heatmap.pdf (ΔRMSF vs APO)\n")
