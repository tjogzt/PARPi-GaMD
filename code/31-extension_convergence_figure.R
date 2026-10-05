#!/usr/bin/env Rscript
# 31-extension_convergence_figure.R — Extension-panel results figure (SI Figure S15)
# =============================================================================
# Three-panel figure for SI section 2.8 (Extension Panel):
#   (A) S1 PMF overlay (C3, 200 ns) for fluzoparib / pamiparib / senaparib
#       Source: results/analysis/sys1_<ligand>_pmf_c3.xvg (min-normalized per curve)
#   (B) S1 (200 ns) and S2 CV2 (22 ns) well depths (C3) ± SD across C1--C3
#       S1 values: SI §2.8 text (54.0±8.0, 48.7±5.8, 52.2±7.8 kcal/mol)
#       S2 CV2: results/analysis/s2_dbe_final.csv (C1..C3 -> mean + SD)
#       Grey band: seven-system S1 well-depth range from
#       results/analysis/s1_dbe_unified_wells.csv (48.0--62.9 kcal/mol)
#   (C) Cumulative S1 well depth (C3) vs production time + 22 ns truncation
#       Source: results/analysis/extension_s1_convergence.csv
# Panel letters sit at the top-left corner of each panel canvas.
# Supersedes the former single-panel convergence-only version (2026-10-05 upgrade,
# adopted from the four-team review; draft approved by user with letter placement).
# =============================================================================
set.seed(49)

PROJ <- "."   # run from repository root
OUT  <- "results/figures/Fig_SI_Extension_Convergence.pdf"

read_pmf <- function(f) {
  L <- readLines(f); L <- L[grepl("^\\s*[0-9]", L)]
  m <- do.call(rbind, lapply(strsplit(trimws(L), "\\s+"), function(x) as.numeric(x[1:2])))
  m <- m[order(m[, 1]), ]
  if (nrow(m) < 5) stop("read_pmf: too few valid rows in ", f)
  data.frame(x = m[, 1], y = m[, 2])
}

ligs <- c("fluzoparib", "pamiparib", "senaparib")
cols <- c(fluzoparib = "#C23531", pamiparib = "#3D6BA8", senaparib = "#9D2933")

pmfs <- lapply(ligs, function(l)
  read_pmf(sprintf("%s/results/analysis/sys1_%s_pmf_c3.xvg", PROJ, l)))
names(pmfs) <- ligs

s1u  <- read.csv(sprintf("%s/results/analysis/s1_dbe_unified_wells.csv", PROJ))
band <- range(s1u$C3)

s2 <- read.csv(sprintf("%s/results/analysis/s2_dbe_final.csv", PROJ))
s2cv2 <- s2[s2$metric == "CV2" & s2$ligand %in% ligs, ]
s2cv2 <- s2cv2[match(ligs, s2cv2$ligand), ]
s2sd  <- apply(s2cv2[, c("C1", "C2", "C3")], 1, sd)

s1_c3 <- c(54.0, 48.7, 52.2)   # SI §2.8 text (C3)
s1_sd <- c(8.0, 5.8, 7.8)      # SI §2.8 text (± SD across C1--C3)

conv <- read.csv(sprintf("%s/results/analysis/extension_s1_convergence.csv", PROJ))

draw_all <- function() {
  par(cex = 8/12, ps = 8, family = "Arial", col.axis = "grey20", col.lab = "grey10",
      tcl = -0.25, mgp = c(1.6, 0.5, 0), xaxs = "i", yaxs = "i")
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(1.12, 1))

  ## ---- Panel A: S1 PMF overlay ----
  par(mar = c(2.7, 3.1, 1.8, 0.7))
  plot(0, 0, type = "n", xlim = c(22.4, 26.6), ylim = c(0, 58), axes = FALSE,
       xlab = "", ylab = "")
  axis(1, at = 23:26); axis(2, at = seq(0, 50, 10))
  for (l in ligs) { p <- pmfs[[l]]; p$y <- p$y - min(p$y)
    lines(p$x, p$y, col = cols[l], lwd = 1.1) }
  mtext("C3 distance (\u00C5)", 1, line = 1.75, cex = 8/12)
  mtext("PMF (kcal/mol)", 2, line = 2.0, cex = 8/12)
  usr <- par("usr")
  mtext("A", 3, at = usr[1], adj = 0, line = 0.55, font = 2, cex = 10/12)
  legend("top", legend = ligs, col = cols[ligs], lwd = 1.8, bty = "n",
         cex = 8/12, seg.len = 1.2, x.intersp = 0.5, y.intersp = 0.8, horiz = TRUE)
  mtext("S1 PMF (C3, 200 ns)", 3, line = 0.55, adj = 0.5, cex = 8.5/12, col = "grey30")

  ## ---- Panel B: wells bars + seven-system reference band ----
  par(mar = c(3.0, 3.1, 1.8, 0.7))
  gx <- c(1.05, 2.55, 4.05)
  plot(0, 0, type = "n", xlim = c(0.1, 5.0), ylim = c(30, 82), axes = FALSE,
       xlab = "", ylab = "")
  rect(0.1, band[1], 5.0, band[2], col = adjustcolor("grey60", 0.10), border = NA)
  abline(h = band, col = "grey55", lty = 3, lwd = 0.7)
  for (i in 1:3) {
    b1 <- gx[i] - 0.38; b2 <- gx[i] + 0.38; w <- 0.30
    rect(b1 - w, 30, b1 + w, s1_c3[i], col = cols[ligs[i]], border = "grey20", lwd = 0.5)
    segments(b1, s1_c3[i] - s1_sd[i], b1, s1_c3[i] + s1_sd[i], lwd = 0.9)
    segments(b1 - 0.09, s1_c3[i] - s1_sd[i], b1 + 0.09, s1_c3[i] - s1_sd[i], lwd = 0.9)
    segments(b1 - 0.09, s1_c3[i] + s1_sd[i], b1 + 0.09, s1_c3[i] + s1_sd[i], lwd = 0.9)
    rect(b2 - w, 30, b2 + w, s2cv2$C3[i], col = adjustcolor(cols[ligs[i]], 0.42),
         border = "grey20", lwd = 0.5)
    segments(b2, s2cv2$C3[i] - s2sd[i], b2, s2cv2$C3[i] + s2sd[i], lwd = 0.9)
    segments(b2 - 0.09, s2cv2$C3[i] - s2sd[i], b2 + 0.09, s2cv2$C3[i] - s2sd[i], lwd = 0.9)
    segments(b2 - 0.09, s2cv2$C3[i] + s2sd[i], b2 + 0.09, s2cv2$C3[i] + s2sd[i], lwd = 0.9)
    if (i < 3) abline(v = (gx[i] + gx[i + 1]) / 2, col = "grey88", lwd = 0.6)
  }
  axis(2, at = seq(30, 80, 10))
  axis(1, at = gx, labels = FALSE, tick = FALSE)
  mtext(c("Fluzoparib", "Pamiparib", "Senaparib"), 1, at = gx, line = 1.35, cex = 7.6/12)
  mtext(c("S1", "S2", "S1", "S2", "S1", "S2"), 1,
        at = as.vector(rbind(gx - 0.38, gx + 0.38)), line = 0.25, cex = 6.8/12, col = "grey35")
  mtext("Well depth (kcal/mol)", 2, line = 2.0, cex = 8/12)
  usr <- par("usr")
  mtext("B", 3, at = usr[1], adj = 0, line = 0.55, font = 2, cex = 10/12)
  legend("topright", legend = c("S1 (200 ns)", "S2 CV2 (22 ns)"),
         fill = c("grey30", adjustcolor("grey30", 0.42)), border = "grey20",
         bty = "n", cex = 7.8/12, seg.len = 0.8, y.intersp = 0.85)
  mtext("wells \u00B1 C1\u2013C3 spread; grey band = 7-system S1 range", 3,
        line = 0.55, adj = 0.5, cex = 7.2/12, col = "grey30")

  ## ---- Panel C: convergence + 22 ns truncation ----
  par(mar = c(2.7, 3.1, 1.8, 0.7))
  plot(0, 0, type = "n", xlim = c(0, 205), ylim = c(35, 95), axes = FALSE,
       xlab = "", ylab = "")
  axis(1); axis(2, at = seq(40, 90, 10))
  for (l in ligs) { cd <- conv[conv$ligand == l, ]; cd <- cd[order(cd$time_ns), ]
    lines(cd$time_ns, cd$well_depth, col = cols[l], lwd = 1.0)
    segments(cd$time_ns[nrow(cd)] - 22, tail(cd$well_depth, 1),
             cd$time_ns[nrow(cd)] + 3, tail(cd$well_depth, 1),
             col = cols[l], lty = 3, lwd = 0.7) }
  abline(v = 22, col = "grey35", lty = 2, lwd = 0.8)
  text(30, 92, "22 ns truncation", adj = 0, cex = 7.2/12, col = "grey35")
  mtext("Production time (ns)", 1, line = 1.75, cex = 8/12)
  mtext("Cumulative S1 well depth (C3)", 2, line = 2.0, cex = 8/12)
  usr <- par("usr")
  mtext("C", 3, at = usr[1], adj = 0, line = 0.55, font = 2, cex = 10/12)
  mtext("S1 convergence & sensitivity", 3, line = 0.55, adj = 0.5, cex = 8.5/12, col = "grey30")
}

cairo_pdf(OUT, width = 160/25.4, height = 115/25.4,
          pointsize = 8, family = "Arial")
draw_all(); invisible(dev.off())

# self-check: all three panels must carry real data
stopifnot(all(vapply(pmfs, function(p) diff(range(p$y)) > 30, TRUE)),
          nrow(s2cv2) == 3, diff(band) > 10, nrow(conv) >= 30)
cat("Saved:", OUT, "\n")
