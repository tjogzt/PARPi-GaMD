#!/usr/bin/env Rscript
# 07-master_figure.R — Combined master figure for PARPi design project
# Panels: (A) sys2 PMF, (B) sys1 PMF, (C) RMSD, (D) RMSF
# SCI figure: cairo_pdf, 7pt base_size

library(ggplot2)
library(data.table)
library(patchwork)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ===== Helper: read XVG PMF file =====
read_xvg <- function(f, label, domain = NULL) {
  lines <- readLines(f)
  data_start <- which(!grepl("^[@#]", lines) & lines != "")[1]
  x <- fread(f, skip = data_start - 1, header = FALSE)
  if (ncol(x) == 2) {
    setnames(x, c("cv", "pmf"))
    x[, pmf_err := 0]
  } else {
    setnames(x, c("cv", "pmf", "pmf_err"))
  }
  x[, system := label]
  if (!is.null(domain)) x[, domain := domain]
  x
}

# ===========================================
# Panel A: Sys2 PMF (APO vs talazoparib, CV1+CV2)
# ===========================================
sys2_cv1_apo  <- read_xvg("results/analysis/pmf-c3-sys2_APO_CV1_cv.dat.xvg", "APO", "CV1: Prot-DNA")
sys2_cv1_tlz  <- read_xvg("results/analysis/pmf-c3-sys2_talazoparib_CV1_cv.dat.xvg", "talazoparib", "CV1: Prot-DNA")
sys2_cv2_apo  <- read_xvg("results/analysis/pmf-c3-sys2_APO_CV2_cv.dat.xvg", "APO", "CV2: HD-ART")
sys2_cv2_tlz  <- read_xvg("results/analysis/pmf-c3-sys2_talazoparib_CV2_cv.dat.xvg", "talazoparib", "CV2: HD-ART")

sys2_pmf <- rbind(sys2_cv1_apo, sys2_cv1_tlz, sys2_cv2_apo, sys2_cv2_tlz)
sys2_pmf[, pmf_norm := pmf - min(pmf), by = .(domain, system)]
sys2_pmf[, system := factor(system, levels = c("APO", "talazoparib"))]

p_a <- ggplot(sys2_pmf, aes(x = cv, y = pmf_norm, color = system)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~domain, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c("APO" = "#2166AC", "talazoparib" = "#B2182B")) +
  labs(x = "CV (Å)", y = "Free Energy (kcal/mol)", color = NULL,
       title = "A  Sys2 PMF (CE-c3, full PARP1+DNA)") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_rect(fill = "grey95"),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 6.5),
    plot.margin = margin(2, 2, 2, 2)
  )

# ===========================================
# Panel B: Sys1 PMF (APO vs AZD5305)
# ===========================================
sys1_apo <- read_xvg("results/analysis/sys1_APO_pmf_c2.xvg", "APO")
sys1_azd <- read_xvg("results/analysis/sys1_AZD5305_pmf_c2.xvg", "AZD5305")
sys1_pmf <- rbind(sys1_apo, sys1_azd)
sys1_pmf[, pmf_norm := pmf - min(pmf), by = system]
sys1_pmf[, system := factor(system, levels = c("APO", "AZD5305"))]

p_b <- ggplot(sys1_pmf, aes(x = cv, y = pmf_norm, color = system)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = c("APO" = "#2166AC", "AZD5305" = "#B2182B")) +
  labs(x = "HD-ART Distance (Å)", y = "Free Energy (kcal/mol)", color = NULL,
       title = "B  Sys1 PMF (CE-c2, CAT domain)") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    plot.margin = margin(2, 2, 2, 2)
  )

# ===========================================
# Panel C: RMSD comparison (APO vs talazoparib)
# ===========================================
read_rmsd <- function(sys, dom) {
  f <- file.path("results/analysis", paste0(sys, "_rmsd_", dom, ".dat"))
  if (!file.exists(f)) return(NULL)
  x <- fread(f, skip = 1, col.names = c("frame", "rmsd"))
  x[, system := fifelse(grepl("APO", sys), "APO", "talazoparib")]
  x[, domain := dom]
  x
}

rmsd_list <- rbindlist(lapply(c("sys2_APO", "sys2_talazoparib"), function(s) {
  rbindlist(lapply(c("all", "HD", "ART", "DNA"), function(d) read_rmsd(s, d)))
}))
rmsd_list[, domain := factor(domain, levels = c("all", "HD", "ART", "DNA"),
                               labels = c("All Cα", "HD", "ART", "DNA"))]
rmsd_list[, time_ns := frame * 0.2]
rmsd_list[, system := factor(system, levels = c("APO", "talazoparib"))]

p_c <- ggplot(rmsd_list, aes(x = time_ns, y = rmsd, color = system)) +
  geom_line(linewidth = 0.3, alpha = 0.9) +
  facet_wrap(~domain, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("APO" = "#2166AC", "talazoparib" = "#B2182B")) +
  labs(x = "Time (ns)", y = "RMSD (Å)", color = NULL,
       title = "C  RMSD Time Series (Full PARP1+DNA)") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_rect(fill = "grey95"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 6.5),
    plot.margin = margin(2, 2, 2, 2)
  )

# ===========================================
# Panel D: RMSF comparison
# ===========================================
read_rmsf <- function(sys, dom) {
  f <- file.path("results/analysis", paste0(sys, "_rmsf_", dom, ".dat"))
  if (!file.exists(f)) return(NULL)
  x <- fread(f, skip = 1, col.names = c("resid", "rmsf"))
  x[, system := fifelse(grepl("APO", sys), "APO", "talazoparib")]
  x[, domain := dom]
  x
}

rmsf_list <- rbindlist(lapply(c("sys2_APO", "sys2_talazoparib"), function(s) {
  rbindlist(lapply(c("all", "HD", "ART", "DNA"), function(d) read_rmsf(s, d)))
}))
rmsf_list[, domain := factor(domain, levels = c("all", "HD", "ART", "DNA"),
                               labels = c("All Cα", "HD", "ART", "DNA"))]
rmsf_list[, system := factor(system, levels = c("APO", "talazoparib"))]

p_d <- ggplot(rmsf_list, aes(x = resid, y = rmsf, color = system)) +
  geom_line(linewidth = 0.3, alpha = 0.9) +
  facet_wrap(~domain, scales = "free", ncol = 4) +
  scale_color_manual(values = c("APO" = "#2166AC", "talazoparib" = "#B2182B")) +
  labs(x = "Residue", y = "RMSF (Å)", color = NULL,
       title = "D  RMSF by Residue (Full PARP1+DNA)") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_rect(fill = "grey95"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 6.5),
    plot.margin = margin(2, 2, 2, 2)
  )

# ===========================================
# Assemble: 4-panel 2x2 grid
# ===========================================
master <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

cairo_pdf("results/figures/07-master_figure.pdf",
          width = 190 / 25.4, height = 210 / 25.4)
print(master)
dev.off()

cat("Saved: results/figures/07-master_figure.pdf\n")
