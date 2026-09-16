#!/usr/bin/env Rscript
# 06-sys1_pmf_figure.R — Sys1 HD-ART PMF: APO vs AZD5305 (corrected weights)
# SCI figure: cairo_pdf, 7pt base_size

library(ggplot2)
library(data.table)
library(patchwork)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

read_xvg <- function(f, label) {
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
  x
}

# ===== Panel A: PMF (CE-c2, corrected weights) =====
apo_c2 <- read_xvg("results/analysis/sys1_APO_pmf_c2.xvg", "APO")
azd_c2 <- read_xvg("results/analysis/sys1_AZD5305_pmf_c2.xvg", "AZD5305")
pmf_data <- rbind(apo_c2, azd_c2)
pmf_data[, system := factor(system, levels = c("APO", "AZD5305"))]
pmf_data[, pmf_norm := pmf - min(pmf), by = system]

# Annotation: dV statistics
dv_anno <- data.table(
  system = factor(c("APO", "AZD5305"), levels = c("APO", "AZD5305")),
  label = c("dV = 3.40 ± 1.94 kcal/mol\naccel = 60,083×",
            "dV = 0.08 ± 0.43 kcal/mol\naccel = 1.49× (near-equilibrium)"),
  x = c(23.8, 22.5), y = c(16, 2)
)

p1 <- ggplot(pmf_data, aes(x = cv, y = pmf_norm, color = system, fill = system)) +
  geom_ribbon(aes(ymin = pmf_norm - pmf_err, ymax = pmf_norm + pmf_err),
              alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 0.5) +
  geom_text(data = dv_anno, aes(x = x, y = y, label = label, color = system),
            hjust = 0, size = 1.8, lineheight = 0.85, fontface = "italic") +
  scale_color_manual(values = c("APO" = "#2166AC", "AZD5305" = "#B2182B")) +
  scale_fill_manual(values = c("APO" = "#2166AC", "AZD5305" = "#B2182B")) +
  labs(x = "HD-ART COM Distance (Å)", y = "Free Energy (kcal/mol)",
       title = "A  PMF (CE-c2, corrected)") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 6.5),
    plot.margin = margin(3, 3, 3, 3)
  )

# ===== Panel B: CV distribution =====
cv_data <- rbind(
  fread("results/analysis/sys1_APO_cv.dat", col.names = "cv")[, system := "APO"],
  fread("results/analysis/sys1_AZD5305_cv.dat", col.names = "cv")[, system := "AZD5305"]
)
cv_data[, system := factor(system, levels = c("APO", "AZD5305"))]

p2 <- ggplot(cv_data, aes(x = cv, fill = system)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = c("APO" = "#2166AC", "AZD5305" = "#B2182B")) +
  labs(x = "HD-ART COM Distance (Å)", y = "Density",
       title = "B  CV Sampling (stride=100)") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = c(0.15, 0.85),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 6.5),
    plot.margin = margin(3, 3, 3, 3)
  )

# ===== Panel C: Time series =====
cv_ts <- copy(cv_data)
cv_ts[, frame := 1:.N, by = system]
cv_ts[, time_ns := frame * 0.2]  # stride=100, 2ps/frame → 200ps = 0.2ns per sampled frame

p3 <- ggplot(cv_ts, aes(x = time_ns, y = cv, color = system)) +
  geom_line(linewidth = 0.2, alpha = 0.85) +
  scale_color_manual(values = c("APO" = "#2166AC", "AZD5305" = "#B2182B")) +
  labs(x = "Time (ns)", y = "HD-ART Distance (Å)",
       title = "C  CV Trajectory") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 6.5),
    plot.margin = margin(3, 3, 3, 3)
  )

# ===== Panel D: Literature comparison =====
# "This work" rows: computed from the current trajectories (cv_summary above);
# "Prior (sys1)" rows: curated CSV with provenance (early trajectory version).
cv_summary <- cv_data[, .(cv_mean = mean(cv), cv_sd = sd(cv)), by = system]
lit_prior <- fread("data/01_curated/prior_sys1_cv.csv")
lit_this  <- cv_summary[, .(system, cv_mean, cv_sd)][, source := "This work"]
lit <- rbind(lit_prior[, .(system, cv_mean, cv_sd, source)], lit_this)
lit[, system := factor(system, levels = c("APO", "veliparib", "olaparib", "talazoparib", "AZD5305"))]

p4 <- ggplot(lit, aes(x = cv_mean, y = reorder(system, cv_mean, decreasing = TRUE),
                       color = source)) +
  geom_point(size = 2.5) +
  geom_linerange(aes(xmin = cv_mean - cv_sd, xmax = cv_mean + cv_sd),
                 linewidth = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("Prior (sys1)" = "grey50", "This work" = "#D94801")) +
  labs(x = "HD-ART Distance (Å)", y = NULL,
       title = "D  CV Across PARP1 Inhibitors") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = c(0.75, 0.20),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 6.5),
    plot.margin = margin(3, 3, 3, 3)
  )

# ===== Assemble =====
combined <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

cairo_pdf("results/figures/06-sys1_pmf_comparison.pdf",
          width = 190 / 25.4, height = 160 / 25.4)
print(combined)
dev.off()

cat("Saved: results/figures/06-sys1_pmf_comparison.pdf\n")
cat(sprintf("  APO: CV=%.2f±%.2f, dV=3.40±1.94\n",
            cv_summary[system=="APO", cv_mean], cv_summary[system=="APO", cv_sd]))
cat(sprintf("  AZD5305: CV=%.2f±%.2f, dV=0.08±0.43\n",
            cv_summary[system=="AZD5305", cv_mean], cv_summary[system=="AZD5305", cv_sd]))
cat(sprintf("  ΔCV: %.2f Å\n",
            cv_summary[system=="AZD5305", cv_mean] - cv_summary[system=="APO", cv_mean]))
