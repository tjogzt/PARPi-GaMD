#!/usr/bin/env Rscript
# 52_fig_c7_benchmark.R — Figure S13: synthetic benchmark of window-noise diagnostics
#
# Purpose:  SI Figure S13 — synthetic C7 benchmark of window-noise scaling and span-estimator
#           tail sensitivity (reads the archived benchmark CSVs).
# Created:  2026-10-10 (header standardised 2026-10-11)
# Inputs:   data/analysis/r6/c7_toy_benchmark.csv ; data/analysis/r6/c7_estimator_stats.csv
# Outputs:  figures/pdf/Fig_SI_ToyBenchmark.pdf
# Depends:  ggplot2, patchwork
# Run:      Rscript code/52_fig_c7_benchmark.R   (from the repository root)
suppressMessages({library(ggplot2); library(patchwork)})

d <- read.csv("data/analysis/r6/c7_toy_benchmark.csv")
st <- read.csv("data/analysis/r6/c7_estimator_stats.csv")

d$process <- factor(d$process,
  levels = c("fast_tau0.5ns", "medium_tau8ns", "slow_tau50ns", "naive_L^-1/2_ref"),
  labels = c("\u03c4 = 0.5 ns (fast)", "\u03c4 = 8 ns", "\u03c4 = 50 ns",
             "naive L\u207b\u00b9\u2044\u00b2 reference"))
cols <- c("0.5" = "#C23531", "8" = "#3D6BA8", "50" = "#177CB0", "ref" = "grey55")

pA <- ggplot(d, aes(L_ns, window_sd_mean)) +
  geom_line(data = subset(d, grepl("naive", process)),
            aes(group = process), linetype = "dashed", color = "grey55", linewidth = 0.6) +
  geom_line(data = subset(d, !grepl("naive", process)),
            aes(color = process), linewidth = 0.8) +
  geom_point(data = subset(d, !grepl("naive", process)),
             aes(color = process), size = 1.6) +
  scale_x_log10(breaks = c(2, 4, 6, 12, 24, 48, 96)) +
  scale_color_manual(values = setNames(c("#C23531", "#3D6BA8", "#177CB0"),
                                       c("\u03c4 = 0.5 ns (fast)", "\u03c4 = 8 ns", "\u03c4 = 50 ns")), name = NULL) +
  labs(x = "Window length L (ns, log scale)", y = "Segment-range SD (\u03c3 units)",
       title = "A  Window noise vs length:\nknown-truth processes") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.97, 1.00), legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 9),
        axis.text = element_text(size = 9), axis.title = element_text(size = 9),
        plot.title = element_text(size = 10),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5))

b <- data.frame(
  est = rep(c("min\u2013max", "5\u201395% quantile"), each = 2),
  kind = rep(c("baseline", "spiked"), 2),
  mean = c(st$mean_base[st$estimator == "minmax_span"], st$mean_spiked[st$estimator == "minmax_span"],
           st$mean_base[st$estimator == "q90_span"], st$mean_spiked[st$estimator == "q90_span"]))
b$est <- factor(b$est, levels = c("min\u2013max", "5\u201395% quantile"))

pB <- ggplot(b, aes(est, mean, fill = kind)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", mean)), position = position_dodge(width = 0.7),
            vjust = -0.6, size = 3.2, family = "Arial") +
  annotate("text", x = 1, y = 6.15, label = "+11.4%", size = 3.2, family = "Arial", fontface = "bold") +
  annotate("text", x = 2, y = 4.15, label = "+2.3%", size = 3.2, family = "Arial", fontface = "bold") +
  scale_fill_manual(values = c("baseline" = "grey55", "spiked" = "#C23531"), name = NULL) +
  coord_cartesian(ylim = c(0, 6.6)) +
  labs(x = NULL, y = "Mean span (\u03c3 units, 200 replicates)",
       title = "B  Tail sensitivity: rare short spikes") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(), legend.position = c(0.80, 0.90),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        axis.text = element_text(size = 9), axis.title = element_text(size = 9),
        plot.title = element_text(size = 10),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5))

## ---- Panel C: positive control (planted-association recovery, tie-aware) ----
c7g <- read.csv("data/analysis/r6/c7g_positive_control.csv")
c7gx <- read.csv("data/analysis/r6/c7g_positive_control_extra.csv")
sig_levels <- c(2, 4.6, 8, 15, 18.4)
c7g$sigma_lab <- factor(sprintf("\u03c3 = %g kcal/mol", c7g$sigma_kcal),
                        levels = sprintf("\u03c3 = %g kcal/mol", sig_levels))
c7gx$sigma_lab <- factor(sprintf("\u03c3 = %g kcal/mol", c7gx$sigma_kcal),
                         levels = sprintf("\u03c3 = %g kcal/mol", sig_levels))
sig_cols <- setNames(c("#C23531", "#9D2933", "#3D6BA8", "#177CB0", "#4DAF4A"),
                     sprintf("\u03c3 = %g kcal/mol", sig_levels))
sig_ltys <- setNames(c("dashed", "solid", "solid", "solid", "solid"),
                     sprintf("\u03c3 = %g kcal/mol", sig_levels))
pC <- ggplot(c7g, aes(delta_kcal, detect_rate, color = sigma_lab, linetype = sigma_lab)) +
  annotate("rect", xmin = 2.6, xmax = 5.3, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.55) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "grey40", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.0) +
  geom_point(data = c7gx, shape = 1, size = 1.7, stroke = 0.8, show.legend = FALSE) +
  scale_color_manual(values = sig_cols, name = NULL) +
  scale_linetype_manual(values = sig_ltys, name = NULL) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 50)) +
  annotate("label", x = 9.5, y = 0.86, label = "this study's separations\n(2.6\u20135.3 kcal/mol)",
           size = 3.0, family = "Arial", lineheight = 0.9, hjust = 0, label.size = 0.12,
           label.padding = unit(0.08, "lines"), fill = "white", alpha = 0.85) +
  labs(x = "planted separation between adjacent ranks (kcal/mol)",
       y = "detection rate (p < 0.05)",
       title = "C  Positive control: planted-association recovery (tie-aware; 1,500 worlds/cell)",
       subtitle = "open circles: extended grid (2.6/3.3/5.3 kcal/mol; 3,000 worlds/cell)") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.direction = "horizontal",
        legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = 7),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "transparent"),
        axis.text = element_text(size = 9), axis.title = element_text(size = 9),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5))

p <- (pA | pB) / pC
cairo_pdf("figures/pdf/Fig_SI_ToyBenchmark.pdf", width = 170/25.4, height = 114.5/25.4, pointsize = 8)
print(p)
invisible(dev.off())
cat("Saved: figures/pdf/Fig_SI_ToyBenchmark.pdf\n")
