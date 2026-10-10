#!/usr/bin/env Rscript
# r6f_fig_c7_benchmark.R — Figure S18: synthetic benchmark of window-noise diagnostics
suppressMessages({library(ggplot2); library(patchwork)})

d <- read.csv("results/analysis/r6/c7_toy_benchmark.csv")
st <- read.csv("results/analysis/r6/c7_estimator_stats.csv")

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
  labs(x = "Window length L (ns, log scale)", y = "Window-to-window SD of segment range (\u03c3 units)",
       title = "A  Window noise vs length: known-truth processes") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.83, 0.22), legend.background = element_rect(fill = "white", linewidth = 0.2),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 7.6),
        axis.text = element_text(size = 8), axis.title = element_text(size = 8))

b <- data.frame(
  est = rep(c("min\u2013max", "5\u201395% quantile"), each = 2),
  kind = rep(c("baseline", "spiked"), 2),
  mean = c(st$mean_base[st$estimator == "minmax_span"], st$mean_spiked[st$estimator == "minmax_span"],
           st$mean_base[st$estimator == "q90_span"], st$mean_spiked[st$estimator == "q90_span"]))
b$est <- factor(b$est, levels = c("min\u2013max", "5\u201395% quantile"))

pB <- ggplot(b, aes(est, mean, fill = kind)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", mean)), position = position_dodge(width = 0.7),
            vjust = -0.6, size = 2.8, family = "Arial") +
  annotate("text", x = 1, y = 5.95, label = "+11.4%", size = 2.9, family = "Arial", fontface = "bold") +
  annotate("text", x = 2, y = 3.95, label = "+2.3%", size = 2.9, family = "Arial", fontface = "bold") +
  scale_fill_manual(values = c("baseline" = "grey55", "spiked" = "#C23531"), name = NULL) +
  coord_cartesian(ylim = c(0, 6.6)) +
  labs(x = NULL, y = "Mean span (\u03c3 units, 200 replicates)",
       title = "B  Tail sensitivity: rare short spikes") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(), legend.position = c(0.80, 0.90),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 7.6),
        legend.background = element_rect(fill = "white", linewidth = 0.2),
        axis.text = element_text(size = 8), axis.title = element_text(size = 8))

p <- pA | pB
cairo_pdf("results/figures/Fig_SI_ToyBenchmark.pdf", width = 170/25.4, height = 72/25.4, pointsize = 8)
print(p)
invisible(dev.off())
cat("Saved: results/figures/Fig_SI_ToyBenchmark.pdf\n")
