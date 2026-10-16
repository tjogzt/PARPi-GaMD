#!/usr/bin/env Rscript
# 54_fig4_window_production.R — composite Fig 4: window variability (A) +
# production-length sensitivity (B), reusing code/33 panels via source().
#
# Purpose:  Composite Figure 4: descriptor–state–estimator variability (A, from code/33) +
#           production-length sensitivity (B), reusing the sourced panels.
# Created:  2026-10-10 (header standardised 2026-10-11)
# Inputs:   data/analysis/r6/r6c_window_timescale.csv ; data/analysis/timescale_scaling.csv (via sourced scripts)
# Outputs:  figures/pdf/Fig4_Window_Production.pdf
# Depends:  ggplot2, patchwork, data.table
# Run:      Rscript code/54_fig4_window_production.R   (from the repository root)
suppressMessages({library(ggplot2); library(patchwork); library(data.table)})

source("code/33_d4_diagnostics_figures.R")
wp <- p + labs(title = "A  Window noise vs between-system separations")
wq <- q + labs(title = "B  Span SD vs window length")

fig4 <- (wp | wq) + plot_layout(widths = c(1, 1))
cairo_pdf("figures/pdf/Fig4_Window_Production.pdf",
          width = 172/25.4, height = 72/25.4, pointsize = 8)
print(fig4)
invisible(dev.off())
ggsave("figures/pdf/Fig4_Window_Production.png", fig4,
       width = 172/25.4, height = 72/25.4, dpi = 220, type = "cairo", bg = "white")
cat("WROTE Fig4_Window_Production.pdf/png\n")
