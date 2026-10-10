#!/usr/bin/env Rscript
# r7_fig4_window_production.R — R7 composite Fig 4: window variability (A) +
# production-length sensitivity (B), reusing code/33 panels via source().
suppressMessages({library(ggplot2); library(patchwork); library(data.table)})

source("code/33-d4_diagnostics_figures.R")
wp <- p + labs(title = "A  Window noise vs between-system separations")
wq <- q + labs(title = "B  Span SD vs window length")

fig4 <- (wp | wq) + plot_layout(widths = c(1, 1))
cairo_pdf("results/figures/Fig4_Window_Production.pdf",
          width = 172/25.4, height = 72/25.4, pointsize = 8)
print(fig4)
invisible(dev.off())
ggsave("results/figures/Fig4_Window_Production.png", fig4,
       width = 172/25.4, height = 72/25.4, dpi = 220, type = "cairo", bg = "white")
cat("WROTE Fig4_Window_Production.pdf/png\n")
