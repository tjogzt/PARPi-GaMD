#!/usr/bin/env Rscript
# 53_fig3_state_matrix.R — composite Fig 3: state panels (A/B, from code/21) +
# descriptor–state–estimator matrix (C, from code/29). Zero-drift reuse via source().
#
# Purpose:  Composite Figure 3: state panels (A/B, from code/21) + descriptor–state–estimator
#           matrix (C, from code/29). Zero-drift reuse via source().
# Created:  2026-10-10 (header standardised 2026-10-11)
# Inputs:   data/analysis/S1_S2_ratios.csv ; data/analysis/r6/r6c_estimatrix.csv (via sourced scripts)
# Outputs:  figures/pdf/Fig3_State_Matrix.pdf
# Depends:  ggplot2, ggrepel, patchwork, data.table
# Run:      Rscript code/53_fig3_state_matrix.R   (from the repository root)
suppressMessages({library(ggplot2); library(ggrepel); library(patchwork); library(data.table)})

source("code/29_d4_consistency_matrix.R")
pm <- p
pm <- pm + labs(title = "C  Descriptor\u2013state\u2013estimator matrix",
                subtitle = "no cell attains nominal significance (all exact p \u2265 0.267); red = inverse, blue = positive")
suppressMessages(library(ggplot2))

source("code/21_spearman_correlation.R")
pa <- p_a + theme(plot.title = element_text(margin = ggplot2::margin(r = 14)))
pb <- p_b + theme(plot.title = element_text(margin = ggplot2::margin(l = 0)))  # 2026-10-13: was l = 14; aligns the B title with its y-axis
pm <- pm + theme(plot.title = element_text(margin = ggplot2::margin(t = 4, r = 0, b = 5.5, l = 0)))  # 2026-10-13: +4pt above the C title (widens the gap to row 1; b = 5.5 preserves title-subtitle spacing)

fig3 <- (pa | pb) / pm + plot_layout(heights = c(1, 1.65))
fig3 <- fig3 & theme(plot.margin = ggplot2::margin(2, 3, 2, 3, "pt"))
cairo_pdf("figures/pdf/Fig3_State_Matrix.pdf",
          width = 172/25.4, height = 163/25.4, pointsize = 8)
print(fig3)
invisible(dev.off())
ggsave("figures/pdf/Fig3_State_Matrix.png", fig3,
       width = 172/25.4, height = 163/25.4, dpi = 220, type = "cairo", bg = "white")
cat("WROTE Fig3_State_Matrix.pdf/png\n")
