#!/usr/bin/env Rscript
# r7_fig3_state_matrix.R — R7 composite Fig 3: state panels (A/B, from code/21) +
# descriptor–state–estimator matrix (C, from code/29). Zero-drift reuse via source().
suppressMessages({library(ggplot2); library(ggrepel); library(patchwork); library(data.table)})

source("code/29-d4_consistency_matrix.R")
pm <- p
pm <- pm + labs(title = "C  Descriptor\u2013state\u2013estimator matrix",
                subtitle = "no cell attains nominal significance (all exact p \u2265 0.27); red = inverse, blue = positive")
suppressMessages(library(ggplot2))

source("code/21-spearman_correlation.R")
pa <- p_a + theme(plot.title = element_text(margin = ggplot2::margin(r = 14)))
pb <- p_b + theme(plot.title = element_text(margin = ggplot2::margin(l = 14)))

fig3 <- (pa | pb) / pm + plot_layout(heights = c(1, 1.30))
cairo_pdf("results/figures/Fig3_State_Matrix.pdf",
          width = 172/25.4, height = 148/25.4, pointsize = 8)
print(fig3)
invisible(dev.off())
ggsave("results/figures/Fig3_State_Matrix.png", fig3,
       width = 172/25.4, height = 148/25.4, dpi = 220, type = "cairo", bg = "white")
cat("WROTE Fig3_State_Matrix.pdf/png\n")
