# 29_d4_consistency_matrix.R -- Estimator x metric consistency matrix across the four
# chemical-state panels (main-text Fig. 6). R6 rebuild: reads the canonical long-format
# grid data/analysis/r6/r6c_estimatrix.csv (written by scripts/s38_estimator_matrix.py);
# the previous wide data/analysis/d4_consistency_matrix.csv (niraparib-corrected only,
# pre-audit values incl. empty-bin "eff" column) is superseded.
#
# Purpose:  Render the descriptor x state x estimator consistency matrix (standalone
#           diagnostic export; embedded by 53_fig3_state_matrix.R).
# Created:  2026-10-01 (header standardised 2026-10-05; R6 rebuild 2026-10-10)
# Inputs:   data/analysis/r6/r6c_estimatrix.csv ; data/analysis/r6/r6_panel_values.csv
# Outputs:  results/figures/Fig_D4_ConsistencyMatrix.pdf (+png; regenerable byproducts)
# Depends:  ggplot2, data.table, common/helpers.R
# Run:      Rscript code/29_d4_consistency_matrix.R   (from the repository root)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
suppressMessages({library(ggplot2); library(data.table)})

d <- fread("data/analysis/r6/r6c_estimatrix.csv")
d[, rho_num := suppressWarnings(as.numeric(rho))]
d[, p_num   := suppressWarnings(as.numeric(p))]
d[, est := factor(estimator, levels = c("mm", "q90", "t95", "eff"))]
d[, metric := factor(metric, levels = c("AAI", "S1", "S2CV2", "S2CV1"))]
d[, panel := factor(panel, levels = c("neutral", "nira_prot", "ruca_prot", "combined"),
                    labels = c("Neutral", "Niraparib-prot.",
                               "Rucaparib-prot.", "Audit-corrected"))]
d[, lab := ifelse(is.na(rho_num), "n.d.",
                  sprintf("%+.2f\n(p %.3f)", ifelse(rho_num >= 0, rho_num + 1e-9, rho_num - 1e-9), p_num))]

p <- ggplot(d, aes(est, metric, fill = rho_num)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = lab), size = 2.8, family = "Arial", lineheight = 0.95) +
  facet_wrap(~panel, nrow = 2) +
  scale_fill_gradient2(low = "#C23531", mid = "#F7F3EE", high = "#177CB0",
                       midpoint = 0, limits = c(-1, 1), na.value = "grey92",
                       name = "Spearman \u03c1") +
  scale_x_discrete(labels = c(mm = "min\u2013max", q90 = "5\u201395%",
                              t95 = "95% trunc", eff = "eff. support")) +
  scale_y_discrete(labels = c(S2CV1 = "S2 CV1", S2CV2 = "S2 CV2",
                              S1 = "S1", AAI = "span ratio")) +
  labs(x = "Span estimator", y = "",
       title = "Estimator \u00d7 metric consistency across chemical-state panels",
       subtitle = "n = 5, tie-aware exact permutation; red = inverse, blue = positive; no cell attains nominal significance (all exact p \u2265 0.267)") +
  theme_minimal(base_family = "Arial") +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 8),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        strip.text = element_text(size = 8, face = "bold"))
ggsave("results/figures/Fig_D4_ConsistencyMatrix.pdf", p,
       width = 180/25.4, height = 100/25.4, device = cairo_pdf, pointsize = 8)
ggsave("results/figures/Fig_D4_ConsistencyMatrix.png", p,
       width = 180/25.4, height = 100/25.4, dpi = 300, type = "cairo", bg = "white")
cat("WROTE Fig_D4_ConsistencyMatrix.pdf/png (4-panel grid)\n")
