# 29-d4_consistency_matrix.R -- Estimator x metric consistency matrix (Fig 9).
# Data source: data/analysis/d4_consistency_matrix.csv, written by
# scripts/d4_consistency_matrix.py (mean-curve rows) plus the *_tab rows
# written by scripts/audit_r1_table_convention.py. All values are read from
# CSV; nothing is hard-coded here except display labels/factors.
#
# Purpose:  Render the estimator x metric consistency matrix (Figure 9) from the derived CSV.
# Author:   Tao Zhu (tjogzt@gmail.com)
# Created:  2026-10-01 (header standardised 2026-10-05)
# Inputs:   data/analysis/d4_consistency_matrix.csv
# Depends:  data.table, ggplot2
# Run:      Rscript code/29_d4_consistency_matrix.R   (from the repository root)
suppressMessages({library(ggplot2); library(data.table)})

d <- fread("data/analysis/d4_consistency_matrix.csv")

vals <- rbind(
  data.frame(metric = d$metric, est = "mm",      rho = d$mm_rho),
  data.frame(metric = d$metric, est = "q90",     rho = d$q90_rho),
  data.frame(metric = d$metric, est = "trunc95", rho = d$trunc95_rho),
  data.frame(metric = d$metric, est = "eff",     rho = d$eff_rho)
)
vals$est <- factor(vals$est, levels=c("mm","q90","trunc95","eff"))
vals$metric <- factor(vals$metric,
  levels=rev(c("S2_CV1_C3","S2_CV2_C3","S1_C3","S1_C3_tab","AAI","AAI_tab")))
p <- ggplot(vals, aes(est, metric, fill=rho)) +
  geom_tile(color="white", size=0.8) +
  geom_text(aes(label=sprintf("%+.2f", rho)), size=4.5, family="Arial") +
  scale_fill_gradient2(low="#C23531", mid="#F7F3EE", high="#177CB0",
                       midpoint=0, limits=c(-1,1), name="Spearman \u03c1") +
  scale_x_discrete(labels=c(mm="min\u2013max", q90="5\u201395%", trunc95="95% trunc",
                            eff="eff. support")) +
  scale_y_discrete(labels=c(S2_CV1_C3="S2 CV1 C3", S2_CV2_C3="S2 CV2 C3",
                            S1_C3="S1 C3 (curve)", S1_C3_tab="S1 C3 (table)",
                            AAI="AAI (curve)", AAI_tab="AAI (table)")) +
  labs(x="Span estimator", y="",
       title="Estimator \u00d7 metric consistency of the trapping association",
       subtitle="n = 5, tie-aware exact permutation, niraparib-corrected panel; red = inverse, blue = positive") +
  theme_minimal(base_family="Arial") +
  theme(axis.text=element_text(size=8, color="black"),
        axis.title=element_text(size=8),
        plot.title=element_text(size=9, face="bold"),
        plot.subtitle=element_text(size=8),
        legend.position="right")
ggsave("figures/pdf/Fig_D4_ConsistencyMatrix.pdf", p, width=140/25.4, height=85.9/25.4, device=cairo_pdf)
ggsave("figures/png/Fig_D4_ConsistencyMatrix.png", p, width=140/25.4, height=85.9/25.4, dpi=300, type="cairo")
cat("WROTE Fig_D4_ConsistencyMatrix.pdf/png\n")
