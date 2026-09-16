#!/usr/bin/env Rscript
# Regenerate Fig_SI_PMF_Convergence.pdf (SI) from results/tables/pmf_convergence.csv.
#
# DATA PROVENANCE (see data_manifest.md): pmf_convergence.csv is LEGACY data — the original
# generating script is lost and the exact recipe is not fully recoverable from surviving
# artifacts (the dist_raw.npy samples do not reproduce the CSV's local-barrier metric).
# Its terminal well depths are independently validated by the pmf.npy / Table 4 chain.
# This script is the reproducible figure layer (CSV -> figure).
suppressMessages({library(data.table); library(ggplot2)})

theme_7pt <- theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid = element_blank(), legend.key.size = unit(0.3, "cm"),
        legend.position = "bottom")

cv <- fread("results/tables/pmf_convergence.csv")
stopifnot(all(c("inhibitor", "replicate", "ns", "well_depth") %in% names(cv)))
cv[, ns := as.numeric(ns)]
cv[, inhibitor := factor(inhibitor, levels = c("Talazoparib", "Veliparib", "AZD5305", "EB-47"))]

inib_colors <- c("Talazoparib" = "#C23531", "Veliparib" = "#3D6BA8",
                 "AZD5305" = "#9D2933", "EB-47" = "#177CB0")

# Mean +/- SEM per (inhibitor, ns), then LOESS smooth for trend + ribbon
sm <- cv[, .(mean_depth = mean(well_depth), se = sd(well_depth) / sqrt(.N)),
         by = .(inhibitor, ns)]

grid_ns <- seq(min(cv$ns), max(cv$ns), length.out = 200)
smooth_rib <- rbindlist(lapply(levels(cv$inhibitor), function(inh) {
  s <- sm[inhibitor == inh]
  if (nrow(s) < 4 || length(unique(s$ns)) < 4) return(NULL)  # EB-47: ns 0-2 only -> raw lines only
  fit_mean <- loess(mean_depth ~ ns, data = s, span = 0.6)
  fit_se   <- loess(se ~ ns, data = s, span = 0.6)
  data.table(inhibitor = inh, ns = grid_ns,
             mean_depth = predict(fit_mean, newdata = data.frame(ns = grid_ns)),
             se         = predict(fit_se,   newdata = data.frame(ns = grid_ns)))
}))

p <- ggplot(cv, aes(x = ns, y = well_depth)) +
  geom_line(aes(group = replicate), linewidth = 0.3, alpha = 0.4, color = "grey40") +
  geom_ribbon(data = smooth_rib, aes(x = ns, ymin = mean_depth - se, ymax = mean_depth + se,
                                     fill = inhibitor), alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = smooth_rib, aes(x = ns, y = mean_depth, color = inhibitor), linewidth = 0.8) +
  geom_vline(xintercept = 80, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  annotate("text", x = 80, y = -Inf, label = "80 ns", vjust = 1.5, hjust = -0.1,
           size = 2, color = "grey40") +
  scale_color_manual(values = inib_colors, name = NULL) +
  scale_fill_manual(values = inib_colors, guide = "none") +
  labs(x = "Cumulative production time (ns)", y = "S1 Well Depth (kcal/mol)") +
  theme_7pt

cairo_pdf("results/figures/Fig_SI_PMF_Convergence.pdf", width = 4.5, height = 3.4, pointsize = 7)
print(p)
dev.off()
cat("Saved: results/figures/Fig_SI_PMF_Convergence.pdf\n")
