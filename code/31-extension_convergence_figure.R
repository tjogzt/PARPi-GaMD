#!/usr/bin/env Rscript
# 31-extension_convergence_figure.R — extension S1 convergence curves (SI figure)
#
# Purpose:   Plot the cumulative well depth vs production time for the three
#            extension inhibitors (fluzoparib, pamiparib, senaparib).
# Inputs:    results/analysis/extension_s1_convergence.csv
# Outputs:   results/figures/Fig_SI_Extension_Convergence.pdf
# Depends:   ggplot2
library(ggplot2)

d <- read.csv("results/analysis/extension_s1_convergence.csv")
d$ligand <- factor(d$ligand, levels = c("fluzoparib", "pamiparib", "senaparib"),
                   labels = c("Fluzoparib", "Pamiparib", "Senaparib"))
cols <- c("Fluzoparib" = "#C23531", "Pamiparib" = "#3D6BA8", "Senaparib" = "#9D2933")

# Final 200-ns well depths per ligand (from the cumulative CSV's last row), used as
# dotted reference lines — no hardcoded values.
final_vals <- aggregate(well_depth ~ ligand, data = d, FUN = function(x) tail(x, 1))
names(final_vals)[2] <- "final"
p <- ggplot(d, aes(x = time_ns, y = well_depth, color = ligand)) +
  geom_hline(data = final_vals, aes(yintercept = final), linetype = "dotted", alpha = 0.35) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_color_manual(values = cols, name = "") +
  scale_x_continuous(breaks = seq(0, 200, 40)) +
  labs(x = "Cumulative production time (ns)",
       y = "S1 HD--ART well depth (kcal/mol, C3)") +
  theme_bw(base_size = 9) +
  theme(legend.position = c(0.82, 0.82),
        legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Arial"))

ggsave("results/figures/Fig_SI_Extension_Convergence.pdf",
       p, width = 90, height = 62, units = "mm", device = cairo_pdf)
cat("Saved: Fig_SI_Extension_Convergence.pdf\n")
