#!/usr/bin/env Rscript
# 30-extension_convergence_figure.R — extension S1 convergence curves (SI figure)
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

p <- ggplot(d, aes(x = time_ns, y = well_depth, color = ligand)) +
  geom_hline(yintercept = c(48.7, 52.2, 54.0), linetype = "dotted", alpha = 0.35) +
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
