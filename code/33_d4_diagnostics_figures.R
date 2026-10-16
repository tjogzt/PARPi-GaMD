# 33_d4_diagnostics_figures.R -- D4 diagnostic figures (Fig 8 / Fig 10).
# Windowed-bootstrap noise panel (Fig D4a) and consistency-timescale panel.
# Reads data/analysis/r6/r6c_window_timescale.csv. The dashed line marks the
# upper edge of the audited panel's neighbouring between-system separations
# (2.6-5.3 kcal/mol across the four most closely spaced systems; R9 convention).
#
# Purpose:  Render the D4 diagnostics figures: windowed-bootstrap noise panel vs the
#           audited between-system signal, and the convergence-timescale scaling panel
#           (standalone diagnostic exports; embedded by 54_fig4_window_production.R).
# Created:  2026-10-01 (header standardised 2026-10-05; R9 convention 2026-10-11)
# Inputs:   data/analysis/r6/r6c_window_timescale.csv ; data/analysis/timescale_scaling.csv ;
#           data/analysis/r6/r6_panel_values.csv
# Outputs:  results/figures/Fig_D4_WindowNoise.pdf ; results/figures/Fig_D4_Timescale.pdf
#           (regenerable byproducts)
# Depends:  data.table, ggplot2
# Run:      Rscript code/33_d4_diagnostics_figures.R   (from the repository root)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
suppressMessages({library(ggplot2); library(data.table)})

# --- Fig D4a: windowed noise vs between-system signal (R6 values) ----------
d <- fread("data/analysis/r6/r6c_window_timescale.csv")
d <- d[window_ns == 2]  # the primary 2-ns window scale
d <- d[tag %in% c("sys2_niraparib", "sys2_olaparib", "sys2_rucaparib", "sys2_nira_prot", "sys2_ruca_prot")]  # five re-simulated systems
d[, system := sub("^sys2_", "", tag)]
d[, state := ifelse(grepl("_prot$", tag), "protonated (+1)", "neutral")]
d[system == "nira_prot", system := "niraparib+1"]
d[system == "ruca_prot", system := "rucaparib+1"]
d[, system := factor(system, levels = d[order(c3_span_sd)]$system)]
# Upper edge of the audited panel's neighbouring separations (2.6-5.3 kcal/mol,
# R9 convention): largest consecutive gap among the four most closely spaced
# S2 CV1 systems (talazoparib, veliparib-protonated, olaparib, niraparib-protonated).
pv_aud <- fread("data/analysis/r6/r6_panel_values.csv")
key4 <- c("talazoparib", "veliparib_prot", "olaparib", "niraparib_prot")
v4 <- sort(pv_aud$S2CV1_C3[match(key4, pv_aud$ligand)])
stopifnot(length(v4) == 4L, !any(is.na(v4)))
gaps <- diff(v4)
stopifnot(abs(min(gaps) - 2.6) < 0.05, abs(max(gaps) - 5.3) < 0.05)
signal <- max(gaps)

d[, xpos := as.numeric(system) - 0.5]
p <- ggplot(d, aes(x = xpos, y = c3_span_sd)) +
  geom_hline(yintercept = signal, linetype = "dashed", color = "#9D2933", linewidth = 0.7) +
  annotate("text", x = -0.015, y = 4.2, label = "separations 2.6-5.3",
           color = "#9D2933", size = 3.0, family = "Arial", angle = 90, hjust = 0, vjust = 0.5) +
  geom_col(aes(fill = state), alpha = 0.85, width = 0.62) +
  scale_fill_manual(values = c("neutral" = "#3D6BA8", "protonated (+1)" = "#C23531"), name = NULL) +
  theme(legend.position = c(0.80, 0.96), legend.text = element_text(size = 8),
        legend.key.size = unit(0.28, "cm"), legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent")) +
  scale_y_continuous(limits = c(0, 18), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(breaks = seq_along(levels(d$system)) - 0.5,
                     labels = levels(d$system),
                     limits = c(-0.12, 4.88),
                     expand = expansion(0)) +
  labs(x = NULL, y = "Windowed C3-span SD (kcal/mol)") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 0.9)
ggsave("results/figures/Fig_D4_WindowNoise.pdf", p, device = cairo_pdf, width = 140/25.4, height = 101.8/25.4, pointsize = 8)
cat("Saved: Fig_D4_WindowNoise.pdf\n")

# --- Fig D4b: convergence-timescale scaling -------------------------------
t <- fread("data/analysis/timescale_scaling.csv")
t[, drug := ifelse(grepl("tala", ligand), "talazoparib", "veliparib")]
q <- ggplot(t, aes(x = window_ns, y = span_sd, group = ligand,
                   color = drug, linetype = drug, shape = drug)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.7) +
  scale_color_manual(values = c("talazoparib" = "#C23531", "veliparib" = "#177CB0"), name = NULL) +
  scale_linetype_manual(values = c("talazoparib" = "solid", "veliparib" = "dashed"), name = NULL) +
  scale_shape_manual(values = c("talazoparib" = 16, "veliparib" = 17), name = NULL) +
  scale_x_log10(breaks = c(2, 4, 6, 12, 24)) +
  labs(x = "Window length (ns), log scale", y = "Reweighted span SD (kcal/mol)",
       color = NULL) +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.position = c(0.15, 0.88), legend.text = element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        aspect.ratio = 0.9)
ggsave("results/figures/Fig_D4_Timescale.pdf", q, device = cairo_pdf, width = 140/25.4, height = 101.8/25.4, pointsize = 8)
cat("Saved: Fig_D4_Timescale.pdf\n")
