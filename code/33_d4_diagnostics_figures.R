# 33-d4_diagnostics_figures.R -- D4 diagnostic figures (Fig 8 / Fig 10).
# Windowed-bootstrap noise panel (Fig D4a) and consistency-timescale panel.
# Reads data/analysis/convergence_diag.csv; the between-system signal (the
# range of the three re-simulated S2 CV1 C3 spans, 62.0 / 58.5 / 59.3
# kcal/mol; cf. Table S5) is derived from data/analysis/s2_dbe_final.csv,
# with a guard assert against the original documented anchor (3.5).
#
# Purpose:  Render the D4 diagnostics figures: windowed-bootstrap noise panel + consistency-timescale panel.
# Created:  2026-10-01 (header standardised 2026-10-05)
# Inputs:   data/analysis/convergence_diag.csv ; data/analysis/timescale_scaling.csv
# Outputs:  figures/pdf/Fig_D4_Timescale.pdf ; figures/pdf/Fig_D4_WindowNoise.pdf
# Depends:  data.table, ggplot2
# Run:      Rscript code/33_d4_diagnostics_figures.R   (from the repository root)
suppressMessages({library(ggplot2); library(data.table)})

# --- Fig D4a: windowed-bootstrap noise vs between-system signal ----------
d <- fread("data/analysis/convergence_diag.csv")
d <- d[cv == "CV1"]  # the primary association metric
d[, system := sub("^s2_", "", system)]
d[, system := factor(system, levels = unique(d[order(span_sd)]$system))]
# Between-system signal: derived from the consolidated table (the three
# re-simulated systems' S2 CV1 C3 spans; niraparib 62.0 = protonated state,
# olaparib 58.5, rucaparib 59.3 kcal/mol; cf. SI Table S5).
# Original hard-coded anchor (3.5) retained as the guard below.
dbe <- fread("data/analysis/s2_dbe_final.csv")
sig3 <- dbe[metric == "CV1" & ligand %in% c("niraparib", "olaparib", "rucaparib"), C3]
signal <- max(sig3) - min(sig3)
stopifnot(abs(signal - 3.5) < 0.05)

d[, xpos := as.numeric(system) - 0.5]
p <- ggplot(d, aes(x = xpos, y = span_sd)) +
  geom_hline(yintercept = signal, linetype = "dashed", color = "#9D2933", linewidth = 0.7) +
  annotate("text", x = -0.015, y = 4.2, label = sprintf("between-system signal %.1f", signal),
           color = "#9D2933", size = 3.0, family = "Arial", angle = 90, hjust = 0, vjust = 0.5) +
  geom_col(fill = "#3D6BA8", alpha = 0.85, width = 0.62) +
  scale_y_continuous(limits = c(0, 18), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(breaks = seq_along(levels(d$system)) - 0.5,
                     labels = levels(d$system),
                     limits = c(-0.12, 3.08),
                     expand = expansion(0)) +
  labs(x = NULL, y = "Windowed C3-span SD (kcal/mol)") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 0.9)
ggsave("figures/pdf/Fig_D4_WindowNoise.pdf", p, device = cairo_pdf, width = 140/25.4, height = 101.8/25.4, pointsize = 8)
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
  labs(x = "Window length (ns) — log scale", y = "Reweighted span SD (kcal/mol)",
       color = NULL) +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.position = c(0.15, 0.88), legend.text = element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        aspect.ratio = 0.9)
ggsave("figures/pdf/Fig_D4_Timescale.pdf", q, device = cairo_pdf, width = 140/25.4, height = 101.8/25.4, pointsize = 8)
cat("Saved: Fig_D4_Timescale.pdf\n")
