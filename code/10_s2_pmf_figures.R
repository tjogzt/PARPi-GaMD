#!/usr/bin/env Rscript
# 10-s2_pmf_figures.R — S2 PMF figures (Protein–DNA + HD–ART, 7 systems)
#
# Purpose:  Render the S2 PMF figures (protein–DNA and HD–ART, seven systems: overlay + facet panels).
# Author:   Tao Zhu (tjogzt@gmail.com)
# Created:  2026-09-15 (header standardised 2026-10-05)
# Depends:  dplyr, ggplot2, patchwork
# Run:      Rscript code/10_s2_pmf_figures.R   (from the repository root)
library(ggplot2)
library(dplyr)
library(patchwork)

# Shared helpers: read_pmf / theme_7pt / extract_features (single source).
args_h <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", args_h))) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
  source(file.path(script_dir, "..", "common", "helpers.R"))
} else {
  source("common/helpers.R")
}


# ---- Config ----------------------------------------------------------------
data_dir  <- "data/analysis"
out_dir   <- "figures/pdf"
stats_dir <- "results/analysis"   # CSV statistics tables (regenerable)
dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# System metadata
sys_meta <- data.frame(
  ligand   = c("APO", "AZD5305", "olaparib", "veliparib", "niraparib", "rucaparib", "talazoparib"),
  label    = c("APO", "AZD5305", "Olaparib", "Veliparib", "Niraparib", "Rucaparib", "Talazoparib"),
  class    = c("APO", "Unknown", "Type_II", "Type_III", "Type_III", "Type_III", "Type_II"),
  color    = c("grey40", "darkorange", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
  lty      = c("dotted", "dashed", "solid", "solid", "solid", "solid", "solid"),
  stringsAsFactors = FALSE
)

# ---- Load PMF data ---------------------------------------------------------


# Load CV1 (Protein–DNA) and CV2 (HD–ART) for all systems
load_s2_pmf <- function(cv) {
  all_pmf <- list()
  for (i in seq_len(nrow(sys_meta))) {
    lig <- sys_meta$ligand[i]
    f <- file.path(data_dir, sprintf("pmf-c3-sys2_%s_%s_cv.dat.xvg", lig, cv))
    if (!file.exists(f)) {
      message("Missing: ", f)
      next
    }
    df <- read_pmf(f)
    df$ligand  <- lig
    df$label   <- sys_meta$label[i]
    df$class   <- sys_meta$class[i]
    df$color   <- sys_meta$color[i]
    all_pmf[[lig]] <- df
  }
  bind_rows(all_pmf)
}

pmf_cv1 <- load_s2_pmf("CV1")
pmf_cv2 <- load_s2_pmf("CV2")

# Factor ordering
ligand_order <- c("APO", "AZD5305", "olaparib", "talazoparib", "niraparib", "rucaparib", "veliparib")
pmf_cv1$ligand <- factor(pmf_cv1$ligand, levels = ligand_order)
pmf_cv2$ligand <- factor(pmf_cv2$ligand, levels = ligand_order)

# ---- Theme -----------------------------------------------------------------
theme_pmf <- theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    plot.margin       = ggplot2::margin(1, 1, 1, 3, unit = "mm"),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(size = 9, face = "bold"),
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 8),
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 8),
    legend.key.size   = unit(0.3, "cm")
  )

# ---- Fig A: CV1 Protein–DNA overlay ----------------------------------------
p_cv1 <- ggplot(pmf_cv1, aes(x = RC, y = PMF_norm, color = label, linetype = label)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label)) +
  scale_linetype_manual(values = setNames(sys_meta$lty, sys_meta$label)) +
  labs(x = "Protein–DNA COM Distance (Å)", y = "PMF (kcal/mol)",
       color = "System", linetype = "System") +
  theme_pmf + theme_open

# ---- Fig B: CV2 HD–ART overlay ---------------------------------------------
p_cv2 <- ggplot(pmf_cv2, aes(x = RC, y = PMF_norm, color = label, linetype = label)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label)) +
  scale_linetype_manual(values = setNames(sys_meta$lty, sys_meta$label)) +
  labs(x = "HD–ART COM Distance (Å)", y = "PMF (kcal/mol)",
       color = "System", linetype = "System") +
  theme_pmf + theme_open

# ---- Fig C: CV1 multi-panel (one per system) -------------------------------
p_cv1_list <- lapply(ligand_order, function(lig) {
  df <- filter(pmf_cv1, ligand == lig)
  meta <- sys_meta[sys_meta$ligand == lig, ]
  pmf_min <- df$RC[which.min(df$PMF)]
  
  ggplot(df, aes(x = RC, y = PMF_norm)) +
    geom_line(color = meta$color, linewidth = 0.4) +
    geom_vline(xintercept = pmf_min, color = meta$color, 
               linetype = "dashed", linewidth = 0.3) +
    labs(title = meta$label, x = "Protein–DNA (Å)", y = "PMF (kcal/mol)") +
    annotate("text", x = pmf_min, y = max(df$PMF_norm) * 0.85,
             label = sprintf("%.1f", pmf_min), 
             hjust = 1.15, size = 2.9, color = meta$color) +
    theme_pmf + theme(legend.position = "none")
})

# ---- Fig D: CV2 multi-panel ------------------------------------------------
p_cv2_list <- lapply(ligand_order, function(lig) {
  df <- filter(pmf_cv2, ligand == lig)
  meta <- sys_meta[sys_meta$ligand == lig, ]
  pmf_min <- df$RC[which.min(df$PMF)]
  
  ggplot(df, aes(x = RC, y = PMF_norm)) +
    geom_line(color = meta$color, linewidth = 0.4) +
    geom_vline(xintercept = pmf_min, color = meta$color, 
               linetype = "dashed", linewidth = 0.3) +
    labs(title = meta$label, x = "HD–ART (Å)", y = "PMF (kcal/mol)") +
    annotate("text", x = pmf_min, y = max(df$PMF_norm) * 0.85,
             label = sprintf("%.1f", pmf_min), 
             hjust = 1.15, size = 2.9, color = meta$color) +
    theme_pmf + theme(legend.position = "none")
})

# ---- Save figures ----------------------------------------------------------
cairo_pdf(file.path(out_dir, "Fig_S2_pmf_prot_dna_overlay.pdf"), width = 149/25.4, height = 99.3/25.4, pointsize = 8)
print(p_cv1)
dev.off()

cairo_pdf(file.path(out_dir, "Fig_S2_pmf_hd_art_overlay.pdf"), width = 149/25.4, height = 99.3/25.4, pointsize = 8)
print(p_cv2)
dev.off()

# Combined 2-panel
p_overlay <- wrap_plots(p_cv1, p_cv2, ncol = 1)
cairo_pdf(file.path(out_dir, "Fig_S2_pmf_overlay_combined.pdf"), width = 4.5, height = 6.5, pointsize = 8)
print(p_overlay)
dev.off()

# Multi-panel CV1
panel_cv1 <- wrap_plots(p_cv1_list, ncol = 4, nrow = 2)
cairo_pdf(file.path(out_dir, "Fig_S2_pmf_prot_dna_panels.pdf"), width = 157/25.4, height = 78.5/25.4, pointsize = 8)
print(panel_cv1)
dev.off()

# Multi-panel CV2
panel_cv2 <- wrap_plots(p_cv2_list, ncol = 4, nrow = 2)
cairo_pdf(file.path(out_dir, "Fig_S2_pmf_hd_art_panels.pdf"), width = 157/25.4, height = 78.5/25.4, pointsize = 8)
print(panel_cv2)
dev.off()

# ---- Summary stats ---------------------------------------------------------
cv1_stats <- pmf_cv1 %>%
  group_by(ligand, label, class) %>%
  summarise(
    n_points    = n(),
    rc_min      = RC[which.min(PMF)],
    pmf_range   = max(PMF_norm, na.rm = TRUE),
    rc_range    = max(RC) - min(RC),
    .groups     = "drop"
  )
cv2_stats <- pmf_cv2 %>%
  group_by(ligand, label, class) %>%
  summarise(
    n_points    = n(),
    rc_min      = RC[which.min(PMF)],
    pmf_range   = max(PMF_norm, na.rm = TRUE),
    rc_range    = max(RC) - min(RC),
    .groups     = "drop"
  )

cat("\n=== S2 CV1 (Protein–DNA) PMF Stats ===\n")
print(cv1_stats)
cat("\n=== S2 CV2 (HD–ART) PMF Stats ===\n")
print(cv2_stats)

write.csv(cv1_stats, file.path(stats_dir, "S2_CV1_stats.csv"), row.names = FALSE)
write.csv(cv2_stats, file.path(stats_dir, "S2_CV2_stats.csv"), row.names = FALSE)

message("===== S2 Figures saved to ", out_dir, " =====")
