#!/usr/bin/env Rscript
# 10_s2_pmf_figures.R — S2 PMF figures (Protein–DNA + HD–ART, 7 systems)
#
# Purpose:  Render the S2 PMF figures (protein–DNA and HD–ART, seven systems: overlay + facet panels).
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

# System metadata: single source of truth (labels / classes / colours / linetypes).
if (length(grep("^--file=", args_h))) {
  source(file.path(script_dir, "..", "common", "sys_meta.R"))
} else {
  source("common/sys_meta.R")
}


# ---- Config ----------------------------------------------------------------
data_dir  <- "data/analysis"
out_dir   <- "figures/pdf"
stats_dir <- "results/analysis"   # CSV statistics tables (regenerable)
dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# System metadata: sys_meta is defined in common/sys_meta.R (single source).

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
ligand_order <- sys_display_order   # panel display order (common/sys_meta.R)
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
# 2026-10-14 (Figure S3 standard): legend removed from the overlay; colour and
# linetype are identified by the shared legend in the panel grid instead.
p_cv1 <- ggplot(pmf_cv1, aes(x = RC, y = PMF_norm, color = label, linetype = label)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label),
                     breaks = sys_meta$label[sys_meta$ligand %in% ligand_order]) +
  scale_linetype_manual(values = setNames(sys_meta$lty, sys_meta$label),
                        breaks = sys_meta$label[sys_meta$ligand %in% ligand_order]) +
  labs(x = "Protein–DNA COM Distance (Å)", y = "PMF (kcal/mol)",
       color = NULL, linetype = NULL) +
  theme_pmf + theme_open + theme(legend.position = "none")

# ---- Fig B: CV2 HD–ART overlay ---------------------------------------------
# 2026-10-14 (Figure S4 standard): legend removed from the overlay; colour and
# linetype are identified by the shared legend in the panel grid instead.
p_cv2 <- ggplot(pmf_cv2, aes(x = RC, y = PMF_norm, color = label, linetype = label)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label),
                     breaks = sys_meta$label[sys_meta$ligand %in% ligand_order]) +
  scale_linetype_manual(values = setNames(sys_meta$lty, sys_meta$label),
                        breaks = sys_meta$label[sys_meta$ligand %in% ligand_order]) +
  labs(x = "HD–ART COM Distance (Å)", y = "PMF (kcal/mol)",
       color = NULL, linetype = NULL) +
  theme_pmf + theme_open + theme(legend.position = "none")

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
    theme_pmf + theme(legend.position = "none", aspect.ratio = 1)   # 2026-10-14: square
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
    theme_pmf + theme(legend.position = "none", aspect.ratio = 1)   # 2026-10-14: square
})

# ---- Save figures ----------------------------------------------------------
# CV1 overlay (Figure S3 block B, 2026-10-14 standard): legend-free, full width;
# stacked in the SI under the grid so rows 1+2 equal the overlay in height.
cairo_pdf(file.path(out_dir, "Fig_S2_pmf_prot_dna_overlay.pdf"), width = 178/25.4, height = 94/25.4, pointsize = 8)
print(p_cv1)
dev.off()

# CV2 overlay (Figure S4 block B, 2026-10-14 standard): legend-free, full width;
# stacked in the SI under the grid so rows 1+2 equal the overlay in height.
cairo_pdf(file.path(out_dir, "Fig_S2_pmf_hd_art_overlay.pdf"), width = 178/25.4, height = 94/25.4, pointsize = 8)
print(p_cv2)
dev.off()

# Combined 2-panel
p_overlay <- wrap_plots(p_cv1, p_cv2, ncol = 1)
cairo_pdf(file.path(out_dir, "Fig_S2_pmf_overlay_combined.pdf"), width = 4.5, height = 6.5, pointsize = 8)
print(p_overlay)
dev.off()

# Multi-panel CV1 -> 4+3 grid with the shared legend cell (Figure S3 blocks A, 2026-10-14)
names(p_cv1_list) <- ligand_order
p_cv1_leg <- p_cv1 + theme(legend.position = "bottom",
                           legend.text       = element_text(size = 8),
                           legend.key.size   = unit(0.30, "cm"),
                           legend.key.width  = unit(0.50, "cm"),
                           legend.spacing.y  = unit(0.02, "cm")) +
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
leg_grob <- {
  g <- ggplot2::ggplotGrob(p_cv1_leg)
  g$grobs[[which(vapply(g$grobs, function(x) x$name, "") == "guide-box")]]
}
row1 <- wrap_plots(p_cv1_list[c("APO", "AZD5305", "olaparib", "talazoparib")], ncol = 4)
row2 <- wrap_plots(c(p_cv1_list[c("niraparib", "rucaparib", "veliparib")], list(leg_grob)), ncol = 4)
panel_grid <- wrap_plots(row1, row2, ncol = 1)
cairo_pdf(file.path(out_dir, "Fig_S2_pmf_prot_dna_grid.pdf"), width = 178/25.4, height = 97/25.4, pointsize = 8)
print(panel_grid)
dev.off()

# Multi-panel CV2 -> 4+3 grid with the shared legend cell (Figure S4 blocks A, 2026-10-14)
names(p_cv2_list) <- ligand_order
p_cv2_leg <- p_cv2 + theme(legend.position = "bottom",
                           legend.text       = element_text(size = 8),
                           legend.key.size   = unit(0.30, "cm"),
                           legend.key.width  = unit(0.50, "cm"),
                           legend.spacing.y  = unit(0.02, "cm")) +
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
leg_grob2 <- {
  g <- ggplot2::ggplotGrob(p_cv2_leg)
  g$grobs[[which(vapply(g$grobs, function(x) x$name, "") == "guide-box")]]
}
row1b <- wrap_plots(p_cv2_list[c("APO", "AZD5305", "olaparib", "talazoparib")], ncol = 4)
row2b <- wrap_plots(c(p_cv2_list[c("niraparib", "rucaparib", "veliparib")], list(leg_grob2)), ncol = 4)
panel_grid2 <- wrap_plots(row1b, row2b, ncol = 1)
cairo_pdf(file.path(out_dir, "Fig_S2_pmf_hd_art_grid.pdf"), width = 178/25.4, height = 97/25.4, pointsize = 8)
print(panel_grid2)
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
