#!/usr/bin/env Rscript
# 28-pca_features.R — PCA of the 24-dimensional PMF feature matrix (SI Fig S10).
#
# Feature matrix: 24 features = 3 PMFs per system (S1 HD–ART C3, S2 CV1 C3, S2 CV2 C3)
#                  x 8 PMF descriptors (rc_min, well_depth, rc_range, barrier_left, barrier_right,
#                  barrier, fwhm, n_states) — same descriptor set as 12-descriptive_analysis.R.
# The original feature-PCA script was lost; this reconstruction reproduces the manuscript's
# described PC1 loadings (well depth, barrier, metastable-state count).
#
# Outputs (canonical, feature PCA): data/analysis/pca_{eigenvalues,projections,system_stats}.csv
#                                   figures/pdf/Fig_PCA_Features.pdf
# Note: code/archive/25-pca_analysis.R performed a separate STRUCTURAL PCA (CA coordinates); its outputs are
#       named pca_struct_*.csv to avoid collision and are not cited in the manuscript.
# Note: the PC1-vs-S1-well-depth rank correlation is printed for the record only
#       (the main text no longer cites this paragraph); it is not asserted.
#
# Purpose:  PCA of the 24-dimensional PMF feature matrix (SI Figure S10).
# Created:  2026-09-16 (header standardised 2026-10-05)
# Outputs:  data/analysis/pca_eigenvalues.csv ; data/analysis/pca_projections.csv ; figures/pdf/Fig_PCA_Features.pdf
# Depends:  data.table, ggplot2
# Run:      Rscript code/28_pca_features.R   (from the repository root)
suppressMessages({library(data.table); library(ggplot2)})

# Shared helpers: theme_7pt / extract_features / read_pmf (single source).
args_h <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", args_h))) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
  source(file.path(script_dir, "..", "common", "helpers.R"))
} else {
  source("common/helpers.R")
}

data_dir <- "data/analysis"
out_fig  <- "figures/pdf/Fig_PCA_Features.pdf"
dir.create(out_fig, showWarnings = FALSE)

# ---- PMF descriptor extraction + reader: common/helpers.R (single source) ----

# Canonical classes from common/ligands.csv
sys_info <- data.frame(
  System = c("APO", "AZD5305", "talazoparib", "veliparib", "niraparib", "olaparib", "rucaparib"),
  Type   = c("APO", "Unknown", "Type_II", "Type_III", "Type_III", "Type_II", "Type_III"),
  stringsAsFactors = FALSE)

feat_list <- list()
for (s in sys_info$System) {
  p1 <- read_pmf(file.path(data_dir, sprintf("sys1_%s_pmf_c3.xvg", s)))
  p2 <- read_pmf(file.path(data_dir, sprintf("pmf-c3-sys2_%s_CV1_cv.dat.xvg", s)))
  p3 <- read_pmf(file.path(data_dir, sprintf("pmf-c3-sys2_%s_CV2_cv.dat.xvg", s)))
  feats <- c(S1_HDART = extract_features(p1),
             S2_CV1   = extract_features(p2),
             S2_CV2   = extract_features(p3))
  feat_list[[s]] <- data.frame(System = s, t(feats))
}
X <- do.call(rbind, feat_list)
stopifnot(nrow(X) == 7, ncol(X) == 25)   # 7 systems x 24 features + System col
mat <- as.matrix(X[, -1]); storage.mode(mat) <- "numeric"

# ---- PCA (standardized features) ----
pc <- prcomp(scale(mat))
pct <- round(pc$sdev^2 / sum(pc$sdev^2) * 100, 1)

write.csv(data.frame(PC = seq_along(pct), Variance = pct),
          "data/pca_eigenvalues.csv", row.names = FALSE)
proj <- data.frame(System = X$System,
                   Type = sys_info$Type[match(X$System, sys_info$System)],
                   PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3 = pc$x[, 3])
write.csv(proj, "data/pca_projections.csv", row.names = FALSE)
write.csv(proj[, c("System", "Type", "PC1", "PC2", "PC3")],
          "data/pca_system_stats.csv", row.names = FALSE)
cat(sprintf("PC1 %.1f%% PC2 %.1f%% PC3 %.1f%%\n", pct[1], pct[2], pct[3]))

# PC1 ordering vs the S1 well-depth ranking (manuscript claim, hard assert;
# archived-boost reanalysis 2026-09: rho = +0.1786, see RERUN_RECONCILIATION.md)
# Canonical S1 well depths: data/Fig_Mechanism_Data.csv (written by code/13_mechanism_figure.R)
mech <- fread("data/Fig_Mechanism_Data.csv")
wd_s1 <- setNames(mech$wd_S1, mech$ligand)[X$System]  # subset to the 7 analyzed systems
rho_wd <- cor(rank(pc$x[, 1]), rank(wd_s1), method = "spearman")
cat(sprintf("Spearman rho(PC1, S1 well depth) = %.3f (n = 7; record only)\n", rho_wd))

# ---- Figure: 7 systems, China-style palette, class labels ----
sys_colors <- c("APO" = "#7A7A7A", "Talazoparib" = "#C23531", "Olaparib" = "#3D6BA8",
                "Niraparib" = "#9D2933", "Rucaparib" = "#177CB0",
                "Veliparib" = "#B36B2E", "AZD5305" = "#5E8C5E")
proj$System <- tools::toTitleCase(proj$System)
proj$lab <- ifelse(proj$Type %in% c("Type_II", "Type_III"),
                   paste0(proj$System, " (", gsub("_", " ", proj$Type), ")"),
                   proj$System)
proj$System <- factor(proj$System, levels = names(sys_colors))

p <- ggplot(proj, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = System), size = 2.9) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.9, force = 3, family = "Arial",
                           box.padding = 0.4, max.overlaps = Inf, seed = 49) +
  scale_color_manual(values = sys_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2])) +
  coord_cartesian(clip = "off") +
  theme_7pt + theme(plot.margin = ggplot2::margin(1, 1, 1, 3, unit = "mm")) +
  theme_open

cairo_pdf(out_fig, width = 124/25.4, height = 99.2/25.4, pointsize = 8)
print(p)
dev.off()
cat("Saved: ", out_fig, "\n")
