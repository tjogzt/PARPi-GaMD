#!/usr/bin/env Rscript
# 25-pca_analysis.R — STRUCTURAL PCA (CA cartesian coordinates) of S2 trajectories
# NOTE: this is NOT the manuscript's "24-dimensional PMF feature matrix" PCA (that analysis is
# code/28-pca_features.R). This script's outputs are named pca_struct_*.csv and are not cited.

library(bio3d)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Config ----------------------------------------------------------------
data_root <- function() {
  r <- Sys.getenv("DATA_ROOT", unset = "")
  if (!nzchar(r)) stop("Set DATA_ROOT env var to the trajectory data directory")
  r
}

data_dir  <- "results/analysis"
dcd_dir   <- file.path(data_root(), "md_analysis")
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_7pt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
        legend.key.size = unit(0.3, "cm"))

# System definitions
SYSTEMS <- c(
  "sys2_APO" = "APO",
  "sys2_AZD5305" = "AZD5305",
  "sys2_talazoparib" = "Talazoparib",
  "sys2_veliparib" = "Veliparib",
  "sys2_niraparib" = "Niraparib",
  "sys2_olaparib" = "Olaparib",
  "sys2_rucaparib" = "Rucaparib"
)

TYPE <- c(
  "APO" = "APO", "AZD5305" = "Unknown", "Talazoparib" = "Type_II",
  "Veliparib" = "Type_III", "Niraparib" = "Type_III",
  "Olaparib" = "Type_II", "Rucaparib" = "Type_III"
)

COLORS <- c(
  "APO" = "#999999", "AZD5305" = "#E41A1C",
  "Talazoparib" = "#377EB8", "Veliparib" = "#FF7F00",
  "Niraparib" = "#984EA3", "Olaparib" = "#4DAF4A",
  "Rucaparib" = "#A65628"
)

# ---- Load reference and all CA trajectories --------------------------------
ref_pdb <- read.pdb(file.path(data_dir, "ref_ca.pdb"))
ca_inds <- atom.select(ref_pdb, elety = "CA")

# Load and concatenate trajectories (subsampled for memory)
cat("Loading CA trajectories for PCA...\n")
all_xyz <- list()
all_labels <- c()
max_frames_per_sys <- 200  # subsample to keep memory manageable

for (sys_name in names(SYSTEMS)) {
  ca_dcd <- file.path(dcd_dir, paste0(sys_name, "_ca.dcd"))
  if (!file.exists(ca_dcd)) {
    cat(sprintf("[SKIP] %s\n", sys_name))
    next
  }
  
  label <- SYSTEMS[sys_name]
  cat(sprintf("  Reading %s...\n", label))
  
  dcd <- read.dcd(ca_dcd)
  n_total <- dim(dcd)[1]
  
  # Subsample
  if (n_total > max_frames_per_sys) {
    stride <- max(floor(n_total / max_frames_per_sys), 1)
    idx <- seq(1, n_total, by = stride)
    if (length(idx) > max_frames_per_sys) idx <- idx[1:max_frames_per_sys]
  } else {
    idx <- seq_len(n_total)
  }
  
  cat(sprintf("    %d frames (of %d)\n", length(idx), n_total))
  
  for (i in seq_along(idx)) {
    all_xyz[[length(all_xyz) + 1]] <- dcd[idx[i], ]
    all_labels <- c(all_labels, label)
  }
}

cat(sprintf("Total frames for PCA: %d\n", length(all_xyz)))

# Convert to matrix
xyz_mat <- do.call(rbind, all_xyz)
cat(sprintf("XYZ matrix: %d x %d\n", nrow(xyz_mat), ncol(xyz_mat)))

# ---- PCA on concatenated ensemble ------------------------------------------
cat("Performing PCA...\n")
pca <- pca.xyz(xyz_mat)
pca_var_pct <- pca$L / sum(pca$L) * 100
cat(sprintf("PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%\n",
            pca_var_pct[1], pca_var_pct[2], pca_var_pct[3]))

# Save full eigenvalues for scree
write.csv(data.frame(PC = seq_along(pca_var_pct), Variance = pca_var_pct),
          file.path(data_dir, "pca_struct_eigenvalues.csv"), row.names = FALSE)

# Project each system's frames onto PC1 and PC2
proj_df <- data.frame(
  PC1 = pca$z[, 1],
  PC2 = pca$z[, 2],
  PC3 = pca$z[, 3],
  System = factor(all_labels, levels = names(COLORS)),
  Type = TYPE[all_labels]
)

write.csv(proj_df, file.path(data_dir, "pca_struct_projections.csv"), row.names = FALSE)

# ---- Compute per-system 2D free energy landscapes --------------------------
compute_fel <- function(pc1, pc2, n_bins = 30, kT = 0.596) {
  # 2D histogram → free energy
  h2d <- matrix(0, n_bins, n_bins)
  
  x_range <- range(pc1)
  y_range <- range(pc2)
  
  x_breaks <- seq(x_range[1], x_range[2], length.out = n_bins + 1)
  y_breaks <- seq(y_range[1], y_range[2], length.out = n_bins + 1)
  
  for (i in seq_along(pc1)) {
    xi <- findInterval(pc1[i], x_breaks, all.inside = TRUE)
    yi <- findInterval(pc2[i], y_breaks, all.inside = TRUE)
    h2d[xi, yi] <- h2d[xi, yi] + 1
  }
  
  # Convert to free energy: F = -kT ln(P)
  h2d[h2d == 0] <- NA
  FEL <- -kT * log(h2d / sum(h2d, na.rm = TRUE))
  FEL <- FEL - min(FEL, na.rm = TRUE)  # zero baseline
  
  list(
    FEL = FEL,
    x_centers = (x_breaks[-1] + x_breaks[-length(x_breaks)]) / 2,
    y_centers = (y_breaks[-1] + y_breaks[-length(y_breaks)]) / 2
  )
}

# ---- FEL heatmap for each system ----
fel_plots <- list()

for (sys_label in names(COLORS)) {
  sys_data <- proj_df[proj_df$System == sys_label, ]
  if (nrow(sys_data) < 10) next
  
  # Use global PC ranges for consistency
  fel <- compute_fel(sys_data$PC1, sys_data$PC2, 
                     kT = 0.596)  # kT at 300K in kcal/mol
  
  # Convert to long format for ggplot
  fel_df <- expand.grid(x = seq_len(nrow(fel$FEL)), y = seq_len(ncol(fel$FEL)))
  fel_df$F <- as.vector(fel$FEL)
  fel_df$PC1 <- fel$x_centers[fel_df$x]
  fel_df$PC2 <- fel$y_centers[fel_df$y]
  
  # Filter to show only low-energy region (within 5 kcal of minimum)
  fel_df <- fel_df[!is.na(fel_df$F), ]
  max_f <- 6  # show up to 6 kcal/mol
  
  p <- ggplot(fel_df, aes(x = PC1, y = PC2, fill = F)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("#440154", "#31688E", "#35B779", "#FDE725"),
      limits = c(0, max_f), oob = scales::squish,
      name = "F\n(kcal/mol)"
    ) +
    coord_fixed() +
    labs(title = sys_label, x = NULL, y = NULL) +
    theme_7pt + 
    theme(axis.text = element_text(size = 5),
          plot.title = element_text(size = 7, face = "bold", hjust = 0.5))
  
  fel_plots[[sys_label]] <- p
}

# Global scatter plot (PC1 vs PC2, all systems)
p_scatter <- ggplot(proj_df, aes(x = PC1, y = PC2, color = System)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = COLORS) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var_pct[2]),
       title = "PCA Projection: All S2 Systems",
       color = NULL) +
  theme_7pt + theme(legend.position = "bottom", 
                     legend.key.size = unit(0.2, "cm"),
                     legend.text = element_text(size = 5))

# ---- Per-system PC1/PC2 variance ----
sys_stats <- proj_df %>%
  group_by(System, Type) %>%
  summarise(
    mean_PC1 = mean(PC1),
    sd_PC1 = sd(PC1),
    mean_PC2 = mean(PC2),
    sd_PC2 = sd(PC2),
    total_var = var(PC1) + var(PC2),
    .groups = "drop"
  )
cat("\n--- Per-System PC Variance ---\n")
print(sys_stats)

write.csv(sys_stats, file.path(data_dir, "pca_struct_system_stats.csv"), row.names = FALSE)

# ---- Variance comparison bar plot ----
p_var <- ggplot(sys_stats, aes(x = reorder(System, total_var), 
                                y = total_var, fill = Type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("APO" = "#999999", "Type_II" = "#2166AC", "Type_III" = "#B2182B")) +
  coord_flip() +
  labs(x = NULL, y = "Total Variance (PC1+PC2)",
       title = "Conformational Sampling Breadth by System",
       fill = NULL) +
  theme_7pt + theme(legend.position = "bottom")

# ---- Scree plot ----
p_scree <- ggplot(data.frame(PC = 1:20, Var = pca_var_pct[1:20]),
                  aes(x = PC, y = Var)) +
  geom_bar(stat = "identity", fill = "#4472C4", width = 0.7) +
  geom_line(aes(y = cumsum(Var)), color = "#B2182B", linewidth = 0.5) +
  geom_point(aes(y = cumsum(Var)), color = "#B2182B", size = 1) +
  scale_x_continuous(breaks = seq(0, 20, by = 5)) +
  labs(x = "Principal Component", y = "Variance Explained (%)",
       title = "PCA Scree Plot (S2 Ensemble)") +
  theme_7pt

# ---- Assemble FEL panels ----
# 3 rows × 3 cols layout
fel_grid <- wrap_plots(fel_plots, ncol = 3) +
  plot_annotation(
    title = "2D Free Energy Landscapes in Common PC Space (S2 Systems)",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

# ---- Figure 1: Scatter + Variance ----
fig_pca1 <- (p_scatter | p_var) +
  plot_layout(widths = c(2, 1))

# ---- Figure 2: FEL panels ----
fig_pca2 <- fel_grid

# ---- Figure 3: Scree + top 3 PCs ----
fig_pca3 <- (p_scree / (p_var)) +
  plot_layout(heights = c(1, 1))

cairo_pdf(file.path(out_dir, "Fig_PCA_Scatter.pdf"),
          width = 190/25.4, height = 120/25.4, pointsize = 7)
print(fig_pca1)
dev.off()
cat("Saved: Fig_PCA_Scatter.pdf\n")

cairo_pdf(file.path(out_dir, "Fig_PCA_FEL.pdf"),
          width = 190/25.4, height = 220/25.4, pointsize = 7)
print(fig_pca2)
dev.off()
cat("Saved: Fig_PCA_FEL.pdf\n")

# ---- Key findings ----
cat("\n========== KEY PCA FINDINGS ==========\n")
cat(sprintf("PC1: %.1f%%, PC2: %.1f%%, Cumulative: %.1f%%\n",
            pca_var_pct[1], pca_var_pct[2], pca_var_pct[1] + pca_var_pct[2]))

# Type comparison
type2 <- sys_stats[sys_stats$Type == "Type_II", ]
type3 <- sys_stats[sys_stats$Type == "Type_III", ]
cat(sprintf("Type II mean variance: %.1f (n=%d)\n",
            mean(type2$total_var), nrow(type2)))
cat(sprintf("Type III mean variance: %.1f (n=%d)\n",
            mean(type3$total_var), nrow(type3)))

# AZD5305 position
azd <- sys_stats[sys_stats$System == "AZD5305", ]
cat(sprintf("AZD5305: PC1=%.1f±%.1f, PC2=%.1f±%.1f, var=%.1f\n",
            azd$mean_PC1, azd$sd_PC1, azd$mean_PC2, azd$sd_PC2, azd$total_var))

# Mean PC positions by type
cat("\nMean PC1 by type:\n")
for (tp in c("APO", "Type_II", "Type_III")) {
  vals <- sys_stats$mean_PC1[sys_stats$Type == tp]
  cat(sprintf("  %s: %.2f (n=%d)\n", tp, mean(vals), length(vals)))
}

message("===== PCA analysis complete =====")
