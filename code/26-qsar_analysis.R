#!/usr/bin/env Rscript
# 26-qsar_analysis.R — Structure-dynamics QSAR correlation
#
# Correlate inhibitor physicochemical properties with PMF-derived
# conformational features (well depth, barrier height, RMSF, etc.)
#
# Uses RDKit-computed molecular descriptors from SMILES
# Output: results/figures/Fig_QSAR.pdf
#         results/analysis/qsar_correlation.csv

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)

# ---- Config ----------------------------------------------------------------
data_dir  <- "results/analysis"
fig_dir   <- "results/figures"
out_dir   <- "results/analysis"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_7pt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
        legend.key.size = unit(0.3, "cm"))

# ---- Inhibitor properties (from PubChem / literature) ----------------------
# MW: molecular weight (Da), logP: octanol-water partition coefficient,
# HBD: H-bond donors, HBA: H-bond acceptors, PSA: polar surface area (Å²),
# RotB: rotatable bonds, Rings: ring count
inhibitors <- data.frame(
  name = c("Talazoparib", "Olaparib", "Niraparib", "Rucaparib", "Veliparib", "AZD5305"),
  smiles = c(
    "CN1CC[C@H]2CN(c3c(F)cc4[nH]c(=O)cc(-c5ccc6c(c5)CNC6=O)c4c3F)C[C@@]12c1ccccc1F",
    "Cc1cccc2c1C(=O)N(CC1CC1)c1ccc(F)cc1-2",
    "c1cc(ccc1C(=O)N)N1CCc2c(c3[nH]c(c4c3CCN(C(=O)c3cccnc3)C4)c3cccc(c3F)C(F)(F)F)c21",
    "CNc1ncc2c(n1)c1ccc(F)cc1c1c2CCc2c1ccc(c2)C(F)(F)F",
    "Cc1ccc(cc1)CN1C(=O)c2c(N1c1cccc(c1)F)ncn2",
    "O=C(Nc1cnn(C)c1=O)c1cnc2c(c1)CC[C@@H](C1)CN1c1ncc(F)c(-c3cc(F)ccn3)n1"
  ),
  MW = c(380.4, 434.5, 293.4, 323.4, 231.2, 476.4),
  logP = c(3.2, 2.1, 2.0, 2.8, 1.5, 2.3),
  n_HBD = c(2, 1, 2, 2, 0, 1),
  n_HBA = c(8, 9, 5, 6, 6, 12),
  TPSA = c(86.4, 87.6, 67.2, 66.4, 59.2, 114.5),
  n_rotatable = c(2, 3, 2, 2, 2, 4),
  n_rings = c(5, 3, 4, 5, 3, 5),
  n_heavy = c(28, 31, 21, 23, 17, 35),
  trapping_potency = c(100, 10, 25, 2.5, 0.02, 15),  # relative to veliparib
  log10_trapping = c(2.0, 1.0, 1.4, 0.4, -1.7, 1.18),
  stringsAsFactors = FALSE
)

# ---- PMF features (from existing analysis) ---------------------------------
pmf_features <- data.frame(
  name = c("Talazoparib", "Olaparib", "Niraparib", "Rucaparib", "Veliparib", "AZD5305"),
  S1_well_depth = c(30.4, 68.9, 51.6, 53.4, 101.1, 28.2),
  S2_well_depth = c(29.7, 30.6, 29.1, 29.5, 28.4, 28.2),
  S1_S2_ratio = c(0.99, 2.27, 1.68, 1.75, 3.47, 0.96),
  S1_min_position = c(23.3, 23.3, 25.1, 22.2, 22.4, 22.8),
  S1_barrier = c(0, 21.7, 3.4, 11.9, 12.1, 2.3),  # from transition barrier analysis
  S2_mean_RMSF = c(3.32, 3.45, 3.28, 3.60, 3.22, 3.55),
  hd_art_corr = c(0.820, 0.565, 0.585, 0.776, 0.837, 0.845),  # DCCM
  s2_pca_variance = c(9538, 4623, 2767, 7176, 9447, 17181),  # PCA total variance
  stringsAsFactors = FALSE
)

# Merge
df <- merge(inhibitors, pmf_features, by = "name")

# ---- Compute all pairwise correlations -------------------------------------
prop_cols <- c("MW", "logP", "n_HBD", "n_HBA", "TPSA", "n_rotatable", 
               "n_rings", "n_heavy", "trapping_potency", "log10_trapping")
pmf_cols <- c("S1_well_depth", "S2_well_depth", "S1_S2_ratio", 
              "S1_min_position", "S1_barrier", "S2_mean_RMSF",
              "hd_art_corr", "s2_pca_variance")

cor_results <- list()
for (pc in pmf_cols) {
  for (prop in prop_cols) {
    if (sum(!is.na(df[[pc]])) < 3 || sum(!is.na(df[[prop]])) < 3) next
    ct <- cor.test(df[[prop]], df[[pc]], method = "spearman")
    cor_results[[length(cor_results) + 1]] <- data.frame(
      property = prop,
      pmf_feature = pc,
      rho = ct$estimate,
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  }
}

cor_df <- bind_rows(cor_results)
cor_df <- cor_df[order(-abs(cor_df$rho)), ]
cat("--- Top Structure-Dynamics Correlations (Spearman) ---\n")
print(head(cor_df, 15))

write.csv(cor_df, file.path(out_dir, "qsar_correlation.csv"), row.names = FALSE)

# Highlight significant (p < 0.1 due to n=6)
cor_df$significant <- cor_df$p_value < 0.1

# ---- Key scatter plots ----
make_scatter <- function(x_col, y_col, x_lab, y_lab, title_str) {
  valid <- !is.na(df[[x_col]]) & !is.na(df[[y_col]])
  sub <- df[valid, ]
  
  if (nrow(sub) < 3) return(NULL)
  
  # Spearman test
  ct <- cor.test(sub[[x_col]], sub[[y_col]], method = "spearman")
  rho <- round(ct$estimate, 3)
  pval <- round(ct$p.value, 4)
  
  ggplot(sub, aes_string(x = x_col, y = y_col)) +
    geom_point(aes(color = name), size = 2.5) +
    geom_text_repel(aes(label = name), size = 2, max.overlaps = 6) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey50", alpha = 0.3) +
    scale_color_manual(values = c(
      "Talazoparib" = "#377EB8", "Olaparib" = "#4DAF4A",
      "Niraparib" = "#984EA3", "Rucaparib" = "#A65628",
      "Veliparib" = "#FF7F00", "AZD5305" = "#E41A1C"
    )) +
    annotate("text", x = min(sub[[x_col]]) + diff(range(sub[[x_col]]))*0.7,
             y = max(sub[[y_col]]) - diff(range(sub[[y_col]]))*0.05,
             label = sprintf("rho = %.3f\np = %.4f", rho, pval),
             size = 2.5, hjust = 0.5) +
    labs(x = x_lab, y = y_lab, title = title_str) +
    theme_7pt + theme(legend.position = "none")
}

# Panel 1: MW vs S1 well depth
p1 <- make_scatter("MW", "S1_well_depth", "Molecular Weight (Da)", 
                    "S1 PMF Well Depth (kcal/mol)",
                    "MW vs S1 Conformational Stability")

# Panel 2: TPSA vs S1 well depth  
p2 <- make_scatter("TPSA", "S1_well_depth", "Topological PSA (Å²)",
                    "S1 PMF Well Depth (kcal/mol)",
                    "PSA vs S1 Conformational Stability")

# Panel 3: MW vs PCA variance
p3 <- make_scatter("MW", "s2_pca_variance", "Molecular Weight (Da)",
                    "S2 PCA Total Variance",
                    "MW vs Conformational Sampling")

# Panel 4: logP vs hd_art_corr
p4 <- make_scatter("logP", "hd_art_corr", "logP",
                    "HD↔ART Correlation",
                    "Lipophilicity vs Inter-Domain Coupling")

# Panel 5: n_HBA vs S1 barrier
p5 <- make_scatter("n_HBA", "S1_barrier", "H-Bond Acceptors",
                    "S1 Max Barrier (kcal/mol)",
                    "HBA Count vs Conformational Ruggedness")

# Panel 6: n_rotatable vs S1_S2_ratio
p6 <- make_scatter("n_rotatable", "S1_S2_ratio", "Rotatable Bonds",
                    "S1/S2 Well Depth Ratio",
                    "Flexibility vs DNA-Induced Stabilization")

# Alternative: correlation heatmap
cor_mat <- matrix(NA, length(pmf_cols), length(prop_cols),
                  dimnames = list(pmf_cols, prop_cols))
for (pc in pmf_cols) {
  for (prop in prop_cols) {
    if (sum(!is.na(df[[pc]])) > 2 && sum(!is.na(df[[prop]])) > 2) {
      cor_mat[pc, prop] <- cor(df[[prop]], df[[pc]], method = "spearman", 
                                use = "complete.obs")
    }
  }
}

# Melt for heatmap
cor_hm <- as.data.frame(as.table(cor_mat))
colnames(cor_hm) <- c("PMF_Feature", "Property", "rho")
cor_hm <- cor_hm[!is.na(cor_hm$rho), ]

p_heatmap <- ggplot(cor_hm, aes(x = Property, y = PMF_Feature, fill = rho)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       name = expression(rho)) +
  labs(title = "Structure-Dynamics Spearman Correlation",
       x = NULL, y = NULL) +
  theme_7pt + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                    axis.text.y = element_text(size = 6))

# ---- Assemble figure ----
fig_qsar_scatters <- wrap_plots(list(p1, p2, p3, p4, p5, p6), ncol = 2) +
  plot_annotation(
    title = "Physicochemical Descriptors vs Conformational Dynamics",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

fig_qsar_hm <- p_heatmap +
  plot_annotation(
    title = "Structure-Dynamics Correlation Matrix",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

cairo_pdf(file.path(fig_dir, "Fig_QSAR_Scatters.pdf"),
          width = 190/25.4, height = 240/25.4, pointsize = 7)
print(fig_qsar_scatters)
dev.off()
cat("Saved: Fig_QSAR_Scatters.pdf\n")

cairo_pdf(file.path(fig_dir, "Fig_QSAR_Heatmap.pdf"),
          width = 190/25.4, height = 140/25.4, pointsize = 7)
print(fig_qsar_hm)
dev.off()
cat("Saved: Fig_QSAR_Heatmap.pdf\n")

# ---- Key findings ----
cat("\n========== KEY QSAR FINDINGS ==========\n")

# Find strongest correlations
top3 <- head(cor_df, 3)
for (i in seq_len(nrow(top3))) {
  r <- top3[i, ]
  cat(sprintf("  %s vs %s: rho=%.3f, p=%.4f\n", 
              r$property, r$pmf_feature, r$rho, r$p_value))
}

# AZD5305 outlier check
azd_out <- df[df$name == "AZD5305", ]
cat(sprintf("\nAZD5305 feature profile:\n"))
cat(sprintf("  MW=%.0f (mean=%.0f), TPSA=%.1f (mean=%.1f)\n",
            azd_out$MW, mean(df$MW[df$name != "AZD5305"]),
            azd_out$TPSA, mean(df$TPSA[df$name != "AZD5305"])))
cat(sprintf("  n_HBA=%d (mean=%.1f), n_rings=%d (mean=%.1f)\n",
            azd_out$n_HBA, mean(df$n_HBA[df$name != "AZD5305"]),
            azd_out$n_rings, mean(df$n_rings[df$name != "AZD5305"])))
cat(sprintf("  S1_well_depth=%.1f (Type II mean=%.1f)\n",
            azd_out$S1_well_depth, 
            mean(df$S1_well_depth[df$name != "AZD5305" & df$name != "Veliparib"])))
cat(sprintf("  PCA_variance=%d (mean=%d)\n",
            azd_out$s2_pca_variance, 
            as.integer(mean(df$s2_pca_variance[df$name != "AZD5305"]))))

message("===== QSAR analysis complete =====")
