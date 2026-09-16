#!/usr/bin/env Rscript
# 26-qsar_analysis.R — Structure-dynamics QSAR correlation
#
# Correlate inhibitor physicochemical properties with PMF-derived
# conformational features (well depth, barrier height, RMSF, etc.)
#
# Uses RDKit-computed molecular descriptors from SMILES
# Output: results/figures/Fig_QSAR.pdf
#         results/analysis/qsar_correlation.csv

library(data.table)
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

# ---- Inhibitor properties (canonical source: seed dataset, RDKit descriptors) ----
six <- c("Talazoparib", "Olaparib", "Niraparib", "Rucaparib", "Veliparib", "AZD5305")
seed <- fread("data/01_curated/trapping_seed_dataset.csv")
# AZD5305 is stored under its INN "saruparib" (synonyms column: AZD5305, CID 155586901)
seed_key <- ifelse(tolower(six) == "azd5305", "saruparib", tolower(six))
seed6 <- seed[match(seed_key, tolower(seed$name))]
stopifnot(all(!is.na(seed6$name)))

# Trapping potency (x olaparib): canonical curated CSV; AZD5305 left-censored <0.01 (Pires 2025)
trap <- fread("data/01_curated/trapping_potency.csv")
trap_map <- setNames(trap$trapping_x_olaparib, trap$inhibitor)
trap_vals <- unname(trap_map[tolower(six)])
trap_vals[is.na(trap_vals)] <- 0.01  # AZD5305: left-censored lower bound, not in the curated table

inhibitors <- data.frame(
  name = six,
  MW = seed6$mw, logP = seed6$clogp,
  n_HBD = seed6$hbd, n_HBA = seed6$hba, TPSA = seed6$tpsa,
  n_rotatable = seed6$rotatable_bonds, n_rings = seed6$n_rings,
  n_heavy = seed6$n_heavy,
  trapping_potency = trap_vals,
  log10_trapping = log10(trap_vals),
  stringsAsFactors = FALSE
)

# ---- PMF features (loaded from canonical analysis outputs) ------------------
mech <- fread("results/figures/Fig_Mechanism_Data.csv")          # wd_S1/wd_S2/wd_ratio/rcmin_S1
feat <- fread("results/analysis/pmf_features_summary.csv")       # S1 barrier (12-descriptive_analysis.R)
rmsf <- fread("results/analysis/rmsf_recomp/s2_rmsf_dccm_uniform.csv")  # RMSF + DCCM
pcs  <- fread("results/analysis/pca_struct_system_stats.csv")    # structural-PCA total variance

mech6 <- mech[match(tolower(six), tolower(mech$ligand))]
feat6 <- feat[system == "S1"][match(tolower(six), tolower(feat$ligand))]
rmsf6 <- rmsf[match(tolower(six), tolower(rmsf$system))]
pcs6  <- pcs[match(six, pcs$System)]
stopifnot(all(!is.na(mech6$ligand)), all(!is.na(rmsf6$system)), all(!is.na(pcs6$System)))

pmf_features <- data.frame(
  name = six,
  S1_well_depth   = mech6$wd_S1,
  S2_well_depth   = mech6$wd_S2,
  S1_S2_ratio     = mech6$wd_ratio,
  S1_min_position = mech6$rcmin_S1,
  S1_barrier      = feat6$barrier,
  S2_mean_RMSF    = rowMeans(cbind(rmsf6$HD_rmsf, rmsf6$ART_rmsf)),
  hd_art_corr     = rmsf6$dccm_abs,   # DCCM mean |r|
  s2_pca_variance = pcs6$total_var,   # structural PCA (25-pca_analysis.R); internal only, not manuscript-cited
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
            mean(df$S1_well_depth[df$name %in% c("Talazoparib", "Olaparib")])))
cat(sprintf("  PCA_variance=%.0f (mean=%.0f)\n",
            azd_out$s2_pca_variance, 
            mean(df$s2_pca_variance[df$name != "AZD5305"])))

message("===== QSAR analysis complete =====")
