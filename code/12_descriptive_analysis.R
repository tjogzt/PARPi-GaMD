#!/usr/bin/env Rscript
# 12_descriptive_analysis.R — Descriptive comparison & AZD5305 anomaly analysis
#
# Purpose:  Descriptor-level descriptive analysis of the seven-system panel (PMF feature
#           comparison + S1/S2 ratio panel); source of the descriptive figures/tables.
# Created:  2026-09-15 (header standardised 2026-10-05; R6 audit extension 2026-10-10)
# Inputs:   data/cumulant_convergence.csv ; data/01_curated/manuscript_anchor_values.csv ;
#           data/analysis/sys*_pmf_c3.xvg ; data/analysis/r6/r6_panel_values.csv
# Outputs:  figures/pdf/Fig_Feature_Comparison.pdf ; figures/pdf/Fig_S1S2_Ratio.pdf ;
#           data/analysis/pmf_features_summary.csv ; data/analysis/S1_S2_ratios.csv
# Depends:  data.table, dplyr, patchwork, ggrepel, common/helpers.R
# Run:      Rscript code/12_descriptive_analysis.R   (from the repository root)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)

# Shared helpers: read_pmf / theme_7pt / extract_features (single source).
args_h <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", args_h))) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
  source(file.path(script_dir, "..", "common", "helpers.R"))
} else {
  source("common/helpers.R")
}


# Shared helpers: theme_7pt / extract_features / read_pmf (single source).
args_h <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", args_h))) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
  source(file.path(script_dir, "..", "common", "helpers.R"))
} else {
  source("common/helpers.R")
}

dir.create("figures/pdf", showWarnings = FALSE, recursive = TRUE)

# ---- 1. PMF feature extraction ---------------------------------------------


# extract_features comes from common/helpers.R (single source).


# ---- 2. Load all systems ----------------------------------------------------
systems <- list(
  # S1 (CAT-only): HD-ART PMF
  S1_APO         = list(path = "data/analysis/sys1_APO_pmf_c3.xvg",        ligand = "APO",         system = "S1"),
  S1_AZD5305     = list(path = "data/analysis/sys1_AZD5305_pmf_c3.xvg",    ligand = "AZD5305",     system = "S1"),
  S1_niraparib   = list(path = "data/analysis/sys1_niraparib_pmf_c3.xvg",  ligand = "Niraparib",   system = "S1"),
  S1_olaparib    = list(path = "data/analysis/sys1_olaparib_pmf_c3.xvg",   ligand = "Olaparib",    system = "S1"),
  S1_rucaparib   = list(path = "data/analysis/sys1_rucaparib_pmf_c3.xvg",  ligand = "Rucaparib",   system = "S1"),
  S1_talazoparib = list(path = "data/analysis/sys1_talazoparib_pmf_c3.xvg",ligand = "Talazoparib", system = "S1"),
  S1_veliparib   = list(path = "data/analysis/sys1_veliparib_pmf_c3.xvg",  ligand = "Veliparib",   system = "S1"),
  # S2 (DNA-bound): HD-ART PMF
  S2_APO         = list(path = "data/analysis/pmf-c3-sys2_APO_CV2_cv.dat.xvg",         ligand = "APO",        system = "S2"),
  S2_AZD5305     = list(path = "data/analysis/pmf-c3-sys2_AZD5305_CV2_cv.dat.xvg",     ligand = "AZD5305",    system = "S2"),
  S2_niraparib   = list(path = "data/analysis/pmf-c3-sys2_niraparib_CV2_cv.dat.xvg",   ligand = "Niraparib",  system = "S2"),
  S2_olaparib    = list(path = "data/analysis/pmf-c3-sys2_olaparib_CV2_cv.dat.xvg",    ligand = "Olaparib",   system = "S2"),
  S2_rucaparib   = list(path = "data/analysis/pmf-c3-sys2_rucaparib_CV2_cv.dat.xvg",   ligand = "Rucaparib",  system = "S2"),
  S2_talazoparib = list(path = "data/analysis/pmf-c3-sys2_talazoparib_CV2_cv.dat.xvg", ligand = "Talazoparib",system = "S2"),
  S2_veliparib   = list(path = "data/analysis/pmf-c3-sys2_veliparib_CV2_cv.dat.xvg",   ligand = "Veliparib",  system = "S2")
)

features <- lapply(systems, function(s) {
  df <- read_pmf(s$path)
  feat <- as.list(extract_features(df))
  feat$ligand <- s$ligand
  feat$system <- s$system
  feat$n_points <- nrow(df)
  feat
})
features_df <- bind_rows(features)

# ---- R6 audit: well_depth = occupied-support C3 spans (Table S1/S5 convention) ----
# Replaces raw min-normalized curve maxima; rc_min / rc_range remain curve descriptors.
panel <- read.csv("data/analysis/r6/r6_panel_values.csv")
key <- tolower(features_df$ligand)
features_df$well_depth <- ifelse(
  features_df$system == "S1",
  panel$S1_C3[match(key, panel$ligand)],
  panel$S2CV2_C3[match(key, panel$ligand)])
stopifnot(all(is.finite(features_df$well_depth)))

# ---- 3. Feature table ------------------------------------------------------
print("=== PMF Feature Summary ===")
print(features_df, n = 50)
write.csv(features_df, "data/analysis/pmf_features_summary.csv", row.names = FALSE)

# ---- 4. S1 vs S2 feature scatter -------------------------------------------
p_well_depth <- ggplot(features_df, aes(x = ligand, y = well_depth, fill = system)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("S1" = "#2166AC", "S2" = "#B2182B")) +
  labs(x = NULL, y = "Sampled-range span (kcal/mol)", title = "Sampled-range span") +
  theme_7pt + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_rc_min <- ggplot(features_df, aes(x = ligand, y = rc_min, fill = system)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("S1" = "#2166AC", "S2" = "#B2182B")) +
  labs(x = NULL, y = "RC at Minimum (Å)", title = "Equilibrium distance") +
  theme_7pt + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_rc_range <- ggplot(features_df, aes(x = ligand, y = rc_range, fill = system)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("S1" = "#2166AC", "S2" = "#B2182B")) +
  labs(x = NULL, y = "RC Range (Å)", title = "Sampling range") +
  theme_7pt + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_features <- (p_well_depth | p_rc_min | p_rc_range) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        # (2026-10-15) full-width layout: sizes tuned for 178-mm canvas
        # displayed at \textwidth in the SI
        axis.text    = element_text(size = 9),
        axis.title   = element_text(size = 9),
        plot.title   = element_text(size = 10, face = "bold"),
        legend.text  = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"))

cairo_pdf("figures/pdf/Fig_Feature_Comparison.pdf", width = 178/25.4, height = 88/25.4, pointsize = 9)
print(p_features)
dev.off()

# ---- 5. AZD5305 anomaly analysis ------------------------------------------
cat("\n========== AZD5305 Anomaly Analysis ==========\n")

s1_others <- features_df %>% filter(system == "S1", ligand != "AZD5305")
s1_azd    <- features_df %>% filter(system == "S1", ligand == "AZD5305")

for (col in c("well_depth", "rc_min", "rc_range", "barrier")) {
  mu   <- mean(s1_others[[col]], na.rm = TRUE)
  sd_v <- sd(s1_others[[col]], na.rm = TRUE)
  azd_v <- s1_azd[[col]]
  z <- (azd_v - mu) / sd_v
  cat(sprintf("  %s: AZD5305 = %.2f, Others (mean±SD) = %.2f±%.2f, Z = %.2f\n",
              col, azd_v, mu, sd_v, z))
}

# Compare S1 veliparib (extreme well depth)
cat("\nVeliparib S1 well_depth =", 
    features_df$well_depth[features_df$system == "S1" & features_df$ligand == "Veliparib"],
    "kcal/mol (2× any other)\n")

# ---- 6. S1 inhibitor ranking ---
cat("\n========== S1 HD-ART Well Depth Ranking ==========\n")
s1_rank <- features_df %>%
  filter(system == "S1") %>%
  arrange(desc(well_depth)) %>%
  select(ligand, well_depth, rc_min, rc_range)
print(s1_rank, n = 10)

cat("\n========== S2 HD-ART Well Depth Ranking ==========\n")
s2_rank <- features_df %>%
  filter(system == "S2") %>%
  arrange(desc(well_depth)) %>%
  select(ligand, well_depth, rc_min, rc_range)
print(s2_rank, n = 10)

# ---- 7. Ratio analysis: S1/S2 fold change ---
cat("\n========== S1/S2 Ratio (DNA-free amplification) ==========\n")
ratio_df <- features_df %>%
  select(ligand, system, well_depth, rc_range, barrier) %>%
  pivot_wider(names_from = system, values_from = c(well_depth, rc_range, barrier)) %>%
  mutate(
    well_depth_ratio = well_depth_S1 / well_depth_S2,
    rc_range_ratio   = rc_range_S1 / rc_range_S2
  ) %>%
  filter(!is.na(well_depth_ratio)) %>%
  arrange(desc(well_depth_ratio))
print(ratio_df, n = 10, width = Inf)

write.csv(ratio_df, "data/analysis/S1_S2_ratios.csv", row.names = FALSE)

# ---- 8. S1 vs S2 span scatter (audited C3 convention) ---
# Seven base conditions (filled) + the two protonated re-simulations (open); the
# dashed line marks equality (AAI = 1). Values from the audited panel table.
pv2 <- read.csv("data/analysis/r6/r6_panel_values.csv")
base <- pv2[!grepl("_prot", pv2$ligand), ]
prot <- pv2[grepl("_prot", pv2$ligand), ]
sc <- rbind(
  data.frame(x = base$S2CV2_C3, y = base$S1_C3, label = base$label, filled = TRUE),
  data.frame(x = prot$S2CV2_C3, y = prot$S1_C3, label = prot$label, filled = FALSE))
sc$fill_key <- ifelse(sc$filled, "base", "prot")
p_ratio <- ggplot(sc, aes(x = x, y = y, label = label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = label, shape = fill_key), size = 3.0, stroke = 0.9) +
  scale_shape_manual(values = c("base" = 16, "prot" = 1), guide = "none") +
  ggrepel::geom_text_repel(size = 3.0, force = 2, box.padding = 0.35, max.overlaps = Inf,
                          seed = 49, family = "Arial") +
  scale_color_manual(values = c("APO" = "grey40", "AZD5305" = "darkorange",
                                 "Olaparib" = "#E41A1C", "Talazoparib" = "#FF7F00",
                                 "Veliparib" = "#377EB8", "Niraparib" = "#4DAF4A",
                                 "Rucaparib" = "#984EA3", "Fluzoparib" = "#C23531",
                                 "Pamiparib" = "#3D6BA8", "Senaparib" = "#9D2933",
                                 "Niraparib (prot.)" = "#4DAF4A",
                                 "Rucaparib (prot.)" = "#984EA3")) +
  labs(x = "S2 CV2 span (kcal/mol, C3)", 
       y = "S1 span (kcal/mol, C3)") +
  theme_7pt + theme(legend.position = "none",
                  plot.margin = ggplot2::margin(1, 1, 1, 3, unit = "mm"))

cairo_pdf("figures/pdf/Fig_S1S2_Ratio.pdf", width = 132/25.4, height = 115.5/25.4, pointsize = 8)
print(p_ratio)
dev.off()

message("===== Descriptive analysis complete =====")
