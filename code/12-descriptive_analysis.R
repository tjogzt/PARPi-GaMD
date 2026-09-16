#!/usr/bin/env Rscript
# 12-descriptive_analysis.R — Descriptive comparison & AZD5305 anomaly analysis
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)

theme_7pt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6), strip.text = element_text(size = 7, face = "bold"))

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ---- 1. PMF feature extraction ---------------------------------------------
read_xvg <- function(path) {
  lines <- readLines(path)
  dstart <- which(grepl("^[0-9]", lines))[1]
  df <- read.table(text = lines[dstart:length(lines)], col.names = c("RC", "PMF"))
  df$PMF_norm <- df$PMF - min(df$PMF, na.rm = TRUE)
  df
}

extract_features <- function(pmf_df) {
  rc  <- pmf_df$RC
  pmf <- pmf_df$PMF_norm
  min_idx <- which.min(pmf)
  
  tibble(
    rc_min       = rc[min_idx],
    well_depth   = max(pmf) - min(pmf),
    rc_range     = max(rc) - min(rc),
    # Barrier: max left or right of minimum
    barrier_left = if (min_idx > 1) max(pmf[1:min_idx]) - pmf[min_idx] else NA,
    barrier_right= if (min_idx < length(rc)) max(pmf[min_idx:length(rc)]) - pmf[min_idx] else NA,
    barrier      = max(barrier_left, barrier_right, na.rm = TRUE),
    # FWHM basin width
    fwhm = {
      hm <- max(pmf) / 2
      idx <- which(pmf <= hm)
      if (length(idx) > 1) rc[tail(idx,1)] - rc[idx[1]] else NA
    },
    # # of local minima > 0.5 kcal depth
    n_states = {
      d1 <- diff(pmf)
      local_mins <- which(diff(sign(d1)) == 2) + 1
      sum(pmf[local_mins] <= (pmf[local_mins] + 0.5))
    }
  )
}

# ---- 2. Load all systems ----------------------------------------------------
systems <- list(
  # S1 (CAT-only): HD-ART PMF
  S1_APO         = list(path = "results/analysis/sys1_APO_pmf_c3.xvg",        ligand = "APO",         system = "S1"),
  S1_AZD5305     = list(path = "results/analysis/sys1_AZD5305_pmf_c3.xvg",    ligand = "AZD5305",     system = "S1"),
  S1_niraparib   = list(path = "results/analysis/sys1_niraparib_pmf_c3.xvg",  ligand = "Niraparib",   system = "S1"),
  S1_olaparib    = list(path = "results/analysis/sys1_olaparib_pmf_c3.xvg",   ligand = "Olaparib",    system = "S1"),
  S1_rucaparib   = list(path = "results/analysis/sys1_rucaparib_pmf_c3.xvg",  ligand = "Rucaparib",   system = "S1"),
  S1_talazoparib = list(path = "results/analysis/sys1_talazoparib_pmf_c3.xvg",ligand = "Talazoparib", system = "S1"),
  S1_veliparib   = list(path = "results/analysis/sys1_veliparib_pmf_c3.xvg",  ligand = "Veliparib",   system = "S1"),
  # S2 (DNA-bound): HD-ART PMF
  S2_APO         = list(path = "results/analysis/pmf-c3-sys2_APO_CV2_cv.dat.xvg",         ligand = "APO",        system = "S2"),
  S2_AZD5305     = list(path = "results/analysis/pmf-c3-sys2_AZD5305_CV2_cv.dat.xvg",     ligand = "AZD5305",    system = "S2"),
  S2_niraparib   = list(path = "results/analysis/pmf-c3-sys2_niraparib_CV2_cv.dat.xvg",   ligand = "Niraparib",  system = "S2"),
  S2_olaparib    = list(path = "results/analysis/pmf-c3-sys2_olaparib_CV2_cv.dat.xvg",    ligand = "Olaparib",   system = "S2"),
  S2_rucaparib   = list(path = "results/analysis/pmf-c3-sys2_rucaparib_CV2_cv.dat.xvg",   ligand = "Rucaparib",  system = "S2"),
  S2_talazoparib = list(path = "results/analysis/pmf-c3-sys2_talazoparib_CV2_cv.dat.xvg", ligand = "Talazoparib",system = "S2"),
  S2_veliparib   = list(path = "results/analysis/pmf-c3-sys2_veliparib_CV2_cv.dat.xvg",   ligand = "Veliparib",  system = "S2")
)

features <- lapply(systems, function(s) {
  df <- read_xvg(s$path)
  feat <- extract_features(df)
  feat$ligand <- s$ligand
  feat$system <- s$system
  feat$n_points <- nrow(df)
  feat
})
features_df <- bind_rows(features)

# ---- 3. Feature table ------------------------------------------------------
print("=== PMF Feature Summary ===")
print(features_df, n = 50)
write.csv(features_df, "results/analysis/pmf_features_summary.csv", row.names = FALSE)

# ---- 4. S1 vs S2 feature scatter -------------------------------------------
p_well_depth <- ggplot(features_df, aes(x = ligand, y = well_depth, fill = system)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("S1" = "#2166AC", "S2" = "#B2182B")) +
  labs(x = NULL, y = "Well Depth (kcal/mol)", title = "HD-ART Free Energy Well Depth") +
  theme_7pt + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_rc_min <- ggplot(features_df, aes(x = ligand, y = rc_min, fill = system)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("S1" = "#2166AC", "S2" = "#B2182B")) +
  labs(x = NULL, y = "RC at Minimum (Å)", title = "HD-ART Equilibrium Distance") +
  theme_7pt + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_rc_range <- ggplot(features_df, aes(x = ligand, y = rc_range, fill = system)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("S1" = "#2166AC", "S2" = "#B2182B")) +
  labs(x = NULL, y = "RC Range (Å)", title = "HD-ART Sampling Range") +
  theme_7pt + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_features <- (p_well_depth | p_rc_min | p_rc_range) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

cairo_pdf("results/figures/Fig_Feature_Comparison.pdf", width = 7.2, height = 3.5, pointsize = 7)
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

write.csv(ratio_df, "results/analysis/S1_S2_ratios.csv", row.names = FALSE)

# ---- 8. Summary figure: S1 well_depth vs S2 well_depth ---
p_ratio <- ggplot(ratio_df, aes(x = well_depth_S2, y = well_depth_S1, label = ligand)) +
  geom_point(aes(color = ligand), size = 2.5) +
  ggrepel::geom_text_repel(size = 2.5, force = 2, box.padding = 0.35, max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("APO" = "grey40", "AZD5305" = "darkorange",
                                 "Olaparib" = "#E41A1C", "Talazoparib" = "#FF7F00",
                                 "Veliparib" = "#377EB8",
                                 "Niraparib" = "#4DAF4A", "Rucaparib" = "#984EA3")) +
  annotate("text", x = 29, y = 85, label = "DNA-free amplifies\ntrapping differences",
           size = 2.5, color = "grey40", hjust = 0) +
  labs(x = "S2 DNA-bound Well Depth (kcal/mol)", 
       y = "S1 CAT-only Well Depth (kcal/mol)") +
  theme_7pt + theme(legend.position = "none")

cairo_pdf("results/figures/Fig_S1S2_Ratio.pdf", width = 4, height = 3.5, pointsize = 7)
print(p_ratio)
dev.off()

message("===== Descriptive analysis complete =====")
