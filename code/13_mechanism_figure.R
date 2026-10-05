#!/usr/bin/env Rscript
# 13-mechanism_figure.R — Core mechanism figure for PARP1 trapping paper
# Panels:
#   A: S1 PMF overlay (Type II vs III, all 7 systems)
#   B: S1 vs S2 well_depth bar chart (DNA-free amplification)
#   C: 2D mechanism scatter — S1/S2 ratio × S1 well_depth × trapping
#   D: Comparison: talazoparib vs veliparib extreme S1 PMF + structural inset
#
# Purpose:  Compose the core mechanism figure (trapping vs allostery master panel).
# Author:   Tao Zhu (tjogzt@gmail.com)
# Created:  2026-09-15 (header standardised 2026-10-05)
# Inputs:   data/cumulant_convergence.csv
# Depends:  data.table, dplyr, ggplot2, patchwork, tibble, tidyr
# Run:      Rscript code/13_mechanism_figure.R   (from the repository root)
library(ggplot2)
library(data.table)
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


# ---- Config ----------------------------------------------------------------
data_dir  <- "data/analysis"
out_dir   <- "figures/pdf"
stats_dir <- "results/analysis"   # CSV statistics tables (regenerable)
dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Theme and shared helpers (single source: common/helpers.R) -------------
# ---- System metadata (10 systems; single source: common/ligands.csv) --------
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)[1]])))
source(file.path(script_dir, "..", "common", "ligands.R"))
source(file.path(script_dir, "..", "common", "helpers.R"))
sys_meta <- load_ligands()
sys_meta$type <- gsub("_", " ", sys_meta$class)  # "Type II" / "Type III"

# ---- Read PMF data ----------------------------------------------------------


# S1 C3 PMF
s1_pmf <- list()
for (i in seq_len(nrow(sys_meta))) {
  lig <- sys_meta$ligand[i]
  f <- file.path(data_dir, sprintf("sys1_%s_pmf_c3.xvg", lig))
  if (!file.exists(f)) next
  df <- read_pmf(f)
  df$ligand <- lig
  df$label  <- sys_meta$label[i]
  df$type   <- sys_meta$type[i]
  df$color  <- sys_meta$color[i]
  s1_pmf[[lig]] <- df
}
s1_all <- bind_rows(s1_pmf)

# S2 CV2 (HD-ART) C3 PMF
s2_pmf <- list()
for (i in seq_len(nrow(sys_meta))) {
  lig <- sys_meta$ligand[i]
  f <- file.path(data_dir, sprintf("pmf-c3-sys2_%s_CV2_cv.dat.xvg", lig))
  if (!file.exists(f)) next
  df <- read_pmf(f)
  df$ligand <- lig
  df$label  <- sys_meta$label[i]
  df$type   <- sys_meta$type[i]
  df$color  <- sys_meta$color[i]
  s2_pmf[[lig]] <- df
}
s2_all <- bind_rows(s2_pmf)

# ---- Extract well_depth and RC_min per system --------------------------------
extract_stats <- function(pmf_list, sys_label) {
  stats <- lapply(pmf_list, function(df) {
    data.frame(
      ligand    = unique(df$ligand),
      label     = unique(df$label),
      type      = unique(df$type),
      color     = unique(df$color),
      well_depth = max(df$PMF_norm, na.rm = TRUE),
      rc_min    = df$RC[which.min(df$PMF)],
      rc_range  = max(df$RC) - min(df$RC),
      n_states  = {
        d1 <- diff(df$PMF)
        sum(diff(sign(d1)) == 2)
      },
      system    = sys_label,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(stats)
}

s1_stats <- extract_stats(s1_pmf, "S1")
s2_stats <- extract_stats(s2_pmf, "S2")

# Merge S1 + S2
combined <- inner_join(
  s1_stats %>% rename(wd_S1 = well_depth, rcmin_S1 = rc_min, nstates_S1 = n_states, rcrange_S1 = rc_range),
  s2_stats %>% rename(wd_S2 = well_depth, rcmin_S2 = rc_min, nstates_S2 = n_states, rcrange_S2 = rc_range),
  by = c("ligand", "label", "type", "color")
) %>% mutate(
  wd_ratio  = wd_S1 / wd_S2,
  delta_wd  = wd_S1 - wd_S2
)
combined <- left_join(combined, sys_meta[, c("ligand", "trap_potency", "shape")], by = "ligand")

# ---- AAI SD (manuscript Table 1 "AAI ±" column) ----------------------------
# Canonical propagation from the C1-C3 cumulant spread: aai_sd = wd_sd(C1-C3) / wd_S2
# (same basis as the extension-panel script s2_new_drugs_aai.py:48).
# Hard asserts pin the manuscript Table 1 values (2026-09-16 recomputation; see data_manifest.md).
cc_spread <- fread("data/cumulant_convergence.csv")
cc_sd <- setNames(cc_spread$wd_sd, cc_spread$ligand)
anch <- read.csv("data/01_curated/manuscript_anchor_values.csv", stringsAsFactors = FALSE)
av <- setNames(anch$value, anch$key)
ms_aai_sd <- av[paste0("aai_sd_", c("APO", "AZD5305", "olaparib", "talazoparib",
                                    "veliparib", "niraparib", "rucaparib"))]
names(ms_aai_sd) <- c("APO", "AZD5305", "olaparib", "talazoparib", "veliparib", "niraparib", "rucaparib")
aai_sd_comp <- cc_sd[names(ms_aai_sd)] /
               combined$wd_S2[match(names(ms_aai_sd), combined$ligand)]
stopifnot(all(abs(round(aai_sd_comp, 2) - ms_aai_sd) < 0.005))
cat("AAI SD asserts passed (Table 1 column pinned).\n")

cat("\n=== Combined S1/S2 Stats with talazoparib ===\n")
print(as.data.frame(combined %>% dplyr::arrange(desc(wd_ratio))))

# Shared Fig-2 panel theme: no top/right borders, left/bottom axes only
theme_f2 <- theme_7pt + theme(panel.border = element_blank(),
                              axis.line = element_line(linewidth = 0.5, color = "black"))

# ---- Panel A: S1 PMF overlay (all systems) ----------------------------------
p_a <- ggplot(s1_all, aes(x = RC, y = PMF_norm, color = label)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label)) +
  labs(x = "HD-ART Distance (Å)", y = "Free Energy (kcal/mol)",
       title = "A  S1: CAT-only HD-ART Free Energy Landscape") +
  theme_f2 +
  theme(legend.position = "bottom", legend.title = element_blank())

# ---- Panel B: S1 vs S2 well depth bar ---------------------------------------
bar_data <- combined %>%
  select(label, wd_S1, wd_S2) %>%
  pivot_longer(cols = c(wd_S1, wd_S2), names_to = "system", values_to = "well_depth") %>%
  mutate(system = factor(system, levels = c("wd_S1", "wd_S2"), 
                         labels = c("S1: CAT-only", "S2: DNA-bound")))

# Panel B shows the seven core systems (APO, AZD5305, five clinical inhibitors);
# the three extension inhibitors remain on display in panels A and C.
bar_data <- bar_data[bar_data$label %in% c("APO", "AZD5305", "Talazoparib", "Olaparib", "Niraparib", "Rucaparib", "Veliparib"), ]
bar_data$label <- factor(bar_data$label, 
                         levels = c("APO", "AZD5305", "Talazoparib", "Olaparib", "Niraparib", "Rucaparib", "Veliparib"))

label_colors <- setNames(
  c("grey40", "darkorange", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#377EB8"),
  c("APO", "AZD5305", "Talazoparib", "Olaparib", "Niraparib", "Rucaparib", "Veliparib")
)

p_b <- ggplot(bar_data, aes(x = label, y = well_depth, fill = system)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("S1: CAT-only" = "#2166AC", "S2: DNA-bound" = "#B2182B")) +
  labs(x = NULL, y = "Well Depth (kcal/mol)",
       title = "B  HD-ART Energy Landscape:\nDNA-free vs DNA-bound") +
  theme_f2 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom", legend.title = element_blank())

# ---- Panel C: 2D mechanism scatter -------------------------------------------
# x = S1/S2 ratio (allosteric amplification), y = S1 well depth
# Shape encodes the HX-MS Type classification (Zandarashvili 2020): the five clinical
# inhibitors are typed (II / III); AZD5305, the extension inhibitors, and APO are unclassified.
combined$type_shape <- ifelse(combined$label %in% c("Olaparib", "Talazoparib"), "Type II",
                       ifelse(combined$label %in% c("Veliparib", "Niraparib", "Rucaparib"), "Type III",
                              "Unclassified / control"))

p_c <- ggplot(combined, aes(x = wd_ratio, y = wd_S1, color = label, shape = type_shape)) +
  geom_point(size = 3, stroke = 1) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.9, show.legend = FALSE,
                           seed = 49,
                           max.overlaps = Inf, min.segment.length = 0.2,
                           box.padding = 0.3, force = 2, family = "Arial") +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label), guide = "none") +
  scale_shape_manual(values = c("Type II" = 16, "Type III" = 17, "Unclassified / control" = 15),
                     name = NULL,
                     guide = guide_legend(order = 2,
                                          override.aes = list(color = "grey30", size = 2.6))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = "S1/S2 Well Depth Ratio (Allosteric Amplification)",
       y = "S1 Well Depth (kcal/mol)",
       title = "C  Two-Dimensional Allosteric Mechanism Map") +
  theme_f2 +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.box = "vertical") +
  xlim(0.35, 1.7)

# ---- Panel D: Extreme comparison (talazoparib vs veliparib) -----------------
extreme_df <- s1_all %>% filter(label %in% c("Talazoparib", "Veliparib", "APO"))

# Add shaded Type regions
p_d <- ggplot(extreme_df, aes(x = RC, y = PMF_norm, color = label)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("APO" = "grey40", "Talazoparib" = "#FF7F00", "Veliparib" = "#377EB8")) +
  labs(x = "HD-ART Distance (Å)", y = "Free Energy (kcal/mol)",
       title = "D  Same CAT Pocket,\nOpposite Allosteric Fate") +
  theme_f2 +
  theme(legend.position = "bottom", legend.title = element_blank())

# ---- Assemble 4-panel master figure -----------------------------------------
# 2×2 layout
master <- (p_a | p_b) / (p_c | p_d)

cairo_pdf(file.path(out_dir, "Fig_Mechanism_Master.pdf"), width = 157/25.4, height = 157/25.4, pointsize = 8)
print(master)
dev.off()

cat("\nSaved: Fig_Mechanism_Master.pdf\n")

# ---- Supplementary: standalone 2D scatter for publications -------------------
# Niraparib two-state overlay: neutral (original parameterization) + protonated.
plot_df2 <- combined %>% filter(type != "APO") %>%
  mutate(state = "base", shape2 = type)
nira_neut <- data.frame(ligand = "niraparib", label = "Niraparib\n(neutral)",
                        type = "Type III", color = "#4DAF4A",
                        wd_ratio = 49.4 / 53.4, wd_S1 = 49.4,
                        state = "nira_neutral", shape2 = "nira")
nira_prot <- combined %>% filter(ligand == "niraparib") %>%
  mutate(label = "Niraparib\n(protonated)", state = "nira_protonated",
         shape2 = "nira")
plot_df2 <- bind_rows(
  plot_df2 %>% filter(ligand != "niraparib"),
  nira_neut, nira_prot
)
p_2d <- ggplot(plot_df2, aes(x = wd_ratio, y = wd_S1)) +
  geom_point(aes(color = type, shape = shape2), size = 4) +
  scale_shape_manual(values = c("Type II" = 16, "Type III" = 17,
                                "Extension" = 18, "Unknown" = 15,
                                "nira" = 1),
                     breaks = c("Type II", "Type III", "Extension", "Unknown"),
                     name = "Classification") +
  ggrepel::geom_text_repel(aes(label = label), size = 3.0, max.overlaps = Inf,
                           seed = 49,
                           min.segment.length = 0.3, box.padding = 0.6,
                           force = 4, max.iter = 5000, family = "Arial") +
  scale_color_manual(values = c("Type II" = "#E41A1C", "Type III" = "#377EB8",
                                "Extension" = "#C23531", "Unknown" = "darkorange"),
                     breaks = c("Type II", "Type III", "Extension", "Unknown"),
                     name = "Classification") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = "S1/S2 Well Depth Ratio (Allosteric Amplification Index)",
       y = "S1 HD-ART Well Depth (kcal/mol)") +
  theme_7pt +
  theme(legend.position = "right",
        plot.margin = ggplot2::margin(1, 4, 1, 3, unit = "mm"))

cairo_pdf(file.path(out_dir, "Fig_2D_Mechanism_Map.pdf"), width = 120/25.4, height = 72.5/25.4, pointsize = 8)
print(p_2d)
dev.off()

cat("Saved: Fig_2D_Mechanism_Map.pdf\n")

# ---- Save combined stats ----------------------------------------------------
write.csv(combined, file.path(stats_dir, "Fig_Mechanism_Data.csv"), row.names = FALSE)
message("===== Mechanism figure complete =====")
