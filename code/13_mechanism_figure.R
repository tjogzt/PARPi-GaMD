#!/usr/bin/env Rscript
# 13_mechanism_figure.R — Core mechanism figure for PARP1 trapping paper
# Panels:
#   A: S1 PMF overlay (Type II vs III, all 7 systems)
#   B: S1 vs S2 well_depth bar chart (DNA-free amplification)
#   C: 2D mechanism scatter — S1/S2 ratio × S1 well_depth × trapping
#   D: Comparison: talazoparib vs veliparib extreme S1 PMF + structural inset
#
# Purpose:  Compose the core mechanism figure pair: Mechanism_Master (panels A-D) and
#           Fig_2D_Mechanism_Map (S1/S2 well-depth vs allosteric amplification map).
# Created:  2026-09-15 (header standardised 2026-10-05; R6 audit extension 2026-10-10)
# Inputs:   data/analysis/sys*_pmf_c3.xvg ; data/analysis/r6/r6_panel_values.csv ;
#           data/01_curated/two_state_panel_values.csv
# Outputs:  figures/pdf/Fig_Mechanism_Master.pdf ; figures/pdf/Fig_2D_Mechanism_Map.pdf ;
#           data/analysis/Fig_Mechanism_Data.csv
# Depends:  ggplot2, data.table, dplyr, patchwork, tibble, tidyr, ggrepel
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
      well_depth = NA_real_,  # filled from the audited occupied-support C3 table below
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

# ---- R6 audit: well depths = occupied-support C3 spans (Table S1/S5 convention) ----
# Replaces the raw min-normalized curve maxima (tail-dominated legacy convention).
pv <- fread("data/analysis/r6/r6_panel_values.csv")
key <- tolower(combined$ligand)
combined$wd_S1 <- pv$S1_C3[match(key, pv$ligand)]
combined$wd_S2 <- pv$S2CV2_C3[match(key, pv$ligand)]
combined$wd_ratio <- combined$wd_S1 / combined$wd_S2
combined$delta_wd  <- combined$wd_S1 - combined$wd_S2
stopifnot(all(is.finite(combined$wd_ratio)),
          all(abs(combined$wd_ratio - pv$AAI_C3[match(key, pv$ligand)]) < 1e-5))

# (legacy AAI-SD pin block removed in R6e: values now live in Table 1 via the R6 rebuild)

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
       title = "A  S1: CAT-only\nHD-ART free-energy landscape") +
  theme_f2 +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.margin = ggplot2::margin(0, 0, 1, 0, unit = "pt"),
        legend.box.spacing = grid::unit(6.5, "pt"),
        legend.key.width = grid::unit(0.7, "lines"),
        legend.spacing.x = grid::unit(0.1, "lines"),
        axis.title.x = element_text(margin = ggplot2::margin(t = 4.5)))

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
  labs(x = NULL, y = "span (kcal/mol)",
       title = "B  HD-ART Energy Landscape:\nDNA-free vs DNA-bound") +
  theme_f2 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom", legend.title = element_blank(),
        legend.margin = ggplot2::margin(0, 0, 1, 0, unit = "pt"),
        legend.box.spacing = grid::unit(2, "pt"))

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
  labs(x = "S1/S2 span ratio",
       y = "S1 sampled-range span (kcal/mol)",
       title = "C  Span ratio versus S1 span") +
  theme_f2 +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.box = "vertical") +
  xlim(0.4, 4.7)

# ---- Panel D: Extreme comparison (talazoparib vs veliparib) -----------------
extreme_df <- s1_all %>% filter(label %in% c("Talazoparib", "Veliparib", "APO"))

# Add shaded Type regions
p_d <- ggplot(extreme_df, aes(x = RC, y = PMF_norm, color = label)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("APO" = "grey40", "Talazoparib" = "#FF7F00", "Veliparib" = "#377EB8")) +
  labs(x = "HD-ART Distance (Å)", y = "Free Energy (kcal/mol)",
       title = "D  S1 profile contrast of trapping extremes") +
  theme_f2 +
  theme(legend.position = "bottom", legend.title = element_blank())

# ---- Assemble 4-panel master figure -----------------------------------------
# 2×2 layout
master <- (free(p_a) | p_b) / (p_c | p_d)  # free(): opt out of patchwork cross-panel section alignment

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
                        wd_ratio = pv$AAI_C3[pv$ligand == "niraparib"],
                        wd_S1 = pv$S1_C3[pv$ligand == "niraparib"],
                        state = "nira_neutral", shape2 = "nira")
nira_prot <- data.frame(ligand = "niraparib", label = "Niraparib\n(protonated)",
                        type = "Type III", color = "#4DAF4A",
                        wd_ratio = pv$AAI_C3[pv$ligand == "niraparib_prot"],
                        wd_S1 = pv$S1_C3[pv$ligand == "niraparib_prot"],
                        state = "nira_protonated", shape2 = "nira")
ruca_prot2 <- data.frame(ligand = "rucaparib", label = "Rucaparib\n(protonated)",
                         type = "Type III", color = "#984EA3",
                         wd_ratio = pv$AAI_C3[pv$ligand == "rucaparib_prot"],
                         wd_S1 = pv$S1_C3[pv$ligand == "rucaparib_prot"],
                         state = "ruca_protonated", shape2 = "ruca")
plot_df2 <- bind_rows(
  plot_df2 %>% filter(ligand != "niraparib"),
  nira_neut, nira_prot, ruca_prot2
)
p_2d <- ggplot(plot_df2, aes(x = wd_ratio, y = wd_S1)) +
  geom_point(aes(color = type, shape = shape2), size = 4) +
  scale_shape_manual(values = c("Type II" = 16, "Type III" = 17,
                                "Extension" = 18, "Unknown" = 15,
                                "nira" = 1, "ruca" = 1),
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
  labs(x = "S1/S2 span ratio",
       y = "S1 HD--ART span (kcal/mol)") +
  theme_7pt +
  theme(legend.position = "right",
        plot.margin = ggplot2::margin(1, 4, 1, 3, unit = "mm"))

cairo_pdf(file.path(out_dir, "Fig_2D_Mechanism_Map.pdf"), width = 120/25.4, height = 72.5/25.4, pointsize = 8)
print(p_2d)
dev.off()

cat("Saved: Fig_2D_Mechanism_Map.pdf\n")

# ---- Save combined stats ----------------------------------------------------
write.csv(combined, "data/analysis/Fig_Mechanism_Data.csv", row.names = FALSE)
message("===== Mechanism figure complete =====")
