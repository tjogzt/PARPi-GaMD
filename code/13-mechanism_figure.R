#!/usr/bin/env Rscript
# 13-mechanism_figure.R — Core mechanism figure for PARP1 trapping paper
# Panels:
#   A: S1 PMF overlay (Type II vs III, all 7 systems)
#   B: S1 vs S2 well_depth bar chart (DNA-free amplification)
#   C: 2D mechanism scatter — S1/S2 ratio × S1 well_depth × trapping
#   D: Comparison: talazoparib vs veliparib extreme S1 PMF + structural inset
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)

# ---- Config ----------------------------------------------------------------
data_dir  <- "results/analysis"
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Theme (SCI-standard: Arial 7pt, cairo_pdf) ----------------------------
theme_7pt <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.2),
    plot.title        = element_text(size = 7, face = "bold"),
    axis.title        = element_text(size = 7),
    axis.text         = element_text(size = 6),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 7),
    legend.key.size   = unit(0.3, "cm"),
    strip.text        = element_text(size = 7, face = "bold"),
    strip.background  = element_rect(fill = "grey95")
  )

# ---- System metadata (7 systems) --------------------------------------------
sys_meta <- data.frame(
  ligand       = c("APO", "AZD5305", "olaparib", "talazoparib", "veliparib", "niraparib", "rucaparib",
                   "fluzoparib", "pamiparib", "senaparib"),
  label        = c("APO", "AZD5305", "Olaparib", "Talazoparib", "Veliparib", "Niraparib", "Rucaparib",
                   "Fluzoparib", "Pamiparib", "Senaparib"),
  type         = c("APO", "Unknown", "Type II", "Type II", "Type III", "Type III", "Type III",
                   "Extension", "Extension", "Extension"),
  trap_potency = c(NA, NA, 1.0, 100, 0.02, 65, 0.8, NA, NA, NA),  # relative to olaparib
  color        = c("grey40", "darkorange", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3",
                   "#C23531", "#3D6BA8", "#9D2933"),
  shape        = c(17, 15, 16, 16, 17, 17, 17, 18, 15, 8),  # extension: diamond/square/star
  stringsAsFactors = FALSE
)

# ---- Read PMF data ----------------------------------------------------------
read_xvg <- function(path) {
  lines <- readLines(path)
  dstart <- which(grepl("^[0-9]", lines))[1]
  df <- read.table(text = lines[dstart:length(lines)], col.names = c("RC", "PMF"))
  df$PMF_norm <- df$PMF - min(df$PMF, na.rm = TRUE)
  df
}

# S1 C3 PMF
s1_pmf <- list()
for (i in seq_len(nrow(sys_meta))) {
  lig <- sys_meta$ligand[i]
  f <- file.path(data_dir, sprintf("sys1_%s_pmf_c3.xvg", lig))
  if (!file.exists(f)) next
  df <- read_xvg(f)
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
  df <- read_xvg(f)
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

cat("\n=== Combined S1/S2 Stats with talazoparib ===\n")
print(as.data.frame(combined %>% dplyr::arrange(desc(wd_ratio))))

# ---- Panel A: S1 PMF overlay (all systems) ----------------------------------
p_a <- ggplot(s1_all, aes(x = RC, y = PMF_norm, color = label)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label)) +
  labs(x = "HD-ART Distance (Å)", y = "Free Energy (kcal/mol)",
       title = "A  S1: CAT-only HD-ART Free Energy Landscape") +
  theme_7pt +
  theme(legend.position = "bottom", legend.title = element_blank())

# ---- Panel B: S1 vs S2 well depth bar ---------------------------------------
bar_data <- combined %>%
  select(label, wd_S1, wd_S2) %>%
  pivot_longer(cols = c(wd_S1, wd_S2), names_to = "system", values_to = "well_depth") %>%
  mutate(system = factor(system, levels = c("wd_S1", "wd_S2"), 
                         labels = c("S1: CAT-only", "S2: DNA-bound")))

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
       title = "B  HD-ART Energy Landscape: DNA-free vs DNA-bound") +
  theme_7pt +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        legend.position = c(0.85, 0.85))

# ---- Panel C: 2D mechanism scatter -------------------------------------------
# x = S1/S2 ratio (allosteric amplification), y = S1 well depth
p_c <- ggplot(combined, aes(x = wd_ratio, y = wd_S1, color = label, shape = factor(shape))) +
  geom_point(size = 3, stroke = 1) +
  geom_text(aes(label = label), hjust = -0.15, vjust = 0.5, size = 2.2, show.legend = FALSE) +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label)) +
  scale_shape_manual(values = c("16" = 16, "17" = 17, "15" = 15), guide = "none") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  annotate("rect", xmin = 0.8, xmax = 1.2, ymin = 20, ymax = 35,
           fill = "grey90", alpha = 0.3) +
  annotate("text", x = 1, y = 38, label = "Neutral\nallostery", size = 2.2, color = "grey40") +
  annotate("text", x = 2.5, y = 95, label = "Pro-release\nallostery →", size = 2.2, color = "grey40") +
  labs(x = "S1/S2 Well Depth Ratio (Allosteric Amplification)",
       y = "S1 Well Depth (kcal/mol)",
       title = "C  Two-Dimensional Allosteric Mechanism Map") +
  theme_7pt +
  theme(legend.position = "none") +
  xlim(0.7, 4.0)

# ---- Panel D: Extreme comparison (talazoparib vs veliparib) -----------------
extreme_df <- s1_all %>% filter(label %in% c("Talazoparib", "Veliparib", "APO"))

# Add shaded Type regions
p_d <- ggplot(extreme_df, aes(x = RC, y = PMF_norm, color = label)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("APO" = "grey40", "Talazoparib" = "#FF7F00", "Veliparib" = "#377EB8")) +
  # Shade the well depth for veliparib
  annotate("rect", xmin = 22, xmax = 25.5, ymin = 95, ymax = 105,
           fill = "#377EB8", alpha = 0.08) +
  annotate("text", x = 23.75, y = 103, label = "ΔG = 101 kcal/mol\n(pro-release strain)",
           size = 2.2, color = "#377EB8", fontface = "italic") +
  annotate("text", x = 23.7, y = 35, label = "ΔG = 30 kcal/mol\n(stable closure)",
           size = 2.2, color = "#FF7F00", fontface = "italic") +
  labs(x = "HD-ART Distance (Å)", y = "Free Energy (kcal/mol)",
       title = "D  Same CAT Pocket, Opposite Allosteric Fate") +
  theme_7pt +
  theme(legend.position = c(0.15, 0.85), legend.background = element_rect(fill = alpha("white", 0.7)))

# ---- Assemble 4-panel master figure -----------------------------------------
# 2×2 layout
master <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = "PARP1 Trapping Mechanism: GaMD Reveals Allosteric Encoding in HD Domain",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

cairo_pdf(file.path(out_dir, "Fig_Mechanism_Master.pdf"), width = 190/25.4, height = 190/25.4, pointsize = 7)
print(master)
dev.off()

cat("\nSaved: Fig_Mechanism_Master.pdf\n")

# ---- Supplementary: standalone 2D scatter for publications -------------------
p_2d <- ggplot(combined %>% filter(type != "APO"), aes(x = wd_ratio, y = wd_S1)) +
  # Classification zones
  annotate("rect", xmin = 0.8, xmax = 1.5, ymin = 25, ymax = 40,
           fill = "#FF7F00", alpha = 0.08) +
  annotate("text", x = 1.15, y = 42, label = "Neutral (Type II)", 
           size = 2.5, color = "#E41A1C", fontface = "italic") +
  annotate("rect", xmin = 1.5, xmax = 4.0, ymin = 45, ymax = 110,
           fill = "#377EB8", alpha = 0.08) +
  annotate("text", x = 2.5, y = 108, label = "Pro-release (Type III)",
           size = 2.5, color = "#377EB8", fontface = "italic") +
  geom_point(aes(color = type, shape = type), size = 4) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                           min.segment.length = 0.3, box.padding = 0.35) +
  scale_color_manual(values = c("Type II" = "#E41A1C", "Type III" = "#377EB8",
                                "Extension" = "#C23531", "Unknown" = "darkorange"),
                     name = "Classification") +
  scale_shape_manual(values = c("Type II" = 16, "Type III" = 17,
                                "Extension" = 18, "Unknown" = 15),
                     name = "Classification") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = "S1/S2 Well Depth Ratio (Allosteric Amplification Index)",
       y = "S1 HD-ART Well Depth (kcal/mol)",
       title = "Allosteric Mechanism Map of PARP1 Inhibitors") +
  theme_7pt +
  theme(legend.position = c(0.85, 0.15))

cairo_pdf(file.path(out_dir, "Fig_2D_Mechanism_Map.pdf"), width = 90/25.4, height = 75/25.4, pointsize = 7)
print(p_2d)
dev.off()

cat("Saved: Fig_2D_Mechanism_Map.pdf\n")

# ---- Save combined stats ----------------------------------------------------
write.csv(combined, file.path(out_dir, "Fig_Mechanism_Data.csv"), row.names = FALSE)
message("===== Mechanism figure complete =====")
