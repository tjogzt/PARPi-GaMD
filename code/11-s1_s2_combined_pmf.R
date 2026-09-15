#!/usr/bin/env Rscript
# 11-s1_s2_combined_pmf.R — S1 vs S2 HD-ART PMF comparison (13 trajectories)
library(ggplot2)
library(dplyr)
library(patchwork)

# ---- Config ----------------------------------------------------------------
data_dir  <- "results/analysis"
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# System metadata with both S1 and S2
sys_meta <- data.frame(
  ligand   = c("APO", "AZD5305", "olaparib", "talazoparib", "veliparib", "niraparib", "rucaparib"),
  label    = c("APO", "AZD5305", "Olaparib", "Talazoparib", "Veliparib", "Niraparib", "Rucaparib"),
  class    = c("APO", "Unknown", "Type_II", "Type_II", "Type_III", "Type_III", "Type_III"),
  color    = c("grey40", "darkorange", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3"),
  stringsAsFactors = FALSE
)

# ---- Load PMF data ---------------------------------------------------------
read_xvg <- function(path) {
  lines <- readLines(path)
  data_start <- which(grepl("^[0-9]", lines))[1]
  if (is.na(data_start)) stop("No data in ", path)
  df <- read.table(text = lines[data_start:length(lines)], 
                   header = FALSE, col.names = c("RC", "PMF"))
  df$PMF_norm <- df$PMF - min(df$PMF, na.rm = TRUE)
  df
}

# S1: CAT-only HD-ART PMF
s1_pmf <- list()
for (i in seq_len(nrow(sys_meta))) {
  lig <- sys_meta$ligand[i]
  f <- file.path(data_dir, sprintf("sys1_%s_pmf_c3.xvg", lig))
  if (!file.exists(f)) next
  df <- read_xvg(f)
  df$ligand  <- lig
  df$label   <- sys_meta$label[i]
  df$class   <- sys_meta$class[i]
  df$color   <- sys_meta$color[i]
  df$system  <- "S1: CAT-only"
  s1_pmf[[lig]] <- df
}
s1_all <- bind_rows(s1_pmf)

# S2: DNA-bound HD-ART PMF
s2_meta <- sys_meta
s2_meta$ligand[s2_meta$ligand == "olaparib"] <- "olaparib"  # same
s2_pmf <- list()
for (i in seq_len(nrow(s2_meta))) {
  lig <- s2_meta$ligand[i]
  f <- file.path(data_dir, sprintf("pmf-c3-sys2_%s_CV2_cv.dat.xvg", lig))
  if (!file.exists(f)) next
  df <- read_xvg(f)
  df$ligand  <- lig
  df$label   <- s2_meta$label[i]
  df$class   <- s2_meta$class[i]
  df$color   <- s2_meta$color[i]
  df$system  <- "S2: DNA-bound"
  s2_pmf[[lig]] <- df
}
s2_all <- bind_rows(s2_pmf)

# Combine
pmf_all <- bind_rows(s1_all, s2_all)
pmf_all$system <- factor(pmf_all$system, levels = c("S1: CAT-only", "S2: DNA-bound"))

# ---- Theme -----------------------------------------------------------------
theme_pmf <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(size = 7, face = "bold"),
    axis.title        = element_text(size = 7),
    axis.text         = element_text(size = 6),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 7),
    legend.key.size   = unit(0.3, "cm"),
    strip.text        = element_text(size = 7, face = "bold"),
    strip.background  = element_rect(fill = "grey95")
  )

# ---- Fig A: Faceted S1 vs S2 overlay per ligand ----------------------------
p_facet <- ggplot(pmf_all, aes(x = RC, y = PMF_norm, color = system)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = c("S1: CAT-only" = "#2166AC", 
                                 "S2: DNA-bound" = "#B2182B")) +
  facet_wrap(~ label, ncol = 3, scales = "free_x") +
  labs(x = "HD-ART COM Distance (Å)", y = "Free Energy (kcal/mol)",
       title = "HD-ART Domain Motion PMF: CAT-only vs DNA-bound PARP1",
       color = NULL) +
  theme_pmf

# ---- Fig B: Overlay all S1 vs all S2 ---------------------------------------
p_overlay <- ggplot(pmf_all, aes(x = RC, y = PMF_norm, 
                                  color = label, linetype = system)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label)) +
  labs(x = "HD-ART COM Distance (Å)", y = "Free Energy (kcal/mol)",
       title = "S1 vs S2 HD-ART PMF Overlay (all systems)",
       color = NULL, linetype = "System") +
  theme_pmf +
  theme(legend.position = "bottom")

# ---- Fig C: Per-ligand S1+S2 panels ----------------------------------------
p_list <- lapply(sys_meta$ligand, function(lig) {
  df <- filter(pmf_all, ligand == lig)
  meta <- sys_meta[sys_meta$ligand == lig, ]
  
  ggplot(df, aes(x = RC, y = PMF_norm, color = system)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("S1: CAT-only" = "#2166AC", 
                                   "S2: DNA-bound" = "#B2182B")) +
    labs(title = meta$label, x = "HD-ART (Å)", y = "PMF") +
    theme_pmf +
    theme(legend.position = "none")
})

panel_c <- wrap_plots(p_list, ncol = 4, nrow = 2)

# ---- Save figures ----------------------------------------------------------
cairo_pdf(file.path(out_dir, "Fig_S1S2_hd_art_facet.pdf"), 
          width = 7.2, height = 4.8, pointsize = 7)
print(p_facet)
dev.off()

cairo_pdf(file.path(out_dir, "Fig_S1S2_hd_art_overlay.pdf"), 
          width = 5, height = 4, pointsize = 7)
print(p_overlay)
dev.off()

cairo_pdf(file.path(out_dir, "Fig_S1S2_hd_art_panels.pdf"), 
          width = 7.2, height = 4.8, pointsize = 7)
print(panel_c)
dev.off()

# Combined master figure
p_master <- wrap_plots(p_facet, p_overlay, ncol = 1, heights = c(1, 0.8)) +
  plot_annotation(title = "PARP1 HD-ART Free Energy Landscape: DNA-free vs DNA-bound",
                  theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)))
cairo_pdf(file.path(out_dir, "Fig_Master_S1S2_HD_ART.pdf"), 
          width = 7.2, height = 8, pointsize = 7)
print(p_master)
dev.off()

# ---- Summary stats ---------------------------------------------------------
stats <- pmf_all %>%
  group_by(system, ligand, label) %>%
  summarise(
    n_points    = n(),
    rc_min      = RC[which.min(PMF)],
    pmf_range   = max(PMF_norm, na.rm = TRUE),
    rc_range    = max(RC) - min(RC),
    .groups     = "drop"
  ) %>%
  arrange(ligand, system)
print(stats)
write.csv(stats, file.path(out_dir, "S1S2_HD_ART_stats.csv"), row.names = FALSE)

message("===== S1+S2 Combined figures saved to ", out_dir, " =====")
