# 09-s1_pmf_figures.R — S1 PMF comparison figures
# HD–ART domain distance PMF curves for 7 systems (APO + 6 inhibitors)
#
# Purpose:  Render the S1 HD–ART PMF comparison figures for the seven systems (overlay + facet panels).
# Created:  2026-09-15 (header standardised 2026-10-05)
# Depends:  dplyr, ggplot2, patchwork
# Run:      Rscript code/09_s1_pmf_figures.R   (from the repository root)
library(ggplot2)
library(dplyr)
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

# ---- System metadata (single source: common/ligands.csv) --------------------
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)[1]])))
source(file.path(script_dir, "..", "common", "ligands.R"))
sys_meta <- load_ligands()
# Display labels for this figure: plain label + class annotation
label_suffix <- c(APO = "(no ligand)", Unknown = "(blind)",
                  Type_II = "(Type II)", Type_III = "(Type III)", Extension = "(ext.)")
sys_meta$label <- paste0(sys_meta$label, "\n", label_suffix[sys_meta$class])

# ---- Load PMF data ---------------------------------------------------------


all_pmf <- list()
for (i in seq_len(nrow(sys_meta))) {
  lig <- sys_meta$ligand[i]
  f <- file.path(data_dir, sprintf("sys1_%s_pmf_c3.xvg", lig))
  if (!file.exists(f)) {
    message("Missing: ", f)
    next
  }
  df <- read_pmf(f)
  df$ligand  <- lig
  df$label   <- sys_meta$label[i]
  df$class   <- sys_meta$class[i]
  df$color   <- sys_meta$color[i]
  all_pmf[[lig]] <- df
}
pmf_all <- bind_rows(all_pmf)
pmf_all$ligand <- factor(pmf_all$ligand, levels = sys_meta$ligand)

# ---- A. Multi-panel: one PMF per system (3×3 grid) -------------------------
theme_pmf <- theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    plot.margin       = ggplot2::margin(1, 1, 1, 3, unit = "mm"),
    panel.grid.minor = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(size = 9, face = "bold"),
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 8)
  )

p_list <- lapply(sys_meta$ligand, function(lig) {
  df <- filter(pmf_all, ligand == lig)
  meta <- sys_meta[sys_meta$ligand == lig, ]
  pmf_min <- df$RC[which.min(df$PMF)]
  
  ggplot(df, aes(x = RC, y = PMF_norm)) +
    geom_line(color = meta$color, linewidth = 0.4) +
    geom_vline(xintercept = pmf_min, color = meta$color, 
               linetype = "dashed", linewidth = 0.3) +
    labs(title = meta$label, x = "HD–ART Distance (Å)", y = "PMF (kcal/mol)") +
    annotate("text", x = pmf_min, y = max(df$PMF_norm) * 0.85,
             label = sprintf("%.1f Å", pmf_min), 
             hjust = -0.15, size = 2.9, color = meta$color) +
    theme_pmf
})
names(p_list) <- sys_meta$ligand

# Arrange in 2×3
wrap_order <- c("APO", "olaparib", "talazoparib",
                "veliparib", "niraparib", "rucaparib", "AZD5305")
panel_a <- wrap_plots(p_list[wrap_order], ncol = 3, nrow = 3)

ggsave(file.path(out_dir, "Fig_S1_pmf_panels.pdf"), panel_a,
       width = 157/25.4, height = 104.7/25.4, device = cairo_pdf)
message("Saved: Fig_S1_pmf_panels.pdf")

# ---- B. Overlay: all systems c3 PMF ----------------------------------------
overlay_known <- filter(pmf_all, class %in% c("Type_II", "Type_III"))
overlay_other <- filter(pmf_all, !class %in% c("Type_II", "Type_III"))
# Short legend labels for the overlay (class is color-coded; long "(Type II)"
# suffixes overflow the 4.5-in figure width)
overlay_known$label <- sub(" \\((Type (II|III)|ext\\.)\\)$", "", overlay_known$label)
overlay_other$label <- sub(" \\((Type (II|III)|ext\\.)\\)$", "", overlay_other$label)

panel_b <- ggplot() +
  # APO + AZD5305 as reference
  geom_line(data = overlay_other, 
            aes(x = RC, y = PMF_norm, color = label, linetype = label),
            linewidth = 0.5) +
  # Known inhibitors
  geom_line(data = overlay_known,
            aes(x = RC, y = PMF_norm, color = label),
            linewidth = 0.6) +
  scale_color_manual(values = setNames(sys_meta$color, sys_meta$label)) +
  scale_linetype_manual(values = c("APO (no ligand)" = "dotted", 
                                   "AZD5305 (blind)" = "dashed")) +
  labs(x = "HD–ART Distance (Å)", y = "PMF (kcal/mol)",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    plot.margin       = ggplot2::margin(1, 1, 1, 3, unit = "mm"),
    legend.position   = "bottom",
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.35, "cm"),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 8)
  ) +
  theme_open +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(out_dir, "Fig_S1_pmf_overlay.pdf"), panel_b,
       width = 140/25.4, height = 108.9/25.4, device = cairo_pdf)
message("Saved: Fig_S1_pmf_overlay.pdf")

# ---- C. C1/C2/C3 comparison per inhibitor (4 panels) -----------------------
# Read all cumulant orders
read_all_cumulants <- function(ligand) {
  orders <- c("c1", "c2", "c3")
  dfs <- lapply(orders, function(o) {
    f <- file.path(data_dir, sprintf("sys1_%s_pmf_%s.xvg", ligand, o))
    if (!file.exists(f)) return(NULL)
    df <- read_pmf(f)
    df$cumulant <- o
    df
  })
  bind_rows(dfs)
}

inhibitors <- c("olaparib", "veliparib", "niraparib", "rucaparib")
p_cum <- lapply(inhibitors, function(lig) {
  df <- read_all_cumulants(lig)
  meta <- sys_meta[sys_meta$ligand == lig, ]
  
  ggplot(df, aes(x = RC, y = PMF_norm, color = cumulant)) +
    geom_line(linewidth = 0.4) +
    scale_color_brewer(palette = "Set1", 
                       labels = c("c1"="1st", "c2"="2nd", "c3"="3rd")) +
    labs(title = meta$label, x = "HD–ART Distance (Å)", 
         y = "PMF (kcal/mol)", color = "Cumulant") +
    theme_pmf +
    theme(legend.position = c(0.78, 0.72),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.35, "cm"))
})

panel_c <- wrap_plots(p_cum, ncol = 4, nrow = 1)

ggsave(file.path(out_dir, "Fig_S1_pmf_cumulants.pdf"), panel_c,
       width = 149/25.4, height = 46.4/25.4, device = cairo_pdf)
message("Saved: Fig_S1_pmf_cumulants.pdf")

# ---- Summary statistics ----------------------------------------------------
pmf_stats <- pmf_all %>%
  group_by(ligand, label, class) %>%
  summarise(
    n_points    = n(),
    rc_min      = RC[which.min(PMF)],
    pmf_range   = max(PMF_norm, na.rm = TRUE),
    rc_range    = max(RC) - min(RC),
    .groups     = "drop"
  )
print(pmf_stats)
write.csv(pmf_stats, file.path(stats_dir, "S1_PMF_c3_stats.csv"), row.names = FALSE)

message("===== All figures saved to ", out_dir, " =====")
