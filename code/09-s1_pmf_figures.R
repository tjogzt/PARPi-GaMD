# 01-s1_pmf_figures.R — S1 PMF comparison figures
# HD-ART domain distance PMF curves for 6 systems
library(ggplot2)
library(dplyr)
library(patchwork)

# ---- Config ----------------------------------------------------------------
data_dir  <- "results/analysis"
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- System metadata (single source: common/ligands.csv) --------------------
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)[1]])))
source(file.path(script_dir, "..", "common", "ligands.R"))
sys_meta <- load_ligands()
# Display labels for this figure: plain label + class annotation
label_suffix <- c(APO = "(no ligand)", Unknown = "(blind)",
                  Type_II = "(Type II)", Type_III = "(Type III)", Extension = "(ext.)")
sys_meta$label <- paste(sys_meta$label, label_suffix[sys_meta$class])

# ---- Load PMF data ---------------------------------------------------------
read_xvg <- function(path) {
  lines <- readLines(path)
  data_start <- which(grepl("^[0-9]", lines))[1]
  if (is.na(data_start)) stop("No data found in ", path)
  df <- read.table(text = lines[data_start:length(lines)], 
                   header = FALSE, col.names = c("RC", "PMF"))
  # Normalize PMF to minimum = 0
  df$PMF_norm <- df$PMF - min(df$PMF, na.rm = TRUE)
  df
}

all_pmf <- list()
for (i in seq_len(nrow(sys_meta))) {
  lig <- sys_meta$ligand[i]
  f <- file.path(data_dir, sprintf("sys1_%s_pmf_c3.xvg", lig))
  if (!file.exists(f)) {
    message("Missing: ", f)
    next
  }
  df <- read_xvg(f)
  df$ligand  <- lig
  df$label   <- sys_meta$label[i]
  df$class   <- sys_meta$class[i]
  df$color   <- sys_meta$color[i]
  all_pmf[[lig]] <- df
}
pmf_all <- bind_rows(all_pmf)
pmf_all$ligand <- factor(pmf_all$ligand, levels = sys_meta$ligand)

# ---- A. Multi-panel: one PMF per system (2×3) ------------------------------
theme_pmf <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(size = 7, face = "bold"),
    axis.title        = element_text(size = 7),
    axis.text         = element_text(size = 6)
  )

p_list <- lapply(sys_meta$ligand, function(lig) {
  df <- filter(pmf_all, ligand == lig)
  meta <- sys_meta[sys_meta$ligand == lig, ]
  pmf_min <- df$RC[which.min(df$PMF)]
  
  ggplot(df, aes(x = RC, y = PMF_norm)) +
    geom_line(color = meta$color, linewidth = 0.4) +
    geom_vline(xintercept = pmf_min, color = meta$color, 
               linetype = "dashed", linewidth = 0.3) +
    labs(title = meta$label, x = "HD-ART Distance (Å)", y = "PMF (kcal/mol)") +
    annotate("text", x = pmf_min, y = max(df$PMF_norm) * 0.85,
             label = sprintf("%.1f Å", pmf_min), 
             hjust = -0.15, size = 2.2, color = meta$color) +
    theme_pmf
})
names(p_list) <- sys_meta$ligand

# Arrange in 2×3
wrap_order <- c("APO", "olaparib", "talazoparib",
                "veliparib", "niraparib", "rucaparib", "AZD5305")
panel_a <- wrap_plots(p_list[wrap_order], ncol = 3, nrow = 3) +
  plot_annotation(
    title = "S1: HD-ART PMF (3rd cumulant expansion)",
    theme = theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
  )

ggsave(file.path(out_dir, "Fig_S1_pmf_panels.pdf"), panel_a,
       width = 7.2, height = 4.8, device = cairo_pdf)
message("Saved: Fig_S1_pmf_panels.pdf")

# ---- B. Overlay: all systems c3 PMF ----------------------------------------
overlay_known <- filter(pmf_all, class %in% c("Type_II", "Type_III"))
overlay_other <- filter(pmf_all, !class %in% c("Type_II", "Type_III"))

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
  labs(x = "HD-ART Distance (Å)", y = "PMF (kcal/mol)",
       title = "S1: C3 PMF Overlay — Type II vs Type III PARP1 Inhibitors",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 7) +
  theme(
    legend.position   = "bottom",
    legend.text       = element_text(size = 6),
    legend.key.size   = unit(0.3, "cm"),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title        = element_text(size = 7),
    axis.text         = element_text(size = 6)
  )

ggsave(file.path(out_dir, "Fig_S1_pmf_overlay.pdf"), panel_b,
       width = 4.5, height = 3.5, device = cairo_pdf)
message("Saved: Fig_S1_pmf_overlay.pdf")

# ---- C. C1/C2/C3 comparison per inhibitor (4 panels) -----------------------
# Read all cumulant orders
read_all_cumulants <- function(ligand) {
  orders <- c("c1", "c2", "c3")
  dfs <- lapply(orders, function(o) {
    f <- file.path(data_dir, sprintf("sys1_%s_pmf_%s.xvg", ligand, o))
    if (!file.exists(f)) return(NULL)
    df <- read_xvg(f)
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
    labs(title = meta$label, x = "HD-ART Distance (Å)", 
         y = "PMF (kcal/mol)", color = "Cumulant") +
    theme_pmf +
    theme(legend.position = c(0.85, 0.7),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.key.size = unit(0.25, "cm"))
})

panel_c <- wrap_plots(p_cum, ncol = 4, nrow = 1) +
  plot_annotation(
    title = "S1: Cumulant expansion convergence (C1→C2→C3)",
    theme = theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
  )

ggsave(file.path(out_dir, "Fig_S1_pmf_cumulants.pdf"), panel_c,
       width = 9, height = 2.8, device = cairo_pdf)
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
write.csv(pmf_stats, file.path(out_dir, "S1_PMF_c3_stats.csv"), row.names = FALSE)

message("===== All figures saved to ", out_dir, " =====")
