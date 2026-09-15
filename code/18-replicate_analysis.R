#!/usr/bin/env Rscript
# 18-replicate_analysis.R — Process GaMD replicate data
# Reads PMF from replicate runs, computes cross-trajectory statistics,
# generates comparison figures, and updates Table 1.
#
# Expected input: runs/rep_{ligand}_{rep}/analysis/HD_ART_dist_pmf_c3.xvg
# (format matches analyze_gamd.py output)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

data_dir   <- "results/analysis"
rep_dir    <- "runs"
out_dir    <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Config ----
theme_7pt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))

# ---- PMF reader ----
read_pmf_xvg <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path)
  dstart <- which(grepl("^[0-9]", lines))[1]
  if (is.na(dstart)) return(NULL)
  df <- read.table(text = lines[dstart:length(lines)], 
                   col.names = c("RC", "PMF"))
  df$PMF_norm <- df$PMF - min(df$PMF, na.rm = TRUE)
  df
}

extract_well_depth <- function(pmf) {
  if (is.null(pmf)) return(c(NA, NA, NA))
  c(well_depth = max(pmf$PMF_norm, na.rm = TRUE),
    rc_min    = pmf$RC[which.min(pmf$PMF)],
    n_bins    = nrow(pmf))
}

# ---- Scan for replicate data ----
# AutoDL will produce data in: runs/rep_{ligand}_{rep}/analysis/
ligands <- c("talazoparib", "veliparib", "AZD5305", "EB47")
rep_ids <- 1:2

results <- data.frame()
for (lig in ligands) {
  for (rep in rep_ids) {
    tag <- sprintf("rep_%s_%d", ifelse(lig=="talazoparib","tala",
                                ifelse(lig=="veliparib","veli",
                                ifelse(lig=="AZD5305","azd","eb47")), rep)
    
    # Multiple possible locations for PMF file
    patterns <- c(
      sprintf("%s/%s/analysis/HD_ART_dist_pmf_c3.xvg", rep_dir, tag),
      sprintf("%s/%s/analysis/pmf_c3.xvg", rep_dir, tag),
      sprintf("%s/%s/analysis/*pmf*c3*.xvg", rep_dir, tag)
    )
    
    pmf <- NULL
    for (p in patterns) {
      if (grepl("\\*", p)) {
        files <- Sys.glob(p)
        if (length(files) > 0) pmf <- read_pmf_xvg(files[1])
      } else {
        pmf <- read_pmf_xvg(p)
      }
      if (!is.null(pmf)) break
    }
    
    wd <- extract_well_depth(pmf)
    results <- rbind(results, data.frame(
      ligand = lig,
      replicate = rep,
      tag = tag,
      well_depth = wd[1],
      rc_min = wd[2],
      n_bins = wd[3],
      stringsAsFactors = FALSE
    ))
  }
}

# ---- Also load original run data for comparison ----
original <- data.frame(
  ligand = c("talazoparib", "olaparib", "niraparib", "rucaparib", "veliparib", "AZD5305", "APO"),
  well_depth_orig = c(30.4, 68.9, 51.6, 53.4, 101.1, 28.2, 42.5),
  wd_sd_orig      = c(0.0, 14.9, 8.0, 7.2, 32.5, 0.7, 5.9),
  type = c("Type_II", "Type_II", "Type_III", "Type_III", "Type_III", "Unknown", "APO"),
  stringsAsFactors = FALSE
)

if (nrow(results) > 0) {
  # ---- Cross-replicate statistics ----
  cross_stats <- results %>%
    group_by(ligand) %>%
    summarise(
      n_reps = n(),
      wd_mean = mean(well_depth, na.rm = TRUE),
      wd_sd   = sd(well_depth, na.rm = TRUE),
      wd_min  = min(well_depth, na.rm = TRUE),
      wd_max  = max(well_depth, na.rm = TRUE),
      wd_range = wd_max - wd_min,
      rc_min_mean = mean(rc_min, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("=== Cross-Replicate Well Depth Statistics ===\n")
  print(as.data.frame(cross_stats), row.names = FALSE)
  
  # ---- Comparison with original ----
  merged <- merge(cross_stats, original, by = "ligand", all = TRUE)
  cat("\n=== Original vs Replicate Comparison ===\n")
  print(merged[, c("ligand", "type", "well_depth_orig", "wd_sd_orig", "wd_mean", "wd_sd")], 
        row.names = FALSE)
  
  # ---- Figure: Replicate comparison ----  
  # Panel A: Per-replicate well_depth scatter
  color_map <- c(talazoparib = "#FF7F00", veliparib = "#377EB8", 
                 AZD5305 = "darkorange", EB47 = "#E41A1C",
                 olaparib = "#E41A1C", niraparib = "#4DAF4A", rucaparib = "#984EA3")
  
  p1 <- ggplot(results, aes(x = ligand, y = well_depth, color = ligand, shape = factor(replicate))) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    scale_color_manual(values = color_map, guide = "none") +
    labs(x = NULL, y = "S1 Well Depth (kcal/mol)",
         title = "A  Replicate Consistency", shape = "Replicate") +
    theme_7pt + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Panel B: Cross-replicate SD vs original cumulant SD
  p2 <- ggplot(merged, aes(x = wd_sd_orig, y = wd_sd, label = ligand, color = type)) +
    geom_point(size = 3) +
    geom_text(hjust = -0.15, vjust = 0.5, size = 2.5, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(Type_II = "#E41A1C", Type_III = "#377EB8", Unknown = "darkorange")) +
    labs(x = "SD across C1-C3 cumulant (original)", y = "SD across replicates",
         title = "B  Cross-Trajectory vs Cumulant Uncertainty") +
    theme_7pt
  
  # Panel C: Combined table-style bar plot
  combined <- results %>%
    mutate(label = paste0(ligand, "_rep", replicate)) %>%
    bind_rows(data.frame(
      ligand = original$ligand,
      replicate = 0,
      tag = paste0(original$ligand, "_orig"),
      well_depth = original$well_depth_orig,
      rc_min = NA, n_bins = NA,
      label = paste0(original$ligand, "_orig"),
      stringsAsFactors = FALSE
    ))
  
  p3 <- ggplot(combined, aes(x = ligand, y = well_depth, fill = factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("0" = "grey50", "1" = "#2166AC", "2" = "#B2182B"),
                      labels = c("0" = "Original", "1" = "Rep 1", "2" = "Rep 2"),
                      name = NULL) +
    labs(x = NULL, y = "S1 Well Depth (kcal/mol)",
         title = "C  Original vs Replicate Well Depths") +
    theme_7pt + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p <- (p1 | p2) / p3 +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(title = "GaMD Replicate Validation",
                    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)))
  
  cairo_pdf(file.path(out_dir, "Fig_Replicate_Validation.pdf"), 
            width = 190/25.4, height = 170/25.4, pointsize = 7)
  print(p)
  dev.off()
  
  write.csv(merged, "results/analysis/replicate_cross_stats.csv", row.names = FALSE)
  write.csv(results, "results/analysis/replicate_raw_data.csv", row.names = FALSE)
  
  cat("\nSaved: Fig_Replicate_Validation.pdf\n")
  cat("Saved: replicate_cross_stats.csv, replicate_raw_data.csv\n")
  
} else {
  cat("No replicate data found yet. This script will work when AutoDL completes.\n")
  cat("Expected paths: runs/rep_{tala,veli,azd,eb47}_{1,2}/analysis/*pmf*c3*.xvg\n")
}

message("===== Replicate analysis pipeline ready =====")
