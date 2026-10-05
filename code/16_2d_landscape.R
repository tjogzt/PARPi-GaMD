#!/usr/bin/env Rscript
# 16-2d_landscape.R — 2D free energy landscapes: talazoparib vs veliparib
# Uses CV1 (Protein-DNA) × CV2 (HD-ART) from .npy files + reweighting weights
# Requires RETICULATE_PYTHON pointing at a numpy-enabled Python, e.g.:
#   RETICULATE_PYTHON=/opt/anaconda3/bin/python3 Rscript code/16-2d_landscape.R
#
# Purpose:  Render the 2D free-energy landscapes (CV1 x CV2): talazoparib-vs-veliparib main panel + all-system SI variant.
# Author:   Tao Zhu (tjogzt@gmail.com)
# Created:  2026-09-15 (header standardised 2026-10-05)
# Depends:  MASS, dplyr, ggplot2, patchwork, reticulate
# Run:      Rscript code/16_2d_landscape.R   (from the repository root)
library(ggplot2)
library(dplyr)
library(MASS)  # for kde2d
library(reticulate)
library(patchwork)

np <- import("numpy")

data_dir  <- "data/analysis"
out_dir   <- "figures/pdf"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Shared helpers: theme_7pt (single source).
args_h <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", args_h))) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
  source(file.path(script_dir, "..", "common", "helpers.R"))
} else {
  source("common/helpers.R")
}

# ---- Load precomputed 2D cumulant-expansion PMF grids (DBE reweighting) ----
# Grids are produced by scripts/compute_2d_dbe.py (PyReweighting-2D,
# amdweight_CE, textbook e^{beta dV_D} weights; consistent with the S2 PMF
# pipeline of the revised manuscript).
load_2d_grid <- function(ligand, ligand_label) {
  z <- np$load(file.path(data_dir, paste0("2d_c3_", ligand, ".npz")))
  xs <- as.numeric(z$f[["X"]])
  ys <- as.numeric(z$f[["Y"]])
  F  <- as.numeric(z$f[["F"]])
  df <- expand.grid(CV1 = xs, CV2 = ys)
  df$PMF <- F - min(F, na.rm = TRUE)
  df$ligand <- ligand_label
  df
}

# NOTE: a legacy weighted-histogram estimator (pre-DBE, 298 K constant) used to
# live here; it was retired 2026-09 and removed during the 2026-10-05 repository
# standardisation (all shipped landscapes use the DBE 300 K grids).

# ---- Compute for talazoparib and veliparib ----------------------------------
cat("Loading 2D C3 grids...\n")
pmf_tala <- load_2d_grid("talazoparib", "Talazoparib")
pmf_veli <- load_2d_grid("veliparib", "Veliparib")
pmf_apo  <- load_2d_grid("APO", "APO")

# ---- Panel A: Talazoparib 2D landscape -------------------------------------
p_a <- ggplot(pmf_tala, aes(x = CV1, y = CV2, fill = PMF, z = PMF)) +
  geom_raster() +
  geom_contour(color = "white", linewidth = 0.3, bins = 8, alpha = 0.6) +
  scale_fill_gradientn(colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B"),
                       name = "PMF (A,B,D)\nkcal/mol", limits = c(0, 8)) +
  labs(x = "CV1: Protein-DNA (Å)", y = "CV2: HD-ART (Å)",
       title = "Talazoparib (Type II)") +
  theme_7pt + coord_fixed()

# ---- Panel B: Veliparib 2D landscape ---------------------------------------
p_b <- ggplot(pmf_veli, aes(x = CV1, y = CV2, fill = PMF, z = PMF)) +
  geom_raster() +
  geom_contour(color = "white", linewidth = 0.3, bins = 8, alpha = 0.6) +
  scale_fill_gradientn(colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B"),
                       name = "PMF (A,B,D)\nkcal/mol", limits = c(0, 8)) +
  labs(x = "CV1: Protein-DNA (Å)", y = "CV2: HD-ART (Å)",
       title = "Veliparib (Type III)") +
  theme_7pt + coord_fixed()

# ---- Panel C: Difference map (Talazoparib - Veliparib) ----------------------
# Align grids
pmf_diff <- pmf_tala
pmf_diff$PMF <- pmf_tala$PMF - pmf_veli$PMF
pmf_diff$ligand <- "ΔPMF"

p_c <- ggplot(pmf_diff, aes(x = CV1, y = CV2, fill = PMF, z = PMF)) +
  geom_raster() +
  geom_contour(color = "grey40", linewidth = 0.3, bins = 8, alpha = 0.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       name = "ΔPMF (C)\nkcal/mol", midpoint = 0) +
  labs(x = "CV1: Protein-DNA (Å)", y = "CV2: HD-ART (Å)",
       title = "ΔPMF: Talazoparib − Veliparib") +
  theme_7pt + coord_fixed()

# ---- Panel D: APO (ligand-free control) 2D landscape ------------------------
p_d <- ggplot(pmf_apo, aes(x = CV1, y = CV2, fill = PMF, z = PMF)) +
  geom_raster() +
  geom_contour(color = "white", linewidth = 0.3, bins = 8, alpha = 0.6) +
  scale_fill_gradientn(colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B"),
                       name = "PMF (A,B,D)\nkcal/mol", limits = c(0, 8)) +
  labs(x = "CV1: Protein-DNA (Å)", y = "CV2: HD-ART (Å)",
       title = "APO (ligand-free control)") +
  theme_7pt + coord_fixed()

# ---- Tighten inter-row spacing (zero bottom margin on top row, zero top on bottom row) ----
p_a <- p_a + theme(plot.margin = ggplot2::margin(7, 7, 0, 7))
p_b <- p_b + theme(plot.margin = ggplot2::margin(7, 7, 0, 7))
p_c <- p_c + theme(plot.margin = ggplot2::margin(0, 7, 7, 7))
p_d <- p_d + theme(plot.margin = ggplot2::margin(0, 7, 7, 7))

# ---- Assemble (2x2: A|B top, C|D bottom) -----------------------------------
fig2d <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(guides = "collect")

cairo_pdf(file.path(out_dir, "Fig_2D_Landscape.pdf"), width = 157/25.4, height = 94/25.4, pointsize = 8)
print(fig2d)
dev.off()

cat("Saved: Fig_2D_Landscape.pdf\n")

# ---- Also compute for all systems (compact panels) --------------------------
all_systems <- c("APO", "AZD5305", "niraparib", "olaparib", "rucaparib", "veliparib")
all_pmf <- list()
for (sys in c("talazoparib", all_systems)) {
  cat(sprintf("Computing 2D PMF for %s...\n", sys))
  sys_name <- paste0("sys2_", sys)
  if (sys %in% c("APO", "talazoparib")) {
    # Already have weights files
  } else {
    # Check if .npy exists
    f <- file.path(data_dir, paste0(sys_name, "_CV1_cv.npy"))
    if (!file.exists(f)) {
      cat(sprintf("  Missing CV1 for %s, skipping\n", sys))
      next
    }
  }
  label <- ifelse(sys == "talazoparib", "Talazoparib",
           ifelse(sys == "APO", "APO",
           ifelse(sys == "AZD5305", "AZD5305",
           ifelse(sys == "niraparib", "Niraparib",
           ifelse(sys == "olaparib", "Olaparib",
           ifelse(sys == "rucaparib", "Rucaparib", "Veliparib"))))))
  tryCatch({
    df <- load_2d_grid(sys, label)
    all_pmf[[sys]] <- df
  }, error = function(e) {
    cat(sprintf("  Error for %s: %s\n", sys, e$message))
  })
}

if (length(all_pmf) > 0) {
  all_df <- bind_rows(all_pmf)
  all_df$ligand <- factor(all_df$ligand, 
                          levels = c("APO", "AZD5305", "Talazoparib", "Olaparib", 
                                     "Niraparib", "Rucaparib", "Veliparib"))
  
  p_all <- ggplot(all_df, aes(x = CV1, y = CV2, fill = PMF, z = PMF)) +
    geom_raster() +
    geom_contour(color = "white", linewidth = 0.2, bins = 6, alpha = 0.4) +
    scale_fill_gradientn(colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B"),
                         name = "kcal/mol") +
    facet_wrap(~ ligand, ncol = 4, scales = "free") +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    labs(x = "CV1: Protein-DNA (Å)", y = "CV2: HD-ART (Å)") +
    theme_7pt + theme(aspect.ratio = 1,
                      strip.text = element_text(size = 8),
                      plot.margin = ggplot2::margin(1, 1, 1, 3, unit = "mm"))
  
  # Square facet panels (aspect.ratio = 1 in the theme above); canvas sized so
  # each panel plots as a square at ~1:1 display in the SI (0.98\textwidth).
  cairo_pdf(file.path(out_dir, "Fig_2D_Landscape_All.pdf"), 
            width = 157/25.4, height = 88/25.4, pointsize = 8)
  print(p_all)
  dev.off()
  cat("Saved: Fig_2D_Landscape_All.pdf\n")
}

message("===== 2D Landscape complete =====")
