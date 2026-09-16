#!/usr/bin/env Rscript
# 16-2d_landscape.R — 2D free energy landscapes: talazoparib vs veliparib
# Uses CV1 (Protein-DNA) × CV2 (HD-ART) from .npy files + reweighting weights
# Requires RETICULATE_PYTHON pointing at a numpy-enabled Python, e.g.:
#   RETICULATE_PYTHON=/opt/anaconda3/bin/python3 Rscript code/16-2d_landscape.R
library(ggplot2)
library(dplyr)
library(MASS)  # for kde2d
library(reticulate)
library(patchwork)

np <- import("numpy")

data_dir  <- "results/analysis"
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_7pt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 7, face = "bold"))

# ---- Load CV data and compute 2D PMF ---------------------------------------
compute_2d_pmf <- function(system_name, ligand_label) {
  # Load CV1 and CV2 numpy arrays
  cv1 <- np$load(file.path(data_dir, paste0(system_name, "_CV1_cv.npy")))
  cv2 <- np$load(file.path(data_dir, paste0(system_name, "_CV2_cv.npy")))
  
  # Load weights
  weights_file <- file.path(data_dir, paste0(system_name, "_cv_weights.dat"))
  if (file.exists(weights_file)) {
    weights <- read.table(weights_file, header = FALSE)[, 1]
    # Subsample if needed to match
    if (length(weights) > length(cv1)) {
      weights <- weights[seq_len(length(cv1))]
    } else if (length(weights) < length(cv1)) {
      cv1 <- cv1[seq_len(length(weights))]
      cv2 <- cv2[seq_len(length(weights))]
    }
  } else {
    weights <- rep(1, length(cv1))
  }
  
  # Normalize weights
  weights <- weights / sum(weights)
  
  # Create 2D histogram with reweighting
  nbins <- 50
  x_breaks <- seq(min(cv1), max(cv1), length.out = nbins + 1)
  y_breaks <- seq(min(cv2), max(cv2), length.out = nbins + 1)
  
  # Compute weighted 2D histogram
  h2d <- matrix(0, nrow = nbins, ncol = nbins)
  for (i in seq_along(cv1)) {
    xi <- findInterval(cv1[i], x_breaks, all.inside = TRUE)
    yi <- findInterval(cv2[i], y_breaks, all.inside = TRUE)
    h2d[xi, yi] <- h2d[xi, yi] + weights[i]
  }
  
  # Convert to PMF (kcal/mol)
  h2d <- h2d / sum(h2d)
  pmf <- -0.592 * log(h2d + 1e-10)  # kT ≈ 0.592 at 298K
  pmf <- pmf - min(pmf, na.rm = TRUE)
  
  # Build data frame
  x_centers <- (x_breaks[-1] + x_breaks[-length(x_breaks)]) / 2
  y_centers <- (y_breaks[-1] + y_breaks[-length(y_breaks)]) / 2
  
  df <- expand.grid(CV1 = x_centers, CV2 = y_centers)
  df$PMF <- as.vector(pmf)
  df$ligand <- ligand_label
  
  df
}

# ---- Compute for talazoparib and veliparib ----------------------------------
cat("Computing 2D PMF for talazoparib...\n")
pmf_tala <- compute_2d_pmf("sys2_talazoparib", "Talazoparib")
cat("Computing 2D PMF for veliparib...\n")
pmf_veli <- compute_2d_pmf("sys2_veliparib", "Veliparib")

# Also APO for reference
cat("Computing 2D PMF for APO...\n")
pmf_apo <- compute_2d_pmf("sys2_APO", "APO")

# ---- Panel A: Talazoparib 2D landscape -------------------------------------
p_a <- ggplot(pmf_tala, aes(x = CV1, y = CV2, fill = PMF, z = PMF)) +
  geom_raster() +
  geom_contour(color = "white", linewidth = 0.3, bins = 8, alpha = 0.6) +
  scale_fill_gradientn(colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B"),
                       name = "kcal/mol", limits = c(0, 8)) +
  labs(x = "CV1: Protein-DNA Distance (Å)", y = "CV2: HD-ART Distance (Å)",
       title = "Talazoparib (Type II, trapping 100× olaparib)") +
  theme_7pt + coord_fixed()

# ---- Panel B: Veliparib 2D landscape ---------------------------------------
p_b <- ggplot(pmf_veli, aes(x = CV1, y = CV2, fill = PMF, z = PMF)) +
  geom_raster() +
  geom_contour(color = "white", linewidth = 0.3, bins = 8, alpha = 0.6) +
  scale_fill_gradientn(colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B"),
                       name = "kcal/mol", limits = c(0, 8)) +
  labs(x = "CV1: Protein-DNA Distance (Å)", y = "CV2: HD-ART Distance (Å)",
       title = "Veliparib (Type III, trapping 0.02× olaparib)") +
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
                       name = "Δ kcal/mol", midpoint = 0) +
  labs(x = "CV1: Protein-DNA Distance (Å)", y = "CV2: HD-ART Distance (Å)",
       title = "ΔPMF: Talazoparib − Veliparib") +
  theme_7pt + coord_fixed()

# ---- Assemble --------------------------------------------------------------
fig2d <- (p_a | p_b | p_c) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "2D Free Energy Landscapes: Same Pocket, Different Encoding",
                  theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)))

cairo_pdf(file.path(out_dir, "Fig_2D_Landscape.pdf"), width = 210/25.4, height = 80/25.4, pointsize = 7)
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
    df <- compute_2d_pmf(sys_name, label)
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
    labs(x = "CV1: Protein-DNA (Å)", y = "CV2: HD-ART (Å)",
         title = "S2: 2D Free Energy Landscapes — All Systems") +
    theme_7pt + theme(strip.text = element_text(size = 6))
  
  cairo_pdf(file.path(out_dir, "Fig_2D_Landscape_All.pdf"), 
            width = 210/25.4, height = 140/25.4, pointsize = 7)
  print(p_all)
  dev.off()
  cat("Saved: Fig_2D_Landscape_All.pdf\n")
}

message("===== 2D Landscape complete =====")
