#!/usr/bin/env Rscript
# 01-pmf_apo.R — Generate PMF plots for sys2_APO GaMD reweighting
#
# Purpose:  Regenerate the S2 APO control PMF panels (protein-DNA CV1 and HD-ART CV2) from the reweighted C3 curves.
# Author:   Tao Zhu (tjogzt@gmail.com)
# Created:  2026-09-15 (header standardised 2026-10-05)
# Inputs:   data/analysis/pmf-c3-sys2_APO_CV1_cv.dat.xvg ; data/analysis/pmf-c3-sys2_APO_CV2_cv.dat.xvg
# Outputs:  figures/pdf/01-pmf_apo_hd_art.pdf ; figures/pdf/01-pmf_apo_prot_dna.pdf
# Depends:  ggplot2
# Run:      Rscript code/01_pmf_apo.R   (from the repository root)
library(ggplot2)

# Helper: read xvg (skip header lines starting with @ or #)
read_xvg <- function(path) {
  lines <- readLines(path)
  # skip header, keep data
  data_lines <- lines[!grepl("^[@#]", lines)]
  # remove empty lines
  data_lines <- data_lines[data_lines != ""]
  writeLines(data_lines, con = tempfile())
  read.table(text = data_lines, header = FALSE)
}

# Read PMF files (c3 = 3rd order cumulant expansion, most accurate)
pmf_cv1 <- read_xvg("data/analysis/pmf-c3-sys2_APO_CV1_cv.dat.xvg")
pmf_cv2 <- read_xvg("data/analysis/pmf-c3-sys2_APO_CV2_cv.dat.xvg")
colnames(pmf_cv1) <- c("RC", "PMF")
colnames(pmf_cv2) <- c("RC", "PMF")

# Set PMF minimum to 0
pmf_cv1$PMF <- pmf_cv1$PMF - min(pmf_cv1$PMF)
pmf_cv2$PMF <- pmf_cv2$PMF - min(pmf_cv2$PMF)

# CV1: Protein-DNA distance
p1 <- ggplot(pmf_cv1, aes(x = RC, y = PMF)) +
  geom_line(linewidth = 0.8, color = "#2166AC") +
  labs(x = "Protein-DNA COM Distance (Å)", 
       y = "Free Energy (kcal/mol)",
       title = "sys2_APO — Protein-DNA Distance PMF") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.minor = element_blank())

# CV2: HD-ART distance
p2 <- ggplot(pmf_cv2, aes(x = RC, y = PMF)) +
  geom_line(linewidth = 0.8, color = "#B2182B") +
  labs(x = "HD-ART COM Distance (Å)", 
       y = "Free Energy (kcal/mol)",
       title = "sys2_APO — HD-ART Distance PMF") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.minor = element_blank())

# Save as PDF (SCI standard: cairo_pdf, 7pt)
dir.create("figures/pdf", showWarnings = FALSE, recursive = TRUE)

cairo_pdf("figures/pdf/01-pmf_apo_prot_dna.pdf", 
          width = 3.5, height = 2.5, pointsize = 8)
print(p1)
dev.off()

cairo_pdf("figures/pdf/01-pmf_apo_hd_art.pdf", 
          width = 3.5, height = 2.5, pointsize = 8)
print(p2)
dev.off()

# Combined panel figure
p_combined <- cowplot::plot_grid(p1, p2, labels = c("A", "B"), 
                                  label_size = 8, ncol = 2)

cairo_pdf("figures/pdf/01-pmf_apo_combined.pdf", 
          width = 7, height = 3, pointsize = 8)
print(p_combined)
dev.off()

cat("[OK] PMF figures saved to figures/pdf/\n")
cat("CV1 (Prot-DNA): range =", range(pmf_cv1$RC), "Å\n")
cat("CV2 (HD-ART): range =", range(pmf_cv2$RC), "Å\n")
