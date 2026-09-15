#!/usr/bin/env Rscript
# 01-pmf_apo.R — Generate PMF plots for sys2_APO GaMD reweighting
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
pmf_cv1 <- read_xvg("results/analysis/pmf-c3-sys2_APO_CV1_cv.dat.xvg")
pmf_cv2 <- read_xvg("results/analysis/pmf-c3-sys2_APO_CV2_cv.dat.xvg")
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
  theme_bw(base_size = 7) +
  theme(plot.title = element_text(size = 7, face = "bold"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.grid.minor = element_blank())

# CV2: HD-ART distance
p2 <- ggplot(pmf_cv2, aes(x = RC, y = PMF)) +
  geom_line(linewidth = 0.8, color = "#B2182B") +
  labs(x = "HD-ART COM Distance (Å)", 
       y = "Free Energy (kcal/mol)",
       title = "sys2_APO — HD-ART Distance PMF") +
  theme_bw(base_size = 7) +
  theme(plot.title = element_text(size = 7, face = "bold"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.grid.minor = element_blank())

# Save as PDF (SCI standard: cairo_pdf, 7pt)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

cairo_pdf("results/figures/01-pmf_apo_prot_dna.pdf", 
          width = 3.5, height = 2.5, pointsize = 7)
print(p1)
dev.off()

cairo_pdf("results/figures/01-pmf_apo_hd_art.pdf", 
          width = 3.5, height = 2.5, pointsize = 7)
print(p2)
dev.off()

# Combined panel figure
p_combined <- cowplot::plot_grid(p1, p2, labels = c("A", "B"), 
                                  label_size = 7, ncol = 2)

cairo_pdf("results/figures/01-pmf_apo_combined.pdf", 
          width = 7, height = 3, pointsize = 7)
print(p_combined)
dev.off()

cat("[OK] PMF figures saved to results/figures/\n")
cat("CV1 (Prot-DNA): range =", range(pmf_cv1$RC), "Å\n")
cat("CV2 (HD-ART): range =", range(pmf_cv2$RC), "Å\n")
