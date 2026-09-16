#!/usr/bin/env Rscript
# 02-pmf_comparison.R — APO vs talazoparib PMF comparison (sys2)
library(ggplot2)

read_xvg <- function(path) {
  lines <- readLines(path)
  data_lines <- lines[!grepl("^[@#]", lines)]
  data_lines <- data_lines[data_lines != ""]
  read.table(text = data_lines, header = FALSE)
}

# Read PMF files (c3 = 3rd order cumulant expansion)
pmf_apo_cv1 <- read_xvg("results/analysis/pmf-c3-sys2_APO_CV1_cv.dat.xvg")
pmf_apo_cv2 <- read_xvg("results/analysis/pmf-c3-sys2_APO_CV2_cv.dat.xvg")
pmf_tlz_cv1 <- read_xvg("results/analysis/pmf-c3-sys2_talazoparib_CV1_cv.dat.xvg")
pmf_tlz_cv2 <- read_xvg("results/analysis/pmf-c3-sys2_talazoparib_CV2_cv.dat.xvg")

colnames(pmf_apo_cv1) <- c("RC", "PMF")
colnames(pmf_apo_cv2) <- c("RC", "PMF")
colnames(pmf_tlz_cv1) <- c("RC", "PMF")
colnames(pmf_tlz_cv2) <- c("RC", "PMF")

# Set PMF minimum to 0
pmf_apo_cv1$PMF <- pmf_apo_cv1$PMF - min(pmf_apo_cv1$PMF)
pmf_apo_cv2$PMF <- pmf_apo_cv2$PMF - min(pmf_apo_cv2$PMF)
pmf_tlz_cv1$PMF <- pmf_tlz_cv1$PMF - min(pmf_tlz_cv1$PMF)
pmf_tlz_cv2$PMF <- pmf_tlz_cv2$PMF - min(pmf_tlz_cv2$PMF)

# Add system labels
pmf_apo_cv1$System <- "APO"
pmf_apo_cv2$System <- "APO"
pmf_tlz_cv1$System <- "Talazoparib"
pmf_tlz_cv2$System <- "Talazoparib"

# Combine
pmf_cv1 <- rbind(pmf_apo_cv1, pmf_tlz_cv1)
pmf_cv2 <- rbind(pmf_apo_cv2, pmf_tlz_cv2)

# Color scheme
sys_colors <- c("APO" = "#2166AC", "Talazoparib" = "#B2182B")

# CV1: Protein-DNA distance
p1 <- ggplot(pmf_cv1, aes(x = RC, y = PMF, color = System)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = sys_colors) +
  labs(x = "Protein-DNA COM Distance (Å)", 
       y = "Free Energy (kcal/mol)",
       title = "PARP1-DNA Retention PMF") +
  theme_bw(base_size = 7) +
  theme(plot.title = element_text(size = 7, face = "bold"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),
        panel.grid.minor = element_blank())

# CV2: HD-ART distance
p2 <- ggplot(pmf_cv2, aes(x = RC, y = PMF, color = System)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = sys_colors) +
  labs(x = "HD-ART COM Distance (Å)", 
       y = "Free Energy (kcal/mol)",
       title = "HD-ART Domain Motion PMF") +
  theme_bw(base_size = 7) +
  theme(plot.title = element_text(size = 7, face = "bold"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),
        panel.grid.minor = element_blank())

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Individual panels
cairo_pdf("results/figures/02-pmf_comparison_prot_dna.pdf", 
          width = 3.5, height = 2.5, pointsize = 7)
print(p1)
dev.off()

cairo_pdf("results/figures/02-pmf_comparison_hd_art.pdf", 
          width = 3.5, height = 2.5, pointsize = 7)
print(p2)
dev.off()

# Combined panel
p_combined <- cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                                  p2 + theme(legend.position = "none"),
                                  labels = c("A", "B"), 
                                  label_size = 7, ncol = 2)
# Extract legend
legend <- cowplot::get_legend(p1)
p_final <- cowplot::plot_grid(p_combined, legend, rel_widths = c(1, 0.15))

cairo_pdf("results/figures/02-pmf_comparison_combined.pdf", 
          width = 7, height = 3, pointsize = 7)
print(p_final)
dev.off()

cat("[OK] PMF comparison figures saved to results/figures/\n")
cat("CV1 (Prot-DNA) APO: range =", range(pmf_apo_cv1$RC), "Å\n")
cat("CV1 (Prot-DNA) TLZ: range =", range(pmf_tlz_cv1$RC), "Å\n")
cat("CV2 (HD-ART) APO: range =", range(pmf_apo_cv2$RC), "Å\n")
cat("CV2 (HD-ART) TLZ: range =", range(pmf_tlz_cv2$RC), "Å\n")
