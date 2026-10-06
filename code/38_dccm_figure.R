#!/usr/bin/env Rscript
# Regenerate Fig_DCCM_Allostery.pdf — 7-system HDxART DCCM block heatmap (whole-protein superposition, unified protocol)
#
# Output: vector PDF (figures/pdf/Fig_DCCM_Allostery.pdf) — canonical, text-searchable.
# Fig_DCCM_Allostery.pdf: 300-dpi raster, compiled for the Supporting Information.
# See the SI figure build section in README for the generation workflow.
#
# Purpose:  Render Fig_DCCM_Allostery.pdf: the seven-system HDxART DCCM block heatmap.
# Created:  2026-09-15 (header standardised 2026-10-05)
# Depends:  data.table, ggplot2, patchwork
# Run:      Rscript code/38_dccm_figure.R   (from the repository root)
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

theme_set(theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        legend.key.height = unit(0.35, "cm"), legend.key.width = unit(0.5, "cm")))

dir <- "data/analysis/rmsf_recomp"
sys_names <- c("APO", "talazoparib", "olaparib", "niraparib", "rucaparib", "veliparib", "AZD5305")
short <- c(APO = "APO", talazoparib = "Talazoparib", olaparib = "Olaparib",
           niraparib = "Niraparib", rucaparib = "Rucaparib",
           veliparib = "Veliparib", AZD5305 = "AZD5305")

plots <- lapply(sys_names, function(s) {
  m <- as.matrix(fread(file.path(dir, paste0(s, "_dccm_block.csv"))))
  art_ids <- as.numeric(colnames(m))
  hd_ids <- 662:787
  dt <- CJ(hd = hd_ids, art = art_ids)
  dt[, r := as.vector(m)]
  ggplot(dt, aes(art, hd, fill = r)) +
    geom_tile() +
    scale_fill_gradient2(low = "#3D6BA8", mid = "white", high = "#C23531",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    scale_y_reverse() +
    labs(x = "ART residue", y = "HD residue", title = short[[s]]) +
    theme(axis.text = element_text(size = 8))
})

p <- wrap_plots(plots, ncol = 3)
cairo_pdf("figures/pdf/Fig_DCCM_Allostery.pdf", width = 149/25.4, height = 107.6/25.4, pointsize = 8)
print(p)
dev.off()
cat("Saved: figures/pdf/Fig_DCCM_Allostery.pdf\n")
