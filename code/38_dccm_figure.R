#!/usr/bin/env Rscript
# Regenerate Fig_DCCM_Allostery.pdf — 7-system HDxART DCCM block heatmap (whole-protein superposition, unified protocol)
#
# Output: vector PDF (figures/pdf/Fig_DCCM_Allostery.pdf) — canonical, text-searchable.
# SI compilation copy: manuscripts/jcim/si_figures/Fig_DCCM_Allostery.pdf is a 300-dpi raster
# version of the same figure (keeps the combined SI PDF under the 2 MB submission limit);
# regenerate it with scripts/rasterize_si_figures.sh after any change to this script.
#
# Purpose:  Render the domain cross-correlation matrix (DCCM) allostery figure (SI Figure S13).
# Created:  2026-10-02 (header standardised 2026-10-05; R6 rebuild 2026-10-10)
# Inputs:   data/analysis/rmsf_recomp/*_dccm_block.csv
# Outputs:  figures/pdf/Fig_DCCM_Allostery.pdf
# Depends:  data.table, ggplot2, patchwork
# Run:      Rscript code/38_dccm_figure.R   (from the repository root)
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

theme_set(theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        # (2026-10-15) full-width layout: sizes tuned for the 178-mm canvas
        # displayed at \textwidth in the SI; narrow color bar (0.25 cm).
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 9),
        legend.key.height = unit(0.35, "cm"), legend.key.width = unit(0.25, "cm"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9)))

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
  # (2026-10-15) right-side band: per-HD-residue mean correlation (mean over
  # the ART axis), drawn in the same coordinate system so it stays row-aligned
  # with the matrix and renders the HD-side modulation directly legible.
  strip <- data.frame(hd = hd_ids, rmean = rowMeans(m))
  ggplot(dt, aes(art, hd)) +
    geom_tile(aes(fill = r)) +
    geom_rect(data = strip, inherit.aes = FALSE,
              aes(xmin = 1018, xmax = 1030, ymin = hd - 0.5, ymax = hd + 0.5,
                  fill = rmean)) +
    scale_fill_gradient2(low = "#3D6BA8", mid = "white", high = "#C23531",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    scale_y_reverse() +
    scale_x_continuous(breaks = seq(800, 1000, 100)) +
    labs(x = "ART residue", y = "HD residue", title = short[[s]]) +
    guides(fill = guide_colorbar(barwidth = unit(2.5, "mm"),
                                 barheight = unit(30, "mm"))) +
    theme(axis.text = element_text(size = 9), aspect.ratio = 1)
})

# (2026-10-15) merged legend: the seven per-panel color bars are identical
# (same scale, same limits) -- a single shared bar is placed in the vacant
# eighth grid cell (row 3, column 2) so that all seven matrices use the full
# available width.
g_ref <- ggplot2::ggplotGrob(plots[[1]])
leg_shared <- g_ref$grobs[[which(vapply(g_ref$grobs, function(x) x$name, "") == "guide-box")]]
plots <- lapply(plots, function(p) p + theme(legend.position = "none"))

p <- wrap_plots(c(plots, list(leg_shared)), ncol = 3)
cairo_pdf("figures/pdf/Fig_DCCM_Allostery.pdf", width = 178/25.4, height = 189/25.4, pointsize = 9)
print(p)
dev.off()
cat("Saved: figures/pdf/Fig_DCCM_Allostery.pdf\n")
