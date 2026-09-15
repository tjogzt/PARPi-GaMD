#!/usr/bin/env Rscript
# 重画 Fig_DCCM_Allostery.pdf — 7 体系 HD×ART DCCM 块热图 (全蛋白超叠, 统一协议)
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

theme_set(theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
        legend.key.height = unit(0.35, "cm"), legend.key.width = unit(0.5, "cm")))

dir <- "results/analysis/rmsf_recomp"
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
    theme(axis.text = element_text(size = 6))
})

p <- wrap_plots(plots, ncol = 3)
cairo_pdf("results/figures/Fig_DCCM_Allostery.pdf", width = 7.2, height = 5.2, pointsize = 8)
print(p)
dev.off()
cat("Saved: results/figures/Fig_DCCM_Allostery.pdf\n")
