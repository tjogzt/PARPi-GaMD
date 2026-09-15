#!/usr/bin/env Rscript
# 重画 Fig_S2_CV1_Retention.pdf — C3 cumulant 管线 CV1 (蛋白-DNA) PMF
# (Top) 7 体系 PMF 叠加  (Bottom) 井深柱状 (几乎不变 27.5-30.0)
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

theme_set(theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank()))

dir <- "results/analysis"
sys_names <- c("APO", "talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib", "AZD5305")
short <- c(APO = "APO", talazoparib = "Talazoparib", niraparib = "Niraparib",
           olaparib = "Olaparib", rucaparib = "Rucaparib",
           veliparib = "Veliparib", AZD5305 = "AZD5305")
cols <- c(APO = "#7A7A7A", talazoparib = "#C23531", niraparib = "#9D2933",
          olaparib = "#3D6BA8", rucaparib = "#177CB0",
          veliparib = "#B36B2E", AZD5305 = "#5E8C5E")

read_pmf <- function(s) {
  f <- file.path(dir, paste0("pmf-c3-sys2_", s, "_CV1_cv.dat.xvg"))
  x <- fread(f, skip = 0, fill = TRUE)
  x <- x[!grepl("^[#@]", V1)]
  x[, c("V3", "V4") := NULL]
  x[, V1 := as.numeric(V1)]
  x[, V2 := as.numeric(V2)]
  setnames(x, c("rc", "pmf"))
  x[, pmf := pmf - min(pmf)]
  x[, system := s]
  x
}
pmf_all <- rbindlist(lapply(sys_names, read_pmf))
pmf_all[, system := factor(system, levels = sys_names)]

wells <- pmf_all[, .(wd = max(pmf)), by = system]
wells[, label := short[as.character(system)]]
wells[, label := factor(label, levels = sapply(sys_names, function(s) short[[s]]))]

pA <- ggplot(pmf_all, aes(rc, pmf, color = system)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = cols, labels = short) +
  labs(x = "CV1: protein\u2013DNA COM distance (\u00C5)", y = "PMF (kcal/mol)",
       title = "A  S2 CV1 PMF overlay (C3 cumulant)", color = NULL) +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 9),
        legend.key.size = unit(0.25, "cm"))

pB <- ggplot(wells, aes(label, wd, fill = system)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = cols, guide = "none") +
  geom_text(aes(y = wd + 0.5, label = sprintf("%.1f", wd)), size = 2.4) +
  coord_cartesian(ylim = c(0, 34)) +
  labs(x = NULL, y = "CV1 well depth (kcal/mol)",
       title = "B  Protein\u2013DNA retention barrier is nearly invariant") +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 9),
        axis.text.x = element_text(angle = 30, hjust = 1))

p <- pA / pB + plot_layout(heights = c(1.3, 1))
cairo_pdf("results/figures/Fig_S2_CV1_Retention.pdf", width = 7.2, height = 4.6, pointsize = 8)
print(p)
dev.off()
cat("Saved: results/figures/Fig_S2_CV1_Retention.pdf\n")
