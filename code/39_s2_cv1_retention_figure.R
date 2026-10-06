#!/usr/bin/env Rscript
# Regenerate Fig_S2_CV1_Retention.pdf — S2 CV1 (protein–DNA) retention summary
# =============================================================================
# (A) 7-system CV1 PMF overlay (C3 cumulant, per-curve min-normalized)
#     Source: data/analysis/pmf-c3-sys2_<sys>_CV1_cv.dat.xvg
#     (niraparib file restored 2026-10-05 from archived neutral per-frame data —
#      validated to the official span 56.64927809233839)
# (B) CV1 sampled-range spans (C1 vs C3) for the seven original-panel systems
#     Values: data/01_curated/s2cv1_spans_table_s5.csv (Table S5 convention).
#     Palette: common/palette_si.csv (single source).
#
# 2026-10-05 fix: the former fread()-based reader silently parsed the xvg files
# into all-NA values, so ggplot dropped every point and rendered an empty frame.
# This version uses a robust line-based reader (same convention as
# common/helpers.R::read_pmf) so a future format drift can never empty the plot
# without an error.
# =============================================================================
#
# Purpose:  Render Fig_S2_CV1_Retention.pdf: S2 CV1 (protein–DNA) retention summary with the C1/C3 span panel.
# Created:  2026-09-15 (header standardised 2026-10-05)
# Outputs:  figures/pdf/Fig_S2_CV1_Retention.pdf
# Depends:  ggplot2, patchwork
# Run:      Rscript code/39_s2_cv1_retention_figure.R   (from the repository root)
set.seed(49)
suppressMessages({library(ggplot2); library(patchwork)})

theme_set(theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
        panel.border = element_blank(),                      # open frame:
        axis.line = element_line(color = "black", linewidth = 0.5),  # left+bottom only
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)))

dir <- "data/analysis"
sys_names <- c("APO", "talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib", "AZD5305")
short <- c(APO = "APO", talazoparib = "Talazoparib", niraparib = "Niraparib",
           olaparib = "Olaparib", rucaparib = "Rucaparib",
           veliparib = "Veliparib", AZD5305 = "AZD5305")
pal <- read.csv("common/palette_si.csv", stringsAsFactors = FALSE)
cols <- setNames(pal$color, pal$system)

read_pmf <- function(s) {
  f <- file.path(dir, paste0("pmf-c3-sys2_", s, "_CV1_cv.dat.xvg"))
  lines <- readLines(f)
  lines <- lines[grepl("^\\s*[0-9]", lines)]
  nums <- strsplit(trimws(lines), "[ \t]+")
  rc  <- as.numeric(vapply(nums, `[`, "", 1))
  pmf <- as.numeric(vapply(nums, `[`, "", 2))
  keep <- is.finite(rc) & is.finite(pmf)
  if (sum(keep) < 5) stop("read_pmf: too few valid rows in ", f)
  data.frame(rc = rc[keep], pmf = pmf[keep] - min(pmf[keep]), system = s)
}
pmf_all <- do.call(rbind, lapply(sys_names, read_pmf))
pmf_all$system <- factor(pmf_all$system, levels = sys_names)

pA <- ggplot(pmf_all, aes(rc, pmf, color = system)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = cols, labels = short) +
  labs(x = "CV1: Protein–DNA COM Distance (\u00C5)", y = "PMF (kcal/mol)",
       title = "A  S2 CV1 PMF overlay (C3 cumulant)", color = NULL) +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 9),
        legend.key.size = unit(0.25, "cm"))

# ---- Panel B: CV1 sampled-range spans (C1 vs C3), Table S5 values ----------
# Values live in data/01_curated/s2cv1_spans_table_s5.csv (Table S5 convention).
sp <- read.csv("data/01_curated/s2cv1_spans_table_s5.csv", stringsAsFactors = FALSE)
bars <- rbind(
  data.frame(system = sp$label, order = "C1", span = sp$C1, stringsAsFactors = FALSE),
  data.frame(system = sp$label, order = "C3", span = sp$C3, stringsAsFactors = FALSE))
bars$system <- factor(bars$system, levels = c("APO", "Talazoparib", "Niraparib", "Olaparib",
                                              "Rucaparib", "Veliparib", "AZD5305"))
sys_cols_bar <- unname(cols[c("APO", "talazoparib", "niraparib", "olaparib",
                              "rucaparib", "veliparib", "AZD5305")])
names(sys_cols_bar) <- levels(bars$system)
bars$fillc <- ifelse(bars$order == "C1", "#B5B5B5", sys_cols_bar[as.character(bars$system)])
# manual dodge: C1 bar at x-0.2, C3 bar at x+0.2 (0.36 width each)
bars$xpos <- as.numeric(bars$system) + ifelse(bars$order == "C1", -0.2, 0.2)
lbl <- bars[bars$order == "C3", ]

pB <- ggplot(bars, aes(xpos, span)) +
  geom_col(aes(fill = fillc), width = 0.36) +
  geom_text(data = lbl, aes(x = xpos, y = span + 1.8, label = sprintf("%.1f", span)),
            size = 2.9, vjust = 0) +
  scale_fill_identity() +
  scale_x_continuous(breaks = 1:7, labels = levels(bars$system)) +
  coord_cartesian(ylim = c(0, 90)) +
  labs(x = NULL, y = "CV1 span (kcal/mol)",
       title = "B  CV1 sampled-range spans: near-uniform C1 (gray), larger C3 variation") +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 9),
        axis.text.x = element_text(angle = 30, hjust = 1))

p <- pA / pB + plot_layout(heights = c(1.25, 1))
cairo_pdf("figures/pdf/Fig_S2_CV1_Retention.pdf",
          width = 149/25.4, height = 95.2/25.4, pointsize = 8)
print(p)
invisible(dev.off())

# self-check: panel A must carry real data for all 7 systems
chk <- tapply(pmf_all$pmf, pmf_all$system, function(v) max(v) - min(v))
stopifnot(length(chk) == 7, all(chk > 30))
cat("Saved: figures/pdf/Fig_S2_CV1_Retention.pdf\n")
cat("Panel A per-system C3 spans rebuilt from files:\n")
print(round(chk, 1))
