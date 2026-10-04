#!/usr/bin/env Rscript
# Fig 1: PARP1 domain architecture + system schematics + inhibitor panel
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Shared helpers: theme_7pt (single source).
args_h <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", args_h))) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
  source(file.path(script_dir, "..", "common", "helpers.R"))
} else {
  source("common/helpers.R")
}

# ---- Panel A: Domain Architecture Schematic ----
domain_data <- data.frame(
  domain   = c("ZnF1", "ZnF2", "ZnF3", "BRCT", "WGR", "HD", "ART"),
  start    = c(1, 97, 215, 384, 518, 662, 788),
  end      = c(96, 214, 383, 517, 661, 787, 1014),
  color    = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"),
  y        = 1,
  stringsAsFactors = FALSE
)
domain_data$domain <- factor(domain_data$domain, levels = domain_data$domain)

# Annotation for binding sites
annotations <- data.frame(
  domain = c("ZnF1", "ZnF3", "WGR", "HD", "ART", "ART"),
  pos    = c(50, 300, 590, 725, 870, 870),
  label  = c("DNA\nbinding", "DNA\nbinding", "DNA\nbinding", "Type I\n(EB-47)", "CAT Pocket\nType II/III", "(this study)"),
  y      = c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3),
  stringsAsFactors = FALSE
)

p_a <- ggplot() +
  geom_rect(data = domain_data, aes(xmin = start, xmax = end, ymin = y - 0.25, ymax = y + 0.25, fill = domain),
            color = "grey30", linewidth = 0.3) +
  geom_text(data = domain_data, aes(x = (start + end) / 2, y = y, label = domain), size = 3.0, fontface = "bold") +
  geom_segment(aes(x = 50, xend = 50, y = 1.25, yend = 1.60), linewidth = 0.25, color = "grey45") +
  geom_segment(aes(x = 300, xend = 300, y = 1.25, yend = 1.60), linewidth = 0.25, color = "grey45") +
  geom_segment(aes(x = 590, xend = 590, y = 1.25, yend = 1.60), linewidth = 0.25, color = "grey45") +
  geom_segment(aes(x = 725, xend = 725, y = 1.25, yend = 1.60), linewidth = 0.25, color = "grey45") +
  geom_segment(aes(x = 870, xend = 870, y = 1.25, yend = 1.60), linewidth = 0.25, color = "grey45") +
  annotate("text", x = 12, y = 1.69, label = "DNA\nbinding", size = 3.0, color = "grey30", lineheight = 0.85, hjust = 0, vjust = 0) +
  annotate("text", x = 210, y = 1.69, label = "DNA\nbinding", size = 3.0, color = "grey30", lineheight = 0.85, hjust = 0, vjust = 0) +
  annotate("text", x = 500, y = 1.69, label = "DNA\nbinding", size = 3.0, color = "grey30", lineheight = 0.85, hjust = 0, vjust = 0) +
  annotate("text", x = 640, y = 1.69, label = "Type I\n(EB-47)", size = 3.0, color = "grey30", lineheight = 0.85, hjust = 0, vjust = 0) +
  annotate("text", x = 780, y = 1.69, label = "CAT Pocket\nType II/III\n(this study)", size = 3.0, color = "grey30", lineheight = 0.85, hjust = 0, vjust = 0) +
  scale_fill_manual(values = setNames(domain_data$color, domain_data$domain), guide = "none") +
  scale_x_continuous(limits = c(0, 1020), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.72, 2.02)) +
  labs(x = "Residue Number", y = NULL,
       title = "A  PARP1 Domain Architecture & Binding Sites") +
  theme_7pt + theme(panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.line.x = element_line(linewidth = 0.5, color = "black"),
                    axis.line.y = element_blank(),
                    axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ---- Panel C (rebuilt, four-team pass): minimal schematic, boxes scaled by system size ----
# Box areas scale with atom counts (~61.5k vs ~286k; width ratio sqrt(286/61.5) = 2.16);
# CV definitions, conditions and replicates live in the caption (no duplication).
p_b <- ggplot() +
  annotate("rect", xmin = 0.07, xmax = 1.93, ymin = 0.35, ymax = 1.45,
           fill = "#B2182B", alpha = 0.08, color = "#B2182B", linewidth = 0.5) +
  annotate("rect", xmin = 0.545, xmax = 1.455, ymin = 1.80, ymax = 2.30,
           fill = "#2166AC", alpha = 0.10, color = "#2166AC", linewidth = 0.5) +
  annotate("text", x = 1, y = 2.17, label = "S1: CAT-only", size = 3.2,
           fontface = "bold", color = "#2166AC", family = "Arial") +
  annotate("text", x = 1, y = 1.99, label = "~61.5k atoms \u00b7 19\u2013200 ns", size = 3.0,
           color = "grey30", family = "Arial") +
  annotate("text", x = 1, y = 1.05, label = "S2: PARP1 + DNA", size = 3.2,
           fontface = "bold", color = "#B2182B", family = "Arial") +
  annotate("text", x = 1, y = 0.83, label = "~286k atoms \u00b7 22\u201326 ns", size = 3.0,
           color = "grey30", family = "Arial") +
  annotate("segment", x = 1, xend = 1, y = 1.78, yend = 1.47,
           arrow = arrow(type = "closed", length = unit(0.08, "inches")),
           color = "grey50", linewidth = 0.5) +
  annotate("text", x = 1.08, y = 1.62, label = "AAI = S1/S2", size = 3.0,
           color = "grey30", hjust = 0, family = "Arial") +
  xlim(0, 2.0) + ylim(0.25, 2.40) +
  labs(title = "B  GaMD Dual-System Design") +
  theme_7pt + theme(panel.border = element_blank(),
                    axis.text = element_blank(), axis.ticks = element_blank(),
                    axis.title = element_blank(), panel.grid = element_blank())

# ---- Panel C: Inhibitor Structures + Trapping Data ----
# Canonical sources: data/01_curated/trapping_potency.csv (x olaparib, PMIDs in file);
# class from common/ligands.csv. AZD5305: potent PARP1-selective trapper measured
# on a different scale (Illuzzi 2022, PMID 35929986) -- not plotted on the
# dual-PARP x olaparib axis.
trap_csv <- fread("data/01_curated/trapping_potency.csv")
cls_csv  <- fread("common/ligands.csv")
cls_map  <- setNames(cls_csv$class, cls_csv$ligand)

inhib_names <- c("Talazoparib", "Niraparib", "Olaparib", "Rucaparib", "Veliparib")
trap_vals <- setNames(trap_csv$trapping_x_olaparib, trap_csv$inhibitor)[tolower(inhib_names)]
cls_vals <- cls_map[c("talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib")]

# Numeric bar heights for the log-scale panel. Niraparib has no quantitative
# fold-change in the open-access sources, so its bar is placed at 10x (rank-only,
# labelled ">1"). AZD5305 is annotated separately (see annotate below).
trap_num <- c(talazoparib = 100, niraparib = 10, olaparib = 1,
              rucaparib = 1, veliparib = 0.1)

inhib_data <- data.frame(
  ligand       = factor(inhib_names, levels = inhib_names),
  trap_potency = unname(trap_num),
  type         = gsub("_", " ", unname(cls_vals)),
  bar_label    = unname(trap_vals),
  stringsAsFactors = FALSE
)

# niraparib enters by rank only: its bar is drawn hollow (non-quantitative
# placeholder height) and labelled accordingly.
inhib_data$bar_label[inhib_data$ligand == "Niraparib"] <- ">1\n(rank only)"
nira <- inhib_data[inhib_data$ligand == "Niraparib", ]

p_c <- ggplot(inhib_data, aes(x = ligand, y = trap_potency)) +
  geom_bar(data = inhib_data[inhib_data$ligand != "Niraparib", ],
           aes(fill = type), stat = "identity", width = 0.7,
           color = "grey30", linewidth = 0.2) +
  geom_bar(data = nira, stat = "identity", width = 0.7,
           fill = "white", color = "#377EB8", linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(label = bar_label, y = trap_potency * 1.5),
            size = 3.0, vjust = 0, lineheight = 0.9) +
  scale_fill_manual(values = c("Type II" = "#E41A1C", "Type III" = "#377EB8", "Unknown" = "darkorange"),
                    breaks = c("Type II", "Type III")) +
  guides(fill = guide_legend(title = "HX-MS Type", nrow = 1)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
  labs(x = NULL, y = "Trapping Potency (× Olaparib, log scale)",
       title = "C  Trapping Potency &\nType Classification",
       subtitle = "AZD5305: PARP1-selective, off the dual-PARP scale") +
  theme_7pt + theme(panel.border = element_blank(),
                    axis.line = element_line(linewidth = 0.5, color = "black"),
                    panel.grid.major.x = element_blank(),
                    plot.margin = margin(0, 2, 2, 2, unit = "mm"),
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                    legend.position = c(0.98, 0.97),
                    legend.justification = c(1, 1),
                    legend.background = element_rect(fill = alpha("white", 0.85),
                                                     color = NA),
                    legend.key.size = unit(0.35, "cm"),
                    plot.subtitle = element_text(size = 8, color = "darkorange",
                                                 hjust = 0, family = "Arial"))

suppressMessages({library(magick); library(shadowtext)})
# ---- Panel B: 3D domain organization (PDB 4DQY; vector labels over raster) ----
struct_raster <- as.raster(magick::image_read("results/figures/Fig1_struct_4DQY_crop.png"))
struct_labels <- data.frame(
  x     = c(56.9, 48.9, 41.9, 29.7, 25.5, 66.1),
  y     = c(29.15, 40.30, 18.43, 24.74, 30.31, 20.76),  # cropped-image coords, y-up
  label = c("Zn1", "Zn3", "WGR", "HD", "ART", "DNA break"),
  stringsAsFactors = FALSE
)
p_struct <- ggplot() +
  annotation_raster(struct_raster, xmin = 0, xmax = 100, ymin = 0, ymax = 61.24) +
  geom_shadowtext(data = struct_labels, aes(x = x, y = y, label = label),
                  size = 3.0, fontface = "bold", family = "Arial",
                  color = "grey10", bg.color = "white", bg.r = 0.13) +
  annotate("text", x = 5.9, y = 59.5, label = "ZnF2 and BRCT not resolved in this crystal",
           size = 3.0, color = "grey40", hjust = 0, family = "Arial") +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 61.24), expand = FALSE) +
  labs(title = "B  Domain Organization in 3D (PDB 4DQY)") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 9, face = "bold", family = "Arial",
                                   margin = margin(t = 0, b = 4, l = 15)),
        plot.margin = margin(0, 2, 2, 2))

# ---- Assemble Fig 1 (2x2) ----
p_b <- p_b + labs(title = sub("^B", "C", p_b$labels$title))
p_c <- p_c + labs(title = sub("^C", "D", p_c$labels$title))
fig1 <- ((p_a | p_struct) / (p_b | p_c)) +
  plot_layout(heights = c(1, 1.1), widths = c(1, 1))

cairo_pdf("results/figures/Fig1_System_Architecture.pdf", width = 157/25.4, height = 148.7/25.4, pointsize = 8)
print(fig1)
dev.off()

message("Saved: Fig1_System_Architecture.pdf")
