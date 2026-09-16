#!/usr/bin/env Rscript
# Fig 1: PARP1 domain architecture + system schematics + inhibitor panel
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

theme_7pt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 6))

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
  geom_text(data = domain_data, aes(x = (start + end) / 2, y = y, label = domain), size = 2.5, fontface = "bold") +
  geom_segment(aes(x = 50, xend = 50, y = 1.25, yend = 1.55), linewidth = 0.2, color = "grey50") +
  geom_segment(aes(x = 300, xend = 300, y = 1.25, yend = 1.55), linewidth = 0.2, color = "grey50") +
  geom_segment(aes(x = 590, xend = 590, y = 1.25, yend = 1.55), linewidth = 0.2, color = "grey50") +
  geom_segment(aes(x = 725, xend = 725, y = 1.25, yend = 1.55), linewidth = 0.2, color = "grey50") +
  geom_segment(aes(x = 870, xend = 870, y = 1.25, yend = 1.55), linewidth = 0.2, color = "grey50") +
  annotate("text", x = 50, y = 1.7, label = "DNA\nbinding", size = 2, color = "grey30", lineheight = 0.85) +
  annotate("text", x = 300, y = 1.7, label = "DNA\nbinding", size = 2, color = "grey30", lineheight = 0.85) +
  annotate("text", x = 590, y = 1.7, label = "DNA\nbinding", size = 2, color = "grey30", lineheight = 0.85) +
  annotate("text", x = 725, y = 1.7, label = "Type I\n(EB-47)", size = 2, color = "grey30", lineheight = 0.85) +
  annotate("text", x = 870, y = 1.7, label = "CAT Pocket\nType II/III\n(this study)", size = 2, color = "grey30", lineheight = 0.85) +
  scale_fill_manual(values = setNames(domain_data$color, domain_data$domain), guide = "none") +
  scale_x_continuous(limits = c(0, 1020), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 2.0)) +
  labs(x = "Residue Number", y = NULL,
       title = "A  PARP1 Domain Architecture & Inhibitor Binding Sites") +
  theme_7pt + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                    panel.grid.major.y = element_blank())

# ---- Panel B: S1/S2 System Schematic ----
# GaMD dual-system schematic
sys_data <- data.frame(
  x = c(1, 1, 2, 2),
  y = c(2, 1, 2, 1),
  system = c("S1: CAT-only", "S2: DNA-bound", "S1: CAT-only", "S2: DNA-bound"),
  label  = c("HD+ART\n+ inhibitor\n~30K atoms\nCV: HD-ART", 
             "PARP1+DNA\n+ inhibitor\n~120K atoms\nCV1: Prot-DNA\nCV2: HD-ART",
             "", ""),
  stringsAsFactors = FALSE
)

p_b <- ggplot(sys_data[1:2,], aes(x = x, y = y)) +
  # S1 box
  annotate("rect", xmin = 0.2, xmax = 1.8, ymin = 1.5, ymax = 2.5,
           fill = "#2166AC", alpha = 0.1, color = "#2166AC", linewidth = 0.5) +
  # S2 box
  annotate("rect", xmin = 0.2, xmax = 1.8, ymin = 0.5, ymax = 1.5,
           fill = "#B2182B", alpha = 0.1, color = "#B2182B", linewidth = 0.5) +
  # S1 label
  annotate("text", x = 1, y = 2.3, label = "S1: CAT-only (HD+ART)", 
           size = 3, fontface = "bold", color = "#2166AC") +
  annotate("text", x = 1, y = 1.95, label = "~61,500 atoms | 10 systems | 19-200 ns\nCV: HD-ART COM distance\nMeasures: allosteric response",
           size = 2.2, color = "grey30", lineheight = 0.9) +
  # S2 label
  annotate("text", x = 1, y = 1.3, label = "S2: Full PARP1 + DNA break",
           size = 3, fontface = "bold", color = "#B2182B") +
  annotate("text", x = 1, y = 0.95, label = "~286,000 atoms | 10 systems | 22-26 ns\nCV1: Protein-DNA COM | CV2: HD-ART COM\nMeasures: complete trapping environment",
           size = 2.2, color = "grey30", lineheight = 0.9) +
  # Arrow between S1 and S2
  annotate("segment", x = 1, xend = 1, y = 1.5, yend = 0.5,
           arrow = arrow(type = "closed", length = unit(0.08, "inches")),
           color = "grey50", linewidth = 0.5) +
  annotate("text", x = 1.3, y = 1, label = "S1/S2 ratio =\nAllosteric\nAmplification\nIndex (AAI)",
           size = 2.2, color = "grey30", lineheight = 0.85, hjust = 0) +
  xlim(0, 3.2) + ylim(0.3, 2.7) +
  labs(title = "B  GaMD Dual-System Design") +
  theme_7pt + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                    axis.title = element_blank(), panel.grid = element_blank())

# ---- Panel C: Inhibitor Structures + Trapping Data ----
# Canonical sources: data/01_curated/trapping_potency.csv (x olaparib, PMIDs in file);
# class from common/ligands.csv. AZD5305: left-censored <0.01x (Pires 2025, PMID 40021124).
trap_csv <- fread("data/01_curated/trapping_potency.csv")
cls_csv  <- fread("common/ligands.csv")
cls_map  <- setNames(cls_csv$class, cls_csv$ligand)

inhib_names <- c("Talazoparib", "Niraparib", "Olaparib", "Rucaparib", "Veliparib", "AZD5305")
trap_vals <- setNames(trap_csv$trapping_x_olaparib, trap_csv$inhibitor)[tolower(inhib_names)]
trap_vals[is.na(trap_vals)] <- 0.01  # AZD5305 left-censored lower bound
cls_vals <- cls_map[c("talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib", "AZD5305")]

inhib_data <- data.frame(
  ligand       = factor(inhib_names, levels = inhib_names),
  trap_potency = unname(trap_vals),
  type         = gsub("_", " ", unname(cls_vals)),
  bar_label    = ifelse(inhib_names == "AZD5305", "<0.01×",
                        sprintf("%g×", unname(trap_vals))),
  stringsAsFactors = FALSE
)

p_c <- ggplot(inhib_data, aes(x = ligand, y = trap_potency, fill = type)) +
  geom_bar(stat = "identity", width = 0.7, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = bar_label, y = trap_potency * 1.5),
            size = 2.5, vjust = 0) +
  scale_fill_manual(values = c("Type II" = "#E41A1C", "Type III" = "#377EB8", "Unknown" = "darkorange"), guide = "none") +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
  labs(x = NULL, y = "Trapping Potency (× Olaparib, log scale)",
       title = "C  Inhibitor Trapping Potency & Type Classification") +
  theme_7pt + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# ---- Assemble Fig 1 ----
fig1 <- (p_a / (p_b | p_c)) +
  plot_layout(heights = c(1, 1.5))

cairo_pdf("results/figures/Fig1_System_Architecture.pdf", width = 190/25.4, height = 180/25.4, pointsize = 7)
print(fig1)
dev.off()

message("Saved: Fig1_System_Architecture.pdf")
