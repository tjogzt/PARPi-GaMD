#!/usr/bin/env Rscript
# 17-toc_graphic.R — TOC graphic: "Same Pocket, Two Fates"
library(ggplot2)
library(dplyr)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ---- TOC: Conceptual schematic ---------------------------------------------
# Left: PARP1+DNA with inhibitor in CAT pocket
# Center: Arrow showing allosteric signal transduction
# Right: Two outcomes — DNA retention (Type II) vs DNA release (Type III)

# Data for the conceptual 2D mechanism map
toc_data <- data.frame(
  x = c(1, 1, 2, 2, 3, 3, 3, 4, 4),
  y = c(3, 1, 3, 1, 3.5, 2, 0.5, 3, 1),
  label = c("PARP1\n+DNA+inhibitor", "CAT pocket\nbinding", 
            "Allosteric\nsignal", "HD domain\nresponse",
            "", "Type II\nNeutral", "Type III\nPro-release",
            "DNA RETENTION\nHigh trapping", "DNA RELEASE\nLow trapping"),
  color = c("grey30", "grey30", "#FFD92F", "#FFD92F",
            "white", "#E41A1C", "#377EB8", "#E41A1C", "#377EB8"),
  size = c(3, 2.5, 3, 2.5, 0, 3.5, 3.5, 4, 4),
  stringsAsFactors = FALSE
)

# Arrows
arrows <- data.frame(
  x = c(1.5, 1.5, 2.5, 2.5),
  xend = c(2.4, 2.4, 3.4, 3.4),
  y = c(2.2, 1.8, 2.8, 1.2),
  yend = c(2.2, 1.8, 2.8, 1.2),
  stringsAsFactors = FALSE
)

p <- ggplot() +
  # Background boxes
  annotate("rect", xmin = 0.3, xmax = 1.7, ymin = 0.3, ymax = 3.7,
           fill = "grey95", color = "grey80", linewidth = 0.3) +
  annotate("rect", xmin = 1.8, xmax = 2.8, ymin = 0.3, ymax = 3.7,
           fill = "#FFF8E1", color = "#FFD92F", linewidth = 0.3) +
  annotate("rect", xmin = 2.9, xmax = 4.7, ymin = 0.3, ymax = 3.7,
           fill = "white", color = "grey80", linewidth = 0.3) +
  # PARP1 schematic (left)
  annotate("text", x = 1, y = 3.2, label = "PARP1 + DNA break", size = 3.5, fontface = "bold") +
  annotate("text", x = 1, y = 2.8, label = "+ inhibitor", size = 3, color = "grey40") +
  # CAT pocket (oval)
  annotate("point", x = 1, y = 2, shape = 21, size = 20, fill = "#E41A1C", alpha = 0.15, color = "#E41A1C") +
  annotate("text", x = 1, y = 2, label = "CAT\npocket", size = 2.5, color = "#E41A1C", fontface = "bold") +
  # Inhibitor in pocket
  annotate("text", x = 1, y = 1.3, label = "Inhibitor", size = 2.5, color = "grey40") +
  # Arrow
  annotate("segment", x = 1.5, xend = 2.3, y = 2.2, yend = 2.2,
           arrow = arrow(type = "closed", length = unit(0.1, "inches")),
           color = "#FFD92F", linewidth = 1) +
  # HD domain (center)
  annotate("text", x = 2.3, y = 3.2, label = "HD Domain", size = 3.5, fontface = "bold", color = "#FF7F00") +
  annotate("text", x = 2.3, y = 2.8, label = "Allosteric\nResponse", size = 3, color = "grey40", lineheight = 0.85) +
  # Branching arrows
  annotate("segment", x = 2.5, xend = 3.4, y = 3, yend = 3,
           arrow = arrow(type = "closed", length = unit(0.1, "inches")),
           color = "#E41A1C", linewidth = 1) +
  annotate("segment", x = 2.5, xend = 3.4, y = 1.2, yend = 1.2,
           arrow = arrow(type = "closed", length = unit(0.1, "inches")),
           color = "#377EB8", linewidth = 1) +
  # Type II outcome (upper right)
  annotate("text", x = 4, y = 3.2, label = "Neutral Allostery", size = 3.5, fontface = "bold", color = "#E41A1C") +
  annotate("text", x = 4, y = 2.8, label = "HD stays closed\nDNA retained\n100× trapping\n(e.g., Talazoparib)", 
           size = 2.5, color = "grey30", lineheight = 0.9) +
  # Type III outcome (lower right)
  annotate("text", x = 4, y = 1.5, label = "Pro-Release Allostery", size = 3.5, fontface = "bold", color = "#377EB8") +
  annotate("text", x = 4, y = 1.1, label = "HD opens\nDNA released\n0.02× trapping\n(e.g., Veliparib)", 
           size = 2.5, color = "grey30", lineheight = 0.9) +
  # Divider
  annotate("segment", x = 3, xend = 3, y = 0.5, yend = 3.5,
           linetype = "dashed", color = "grey70", linewidth = 0.3) +
  # Title
  annotate("text", x = 2.5, y = 3.9, label = "Same Pocket, Two Fates",
           size = 5, fontface = "bold", color = "grey20") +
  xlim(0, 5) + ylim(0, 4.2) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

cairo_pdf("results/figures/TOC_Graphic.pdf", width = 120/25.4, height = 54/25.4, pointsize = 7)
print(p)
dev.off()

cat("Saved: TOC_Graphic.pdf\n")

# ---- Alternative: minimalist version for square format ----
p2 <- ggplot() +
  annotate("rect", xmin = 0, xmax = 3, ymin = 0, ymax = 3, fill = "white", color = NA) +
  # Center circle = CAT pocket
  annotate("point", x = 1, y = 1.5, shape = 21, size = 15, fill = "#E41A1C", alpha = 0.1, color = "#E41A1C", stroke = 0.5) +
  annotate("text", x = 1, y = 1.5, label = "CAT\nPocket", size = 2.5, color = "#E41A1C", fontface = "bold") +
  # Left: DNA complex
  annotate("text", x = 0.2, y = 2.7, label = "PARP1+DNA\n+Inhibitor", size = 2.5, fontface = "bold") +
  # Arrow to branch
  annotate("segment", x = 1.3, xend = 1.8, y = 1.5, yend = 2.2,
           arrow = arrow(type = "closed", length = unit(0.06, "inches")),
           color = "#E41A1C", linewidth = 0.8) +
  annotate("segment", x = 1.3, xend = 1.8, y = 1.5, yend = 0.8,
           arrow = arrow(type = "closed", length = unit(0.06, "inches")),
           color = "#377EB8", linewidth = 0.8) +
  # Outcomes
  annotate("text", x = 2.5, y = 2.3, label = "Retention\n(High Trap)", size = 2.5, color = "#E41A1C", fontface = "bold") +
  annotate("text", x = 2.5, y = 0.7, label = "Release\n(Low Trap)", size = 2.5, color = "#377EB8", fontface = "bold") +
  # HD label
  annotate("text", x = 1.5, y = 1.5, label = "HD Allostery", size = 2, color = "#FF7F00", angle = 90) +
  xlim(-0.2, 3.5) + ylim(0, 3) +
  theme_void()

cairo_pdf("results/figures/TOC_Graphic_Minimal.pdf", width = 80/25.4, height = 80/25.4, pointsize = 7)
print(p2)
dev.off()

cat("Saved: TOC_Graphic_Minimal.pdf\n")
message("===== TOC complete =====")
