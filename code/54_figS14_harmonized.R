#!/usr/bin/env Rscript
# 54_figS14_harmonized.R — SI Figure S14: harmonized-convention S2 CV1 span vs
# trapping rank (50-ps cadence, common 14-ns production window; Table S18).
#
# Purpose:  Harmonized-convention association figure (neutral vs fully protonated).
# Inputs:   data/01_curated/harmonized_panel_values.csv
# Outputs:  figures/pdf/Fig_SI_Harmonized.pdf / .png
# Depends:  ggplot2, ggrepel
# Run:      Rscript code/54_figS14_harmonized.R   (from the repository root)
suppressMessages({library(ggplot2); library(ggrepel)})

d <- read.csv("data/01_curated/harmonized_panel_values.csv", stringsAsFactors = FALSE)
stopifnot(nrow(d) == 5)

pal <- c(neutral = "#3D6BA8", fully_protonated = "#C23531")
d$state_label <- factor(c("neutral")[1], levels = c("neutral"))

long <- rbind(
  data.frame(inhibitor = d$inhibitor, rank = d$trapping_rank,
             span = d$s2cv1_harmonized_neutral, panel = "neutral"),
  data.frame(inhibitor = d$inhibitor, rank = d$trapping_rank,
             span = d$s2cv1_harmonized_fully_protonated, panel = "fully protonated")
)
long$panel <- factor(long$panel, levels = c("neutral", "fully protonated"))

lab <- data.frame(
  x = c(4.95, 4.95), y = c(64, 60),
  colour_key = c("neutral", "fully protonated"),
  label = c("neutral: \u03c1 = -0.15 (p = 0.833)",
            "protonated: \u03c1 = -0.67 (p = 0.267)")
)

# connect the three state-corrected ligands between panels
segs <- merge(long[long$panel == "neutral", c("inhibitor", "rank", "span")],
              long[long$panel == "fully protonated", c("inhibitor", "rank", "span")],
              by = c("inhibitor", "rank"), suffixes = c("_n", "_p"))
segs <- segs[segs$span_n != segs$span_p, ]

p <- ggplot(long, aes(x = rank, y = span)) +
  geom_segment(data = segs, aes(x = rank, xend = rank, y = span_n, yend = span_p),
               colour = "grey55", linewidth = 0.35, linetype = "dashed") +
  geom_point(aes(colour = panel, shape = panel), size = 2.6, stroke = 0.9) +
  geom_text_repel(aes(label = inhibitor, colour = panel), size = 3.0,
                  family = "Arial", show.legend = FALSE, max.overlaps = 20,
                  nudge_x = 0.10, segment.size = 0.2, min.segment.length = 0.6) +
  geom_text(data = lab, aes(x = x, y = y, label = label, colour = colour_key),
            hjust = 1, size = 3.0, family = "Arial", show.legend = FALSE) +
  scale_colour_manual(values = c("neutral" = pal[["neutral"]],
                                 "fully protonated" = pal[["fully_protonated"]]),
                      name = NULL) +
  scale_shape_manual(values = c("neutral" = 16, "fully protonated" = 17), name = NULL) +
  scale_x_continuous(limits = c(0.6, 5.4), breaks = c(1, 2.5, 3, 4, 5),
                     labels = c("1\n(weakest)", "2.5", "3", "4", "5\n(most potent)")) +
  labs(x = "Trapping rank", y = "S2 CV1 span, C3 (kcal/mol)",
       title = "S14  Harmonized convention (50 ps, 14-ns window)",
       subtitle = "olaparib and rucaparib tied at midrank 2.5; dashed links mark state corrections") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        plot.margin = margin(4, 6, 4, 4, "pt"))

cairo_pdf("figures/pdf/Fig_SI_Harmonized.pdf", width = 118/25.4, height = 78/25.4, pointsize = 8)
print(p)
invisible(dev.off())
ggsave("figures/pdf/Fig_SI_Harmonized.png", p, width = 118/25.4, height = 78/25.4,
       dpi = 220, type = "cairo", bg = "white")
cat("WROTE Fig_SI_Harmonized.pdf/png\n")
