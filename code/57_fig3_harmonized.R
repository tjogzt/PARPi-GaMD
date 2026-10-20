#!/usr/bin/env Rscript
# 57_fig3_harmonized.R — main-text Figure 3 under the harmonized convention
# (50-ps cadence, common 14-ns production window). Replaces the archival
# three-panel figure as the primary display; the archival version is retained
# as SI Figure S14.
#
# Panels
#   A  S2 CV1 span (C3) vs trapping rank for the four chemical-state panels
#   B  span ratio (S1 / S2 CV2) vs trapping rank for the same four panels
#   C  estimator-consistency matrix (4 metrics x {3 span estimators or 1} x 4 panels)
#
# Inputs: data/analysis/r6/s63_harmonized_spans.csv
#         data/analysis/r6/s63_harmonized_matrix.csv
# Outputs: figures/pdf/Fig3_Harmonized.pdf / .png
suppressMessages({library(ggplot2); library(patchwork); library(ggrepel)})

sp <- read.csv("data/analysis/r6/s63_harmonized_spans.csv", stringsAsFactors = FALSE)
mx <- read.csv("data/analysis/r6/s63_harmonized_matrix.csv", stringsAsFactors = FALSE)

S1 <- c(tala = 42.2, nira = 50.8, nira_prot = 59.2, ola = 55.4,
        ruca = 27.6, ruca_prot = 27.7, veli = 40.3, veli_prot = 23.7)
Panels <- list(
  neutral = c(tala = "tala", nira = "nira", ola = "ola", ruca = "ruca", veli = "veli"),
  `niraparib +1` = c(tala = "tala", nira = "nira_prot", ola = "ola", ruca = "ruca", veli = "veli"),
  `rucaparib +1` = c(tala = "tala", nira = "nira", ola = "ola", ruca = "ruca_prot", veli = "veli"),
  `all three +1` = c(tala = "tala", nira = "nira_prot", ola = "ola", ruca = "ruca_prot", veli = "veli_prot")
)
pnames <- names(Panels)
pcol <- c("#3D6BA8", "#177CB0", "#9D2933", "#C23531")
names(pcol) <- pnames
off <- c(-0.27, -0.09, 0.09, 0.27)
names(off) <- pnames
Y <- c(tala = 5, nira = 4, ola = 2.5, ruca = 2.5, veli = 1)
lab <- c(tala = "tala", nira = "nira", ola = "ola", ruca = "ruca", veli = "veli")

getspan <- function(sys, metric, est) {
  sp$span_kcal[sp$system == sys & sp$metric == metric & sp$estimator == est]
}

rows <- list()
for (p in pnames) for (k in names(Panels[[p]])) {
  sys <- Panels[[p]][[k]]
  rows[[length(rows) + 1]] <- data.frame(
    panel = p, inhibitor = k, rank = Y[[k]] + off[[p]],
    cv1 = getspan(sys, "CV1", "mm"),
    ratio = S1[[sys]] / getspan(sys, "CV2", "mm"))
}
d <- do.call(rbind, rows)
d$panel <- factor(d$panel, levels = pnames)

recode_panel <- function(x) {
  out <- as.character(x)
  out[out == "niraparib-protonated"] <- "niraparib +1"
  out[out == "rucaparib-protonated"] <- "rucaparib +1"
  out[out == "fully protonated"] <- "all three +1"
  factor(out, levels = pnames)
}

ann <- function(metric, ys) {
  a <- mx[mx$metric == metric & mx$estimator == "mm", c("panel", "rho", "p_exact")]
  a$panel <- recode_panel(a$panel)
  a$label <- sprintf("%s  \u03c1 = %+.2f (%.3f)", a$panel, a$rho, a$p_exact)
  a$y <- ys
  a$x <- 0.62
  a
}
aA <- ann("S2CV1", c(22.4, 20.9, 19.5, 18.1))
aB <- ann("ratio", c(0.60, 0.70, 0.81, 0.94))

pa <- ggplot(d, aes(rank, cv1, colour = panel)) +
  geom_point(size = 2.1, stroke = 0.9) +
  scale_colour_manual(values = pcol, name = NULL) +
  scale_x_continuous(limits = c(0.55, 5.45), breaks = c(1, 2.5, 3, 4, 5),
                     labels = c("1\n(weakest)", "2.5", "3", "4", "5\n(most potent)")) +
  scale_y_continuous(limits = c(17, 59.5)) +
  labs(x = "Trapping rank", y = "S2 CV1 span, C3 (kcal/mol)",
       title = "A  CV1 span vs trapping rank (harmonized)",
       subtitle = "\u03c1 (p): neutral \u22120.15 (.833) \u00b7 nira+1 +0.36 (.633)\nruca+1 \u22120.15 (.833) \u00b7 all+1 \u22120.67 (.267)") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.text = element_text(size = 7),
        legend.key.size = unit(0.26, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.2, colour = "grey25", lineheight = 1.05),
        plot.margin = margin(4, 6, 4, 4, "pt"))

pb <- ggplot(d, aes(rank, ratio, colour = panel)) +
  geom_point(size = 2.1, stroke = 0.9) +
  scale_colour_manual(values = pcol, name = NULL) +
  scale_y_log10(limits = c(0.48, 4.8), breaks = c(0.5, 1, 2, 4)) +
  scale_x_continuous(limits = c(0.55, 5.45), breaks = c(1, 2.5, 3, 4, 5),
                     labels = c("1\n(weakest)", "2.5", "3", "4", "5\n(most potent)")) +
  labs(x = "Trapping rank", y = "Span ratio, S1 / S2 CV2 (log scale)",
       title = "B  Span ratio vs trapping rank (harmonized)",
       subtitle = "\u03c1 (p): neutral \u22120.41 (.500) \u00b7 nira+1 +0.05 (1.000)\nruca+1 \u22120.41 (.500) \u00b7 all+1 +0.05 (1.000)") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.text = element_text(size = 7),
        legend.key.size = unit(0.26, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.2, colour = "grey25", lineheight = 1.05),
        plot.margin = margin(4, 6, 4, 4, "pt"))

# panel C heatmap: collapse S1 (estimator-independent) to one row
mk <- mx
mk$metric2 <- factor(mk$metric, levels = c("S2CV1", "S2CV2", "S1", "ratio"))
mk <- mk[!(mk$metric == "S1" & mk$estimator != "mm"), ]
mk$rowlab <- ifelse(mk$metric == "S1", "S1",
                    paste0(mk$metric, " \u00b7 ", mk$estimator))
mk$row <- factor(mk$rowlab,
                 levels = rev(c(paste0("S2CV1 \u00b7 ", c("mm", "q90", "t95")),
                                paste0("S2CV2 \u00b7 ", c("mm", "q90", "t95")),
                                "S1",
                                paste0("ratio \u00b7 ", c("mm", "q90", "t95")))))
mk$panel <- recode_panel(mk$panel)

pc <- ggplot(mk, aes(panel, row, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.2f", rho),
                colour = abs(rho) > 0.55), size = 2.6, family = "Arial") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
  scale_fill_gradient2(low = "#3D6BA8", mid = "#F7F7F7", high = "#C23531",
                       midpoint = 0, limits = c(-1, 1), name = expression(rho)) +
  labs(x = NULL, y = NULL, title = "C  Estimator--state matrix (harmonized)") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 7.5),
        axis.text.y = element_text(size = 7.2),
        plot.title = element_text(size = 9, face = "bold"),
        legend.key.height = unit(0.9, "cm"), legend.key.width = unit(0.28, "cm"),
        legend.text = element_text(size = 7),
        plot.margin = margin(4, 6, 4, 4, "pt"))

fig <- (pa | pb) / pc + plot_layout(heights = c(1.08, 1.30))

cairo_pdf("figures/pdf/Fig3_Harmonized.pdf", width = 172/25.4, height = 168/25.4,
          pointsize = 8)
print(fig)
invisible(dev.off())
ggsave("figures/pdf/Fig3_Harmonized.png", fig, width = 172/25.4, height = 168/25.4,
       dpi = 220, type = "cairo", bg = "white")
cat("WROTE Fig3_Harmonized.pdf/png\n")
