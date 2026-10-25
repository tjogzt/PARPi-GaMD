#!/usr/bin/env Rscript
# 58-trajectory_qc_figures.R — SI figures S15/S16: trajectory QC and Zn coordination
# =============================================================================
# Inputs:  PARPi-GaMD/data/analysis/s68_cv_timeseries.csv
#          PARPi-GaMD/data/analysis/s68_trajectory_qc.csv
# Outputs: results/figures/Fig_SI_Trajectory_QC.pdf/.png
#          results/figures/Fig_SI_Zn_Coordination.pdf/.png
# Run:     R --file=code/58-trajectory_qc_figures.R --no-save --no-restore
# =============================================================================
set.seed(49)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

PROJ <- "."
OUT  <- "figures/pdf"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cv <- read.csv(file.path(PROJ, "data/analysis/s68_cv_timeseries.csv"))
qc <- read.csv(file.path(PROJ, "data/analysis/s68_trajectory_qc.csv"))

lev <- c("talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib",
         "niraparib_prot", "rucaparib_prot", "veliparib_prot")
lab <- c(talazoparib = "Talazoparib", niraparib = "Niraparib", olaparib = "Olaparib",
         rucaparib = "Rucaparib", veliparib = "Veliparib",
         niraparib_prot = "Niraparib (+1)", rucaparib_prot = "Rucaparib (+1)",
         veliparib_prot = "Veliparib (+1)")
cols <- c(talazoparib = "#FF7F00", niraparib = "#4DAF4A", olaparib = "#E41A1C",
          rucaparib = "#984EA3", veliparib = "#377EB8",
          niraparib_prot = "#4DAF4A", rucaparib_prot = "#984EA3", veliparib_prot = "#377EB8")
for (d in list(cv, qc)) {
  d$system <- factor(d$system, levels = lev)
  d$state <- factor(ifelse(grepl("_prot$", as.character(d$system)), "protonated", "neutral"),
                    levels = c("neutral", "protonated"))
}
cv$system <- factor(cv$system, levels = lev)
qc$system <- factor(qc$system, levels = lev)
cv$state <- factor(ifelse(grepl("_prot$", as.character(cv$system)), "protonated", "neutral"),
                   levels = c("neutral", "protonated"))
qc$state <- factor(ifelse(grepl("_prot$", as.character(qc$system)), "protonated", "neutral"),
                   levels = c("neutral", "protonated"))

smooth_vec <- function(x, k = 21) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
cv$cv1_s <- NA_real_; cv$cv2_s <- NA_real_
for (s in levels(cv$system)) {
  i <- which(cv$system == s)
  cv$cv1_s[i] <- smooth_vec(cv$cv1[i])
  cv$cv2_s[i] <- smooth_vec(cv$cv2[i])
}
qc$rmsd_s <- NA_real_
for (s in levels(qc$system)) {
  i <- which(qc$system == s)
  qc$rmsd_s[i] <- smooth_vec(qc$rmsd_ca[i])
}

base_theme <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "grey20"),
    axis.ticks = element_line(linewidth = 0.35, colour = "grey20"),
    axis.text = element_text(size = 7.2, colour = "grey20"),
    axis.title = element_text(size = 8, colour = "grey10"),
    plot.tag = element_text(size = 8, face = "bold", family = "Arial"),
    plot.tag.position = c(0.01, 0.99),
    legend.text = element_text(size = 7.2),
    legend.title = element_blank(),
    legend.key.width = unit(0.9, "lines"),
    legend.key.height = unit(0.8, "lines"),
    panel.grid = element_blank()
  )

band <- annotate("rect", xmin = 0, xmax = 14, ymin = -Inf, ymax = Inf,
                 fill = "grey90", alpha = 0.55, colour = NA)
vline <- geom_vline(xintercept = 14, linetype = 2, linewidth = 0.35, colour = "grey45")

cv1 <- ggplot(cv, aes(time_ns, cv1, colour = system, linetype = state)) +
  band + vline +
  geom_line(aes(y = cv1), linewidth = 0.18, alpha = 0.16) +
  geom_line(aes(y = cv1_s), linewidth = 0.55, alpha = 0.95) +
  scale_colour_manual(values = cols, labels = lab) +
  scale_linetype_manual(values = c(neutral = "solid", protonated = "22")) +
  scale_x_continuous(limits = c(0, 26), breaks = seq(0, 25, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = "time (ns)", y = "CV1 protein–DNA COM (Å)") +
  base_theme + theme(legend.position = "none") +
  annotate("text", x = 0.6, y = Inf, label = "A", vjust = 1.3, hjust = 0,
           size = 3.0, fontface = "bold", family = "Arial")

cv2 <- ggplot(cv, aes(time_ns, cv2, colour = system, linetype = state)) +
  band + vline +
  geom_line(aes(y = cv2), linewidth = 0.18, alpha = 0.16) +
  geom_line(aes(y = cv2_s), linewidth = 0.55, alpha = 0.95) +
  scale_colour_manual(values = cols, labels = lab) +
  scale_linetype_manual(values = c(neutral = "solid", protonated = "22")) +
  scale_x_continuous(limits = c(0, 26), breaks = seq(0, 25, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = "time (ns)", y = "CV2 HD–ART COM (Å)") +
  base_theme + theme(legend.position = "none") +
  annotate("text", x = 0.6, y = Inf, label = "B", vjust = 1.3, hjust = 0,
           size = 3.0, fontface = "bold", family = "Arial")

rmsd <- ggplot(qc, aes(time_ns, rmsd_ca, colour = system, linetype = state)) +
  band + vline +
  geom_line(aes(y = rmsd_ca), linewidth = 0.18, alpha = 0.16) +
  geom_line(aes(y = rmsd_s), linewidth = 0.55, alpha = 0.95) +
  scale_colour_manual(values = cols, labels = lab) +
  scale_linetype_manual(values = c(neutral = "solid", protonated = "22")) +
  scale_x_continuous(limits = c(0, 26), breaks = seq(0, 25, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = "time (ns)", y = "backbone RMSD (Å)") +
  base_theme + theme(legend.position = "none") +
  annotate("text", x = 0.6, y = Inf, label = "C", vjust = 1.3, hjust = 0,
           size = 3.0, fontface = "bold", family = "Arial")

fig1 <- (cv1 | cv2 | rmsd) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(0.1, "cm"))
ggsave(file.path(OUT, "Fig_SI_Trajectory_QC.pdf"), fig1,
       width = 178 / 25.4, height = 82 / 25.4, device = cairo_pdf, family = "Arial")
ggsave(file.path(OUT, "Fig_SI_Trajectory_QC.png"), fig1,
       width = 178 / 25.4, height = 82 / 25.4, dpi = 600)

# --- Zn coordination: minimum Zn–Sγ distance per structural Zn -----------------
zn <- rbind(
  data.frame(system = qc$system, state = qc$state, time_ns = qc$time_ns,
             site = "Zn1", dist = qc$zn1_min),
  data.frame(system = qc$system, state = qc$state, time_ns = qc$time_ns,
             site = "Zn2", dist = qc$zn2_min),
  data.frame(system = qc$system, state = qc$state, time_ns = qc$time_ns,
             site = "Zn3", dist = qc$zn3_min)
)
zn$site <- factor(zn$site, levels = c("Zn1", "Zn2", "Zn3"))

zn_panel <- function(site_name, tag) {
  d <- zn[zn$site == site_name, ]
  ggplot(d, aes(time_ns, dist, colour = system, linetype = state)) +
    geom_hline(yintercept = 2.8, linetype = 3, linewidth = 0.35, colour = "grey45") +
    geom_line(linewidth = 0.38, alpha = 0.9) +
    scale_colour_manual(values = cols, labels = lab) +
    scale_linetype_manual(values = c(neutral = "solid", protonated = "22")) +
    scale_x_continuous(limits = c(0, 26), breaks = seq(0, 25, 5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(1.75, 3.0), breaks = seq(1.8, 3.0, 0.2),
                       expand = c(0, 0)) +
    labs(x = "time (ns)", y = "min Zn–Sγ (Å)") +
    base_theme + theme(legend.position = "none") +
    annotate("text", x = 0.6, y = 2.97, label = tag, vjust = 1, hjust = 0,
             size = 3.0, fontface = "bold", family = "Arial") +
    annotate("text", x = 25.6, y = 2.80, label = "2.8 Å", hjust = 1, vjust = -0.25,
             size = 2.4, colour = "grey35", family = "Arial")
}

fig2 <- (zn_panel("Zn1", "A") | zn_panel("Zn2", "B") | zn_panel("Zn3", "C")) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(0.1, "cm"))
ggsave(file.path(OUT, "Fig_SI_Zn_Coordination.pdf"), fig2,
       width = 178 / 25.4, height = 70 / 25.4, device = cairo_pdf, family = "Arial")
ggsave(file.path(OUT, "Fig_SI_Zn_Coordination.png"), fig2,
       width = 178 / 25.4, height = 70 / 25.4, dpi = 600)

cat("WROTE", file.path(OUT, "Fig_SI_Trajectory_QC.pdf"), "and",
    file.path(OUT, "Fig_SI_Zn_Coordination.pdf"), "\n")
