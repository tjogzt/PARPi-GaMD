#!/usr/bin/env Rscript
# 04-rmsd_rmsf_figure.R — RMSD/RMSF comparison (APO vs talazoparib)
# Domain-level decomposition: All CA, HD, ART, DNA
# SCI figure: cairo_pdf, 7pt base_size, direct output dimensions

library(ggplot2)
library(data.table)
library(patchwork)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- load data ---
domains <- c("all", "HD", "ART", "DNA")
systems <- c("sys2_APO", "sys2_talazoparib")

read_rmsd <- function(sys, dom) {
  f <- file.path("results/analysis", paste0(sys, "_rmsd_", dom, ".dat"))
  if (!file.exists(f)) return(NULL)
  x <- fread(f, skip = 1, col.names = c("frame", "rmsd"))
  x[, system := fifelse(grepl("APO", sys), "APO", "talazoparib")]
  x[, domain := dom]
  x[, time_ns := frame * 0.002]  # 2fs timestep * 100 stride = 0.2 ps per frame
  x
}

read_rmsf <- function(sys, dom) {
  f <- file.path("results/analysis", paste0(sys, "_rmsf_", dom, ".dat"))
  if (!file.exists(f)) return(NULL)
  x <- fread(f, skip = 1, col.names = c("resid", "rmsf"))
  x[, system := fifelse(grepl("APO", sys), "APO", "talazoparib")]
  x[, domain := dom]
  x
}

# --- Panel A: RMSD over time ---
rmsd_list <- rbindlist(lapply(systems, function(s) {
  rbindlist(lapply(domains, function(d) read_rmsd(s, d)))
}))

rmsd_list[, domain := factor(domain, levels = domains,
                               labels = c("All Cα", "HD", "ART", "DNA"))]
rmsd_list[, system := factor(system, levels = c("APO", "talazoparib"))]

# Convert to ns for time axis
rmsd_list[, time_ns := frame * 0.2]  # stride 100 * 2fs = 0.2 ps, /1000 = ns

p_rmsd <- ggplot(rmsd_list, aes(x = time_ns, y = rmsd,
                                 color = system, group = system)) +
  geom_line(linewidth = 0.3, alpha = 0.9) +
  facet_wrap(~domain, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("APO" = "#2166AC", "talazoparib" = "#B2182B")) +
  labs(x = "Time (ns)", y = "RMSD (Å)", color = NULL,
       title = "A  RMSD Time Series") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    strip.background = element_rect(fill = "grey95"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    plot.margin = margin(3, 3, 3, 3)
  )

# --- Panel B: RMSF by residue ---
rmsf_list <- rbindlist(lapply(systems, function(s) {
  rbindlist(lapply(domains, function(d) read_rmsf(s, d)))
}))

rmsf_list[, domain := factor(domain, levels = domains,
                               labels = c("All Cα", "HD", "ART", "DNA"))]
rmsf_list[, system := factor(system, levels = c("APO", "talazoparib"))]

p_rmsf <- ggplot(rmsf_list, aes(x = resid, y = rmsf,
                                 color = system, group = system)) +
  geom_line(linewidth = 0.3, alpha = 0.9) +
  facet_wrap(~domain, scales = "free", ncol = 2) +
  scale_color_manual(values = c("APO" = "#2166AC", "talazoparib" = "#B2182B")) +
  labs(x = "Residue", y = "RMSF (Å)", color = NULL,
       title = "B  RMSF by Residue") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    strip.background = element_rect(fill = "grey95"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    plot.margin = margin(3, 3, 3, 3)
  )

# --- Panel C: Domain RMSF summary (bar plot) ---
rmsf_sum <- rmsf_list[, .(
  mean_rmsf = mean(rmsf),
  sd_rmsf = sd(rmsf)
), by = .(system, domain)]

p_rmsf_bar <- ggplot(rmsf_sum, aes(x = domain, y = mean_rmsf, fill = system)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_rmsf - sd_rmsf, ymax = mean_rmsf + sd_rmsf),
    position = position_dodge(width = 0.8), width = 0.15, linewidth = 0.3
  ) +
  scale_fill_manual(values = c("APO" = "#2166AC", "talazoparib" = "#B2182B")) +
  labs(x = NULL, y = "Mean RMSF (Å)", fill = NULL,
       title = "C  Domain RMSF Summary") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    plot.margin = margin(3, 3, 3, 3)
  )

# --- Panel D: RMSF difference (talazoparib - APO) ---
rmsf_diff <- rmsf_list[, .(
  rmsf = mean(rmsf)
), by = .(system, domain, resid)]
rmsf_diff <- dcast(rmsf_diff, domain + resid ~ system, value.var = "rmsf")
rmsf_diff[, diff := talazoparib - APO]
rmsf_diff[, domain := factor(domain, levels = c("All Cα", "HD", "ART", "DNA"))]

p_rmsf_diff <- ggplot(rmsf_diff, aes(x = resid, y = diff, fill = diff > 0)) +
  geom_col(linewidth = 0) +
  facet_wrap(~domain, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                    labels = c("TRUE" = "Increased", "FALSE" = "Decreased")) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  labs(x = "Residue", y = "ΔRMSF (TLZ − APO) (Å)", fill = NULL,
       title = "D  RMSF Change upon Talazoparib Binding") +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    strip.background = element_rect(fill = "grey95"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    plot.margin = margin(3, 3, 3, 3)
  )

# --- Assemble ---
combined <- (p_rmsd | p_rmsf) / (p_rmsf_bar | p_rmsf_diff) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Output: 190mm × 200mm (fitting SCI double column)
cairo_pdf("results/figures/04-rmsd_rmsf_comparison.pdf",
          width = 190 / 25.4, height = 200 / 25.4)
print(combined)
dev.off()

cat("Saved: results/figures/04-rmsd_rmsf_comparison.pdf\n")
