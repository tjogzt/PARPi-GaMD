#!/usr/bin/env Rscript
# 24-transition_barriers.R — Extract transition barriers from PMF curves
#
# For each S1 and S2 system, identify local minima and the barrier heights
# between them. Compute Kramers approximate rates.
#
# Input:  results/analysis/sys1_*_pmf_c3.xvg, sys2 pmf files
# Output: results/figures/Fig_Transition_Barriers.pdf
#         results/analysis/barrier_stats.csv

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Config ----------------------------------------------------------------
data_dir  <- "results/analysis"
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_7pt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
        legend.key.size = unit(0.3, "cm"))

# ---- Read XVG --------------------------------------------------------------
read_xvg <- function(path) {
  lines <- readLines(path)
  data_lines <- lines[!grepl("^[@#]", lines)]
  data_lines <- data_lines[data_lines != ""]
  x <- read.table(text = data_lines, header = FALSE)
  colnames(x) <- c("RC", "PMF")
  x$PMF <- x$PMF - min(x$PMF)  # zero baseline
  x
}

# ---- Find local extrema ----------------------------------------------------
find_extrema <- function(x, y, smooth_span = 0.05) {
  # Smooth first
  n <- length(y)
  span <- max(3, floor(n * smooth_span))
  y_smooth <- stats::filter(y, rep(1/span, span), sides = 2)
  y_smooth[is.na(y_smooth)] <- y[is.na(y_smooth)]
  
  # Find minima
  minima <- c()
  for (i in 2:(n-1)) {
    if (!is.na(y_smooth[i]) && !is.na(y_smooth[i-1]) && !is.na(y_smooth[i+1])) {
      if (y_smooth[i] < y_smooth[i-1] && y_smooth[i] < y_smooth[i+1]) {
        minima <- c(minima, i)
      }
    }
  }
  
  # Find maxima (barriers between minima)
  barriers <- list()
  if (length(minima) >= 2) {
    for (m in 1:(length(minima)-1)) {
      seg <- minima[m]:minima[m+1]
      max_idx <- seg[which.max(y_smooth[seg])]
      barriers[[length(barriers) + 1]] <- data.frame(
        from_min_idx = minima[m],
        to_min_idx = minima[m+1],
        from_min_RC = x[minima[m]],
        to_min_RC = x[minima[m+1]],
        barrier_idx = max_idx,
        barrier_RC = x[max_idx],
        barrier_height = y[max_idx] - min(y[minima[m]], y[minima[m+1]]),
        from_min_energy = y[minima[m]],
        to_min_energy = y[minima[m+1]],
        delta_energy = y[minima[m+1]] - y[minima[m]]
      )
    }
  }
  
  if (length(barriers) == 0) {
    return(data.frame(
      from_min_RC = NA, to_min_RC = NA, barrier_RC = NA,
      barrier_height = NA, delta_energy = NA
    ))
  }
  
  barriers_df <- bind_rows(barriers)
  
  # Also add the highest barrier (global maximum between global minimum and others)
  global_min_idx <- minima[which.min(y[minima])]
  highest_barrier <- barriers_df[which.max(barriers_df$barrier_height), ]
  
  barriers_df
}

# ---- Kramers rate (approximate) --------------------------------------------
kramers_rate <- function(barrier_height_kcal, temperature_K = 310) {
  # k = (kT/h) * exp(-ΔG‡/kT)
  # Returns rate in s^-1 (approximate, without prefactor)
  kT <- 0.001987 * temperature_K  # kcal/mol
  h <- 1.583e-37  # kcal·s (Planck constant in kcal·s)
  prefactor <- kT / h
  prefactor * exp(-barrier_height_kcal / kT)
}

# ---- Process all S1 systems (c3 cumulant) ----------------------------------
s1_systems <- c("APO", "AZD5305", "olaparib", "veliparib", "niraparib", "rucaparib", "talazoparib")
all_barriers <- list()
pmf_list <- list()

for (sys in s1_systems) {
  fname <- file.path(data_dir, paste0("sys1_", sys, "_pmf_c3.xvg"))
  if (!file.exists(fname)) {
    cat(sprintf("[SKIP] %s: no S1 C3 PMF\n", sys))
    next
  }
  
  pmf <- read_xvg(fname)
  pmf$system <- sys
  pmf$arch <- "S1"
  pmf_list[[sys]] <- pmf
  
  barriers <- find_extrema(pmf$RC, pmf$PMF)
  
  if (nrow(barriers) > 0 && !all(is.na(barriers$barrier_height))) {
    barriers$system <- sys
    barriers$architecture <- "S1"
    barriers$well_depth <- max(pmf$PMF)  # total well depth
    
    # Find max barrier
    max_barrier <- max(barriers$barrier_height, na.rm = TRUE)
    barriers$max_barrier <- max_barrier
    
    # Kramers rate for the highest barrier
    barriers$kramers_rate <- kramers_rate(max_barrier)
    barriers$log10_rate <- log10(barriers$kramers_rate)
    
    all_barriers[[length(all_barriers) + 1]] <- barriers
    
    cat(sprintf("%s S1: %d barriers, max=%.1f kcal/mol, rate=%.2e s^-1\n",
                sys, nrow(barriers), max_barrier, barriers$kramers_rate[1]))
  } else {
    cat(sprintf("%s S1: single-well (no barriers)\n", sys))
  }
}

# ---- Process all S2 systems (CV1 = HD-ART distance) ------------------------
s2_systems <- c("APO", "AZD5305", "talazoparib", "veliparib", "niraparib", "olaparib", "rucaparib")

for (sys in s2_systems) {
  # Try c3 first, then c2, then c1
  found <- FALSE
  for (cum in c("c3", "c2", "c1")) {
    fname <- file.path(data_dir, paste0("pmf-", cum, "-sys2_", sys, "_CV1_cv.dat.xvg"))
    if (file.exists(fname)) {
      pmf <- read_xvg(fname)
      pmf$system <- sys
      pmf$arch <- "S2"
      pmf_list[[paste0("S2_", sys)]] <- pmf
      
      barriers <- find_extrema(pmf$RC, pmf$PMF)
      
      if (nrow(barriers) > 0 && !all(is.na(barriers$barrier_height))) {
        barriers$system <- sys
        barriers$architecture <- "S2"
        barriers$well_depth <- max(pmf$PMF)
        max_barrier <- max(barriers$barrier_height, na.rm = TRUE)
        barriers$max_barrier <- max_barrier
        barriers$kramers_rate <- kramers_rate(max_barrier)
        barriers$log10_rate <- log10(barriers$kramers_rate)
        all_barriers[[length(all_barriers) + 1]] <- barriers
        cat(sprintf("%s S2 (CV1, %s): %d barriers, max=%.1f kcal/mol\n",
                    sys, cum, nrow(barriers), max_barrier))
      } else {
        cat(sprintf("%s S2 (CV1, %s): single-well\n", sys, cum))
      }
      found <- TRUE
      break
    }
  }
  if (!found) cat(sprintf("[SKIP] %s: no S2 PMF\n", sys))
}

# ---- Aggregate results -----------------------------------------------------
barriers_all <- bind_rows(all_barriers)
cat(sprintf("\nTotal barriers identified: %d\n", nrow(barriers_all)))

# Per-system summary
barrier_summary <- barriers_all %>% 
  group_by(system, architecture) %>%
  summarise(
    well_depth = unique(well_depth)[1],
    n_barriers = n(),
    max_barrier = max(barrier_height, na.rm = TRUE),
    mean_barrier = mean(barrier_height, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    max_barrier_frac = max_barrier / well_depth,
    log10_rate = log10(kramers_rate(max_barrier))
  )

cat("\n--- Barrier Summary ---\n")
print(barrier_summary)

write.csv(barrier_summary, file.path(data_dir, "barrier_stats.csv"), row.names = FALSE)
cat("\nSaved: barrier_stats.csv\n")

# ---- PMF overlay with barriers marked (S1) ---------------------------------
# S1 overlay plot — mark the global minimum and highest barrier
s1_colors <- c(
  "APO" = "#666666", "AZD5305" = "#E41A1C", "olaparib" = "#4DAF4A",
  "veliparib" = "#FF7F00", "niraparib" = "#984EA3",
  "rucaparib" = "#A65628", "talazoparib" = "#377EB8"
)

s1_pmf <- bind_rows(pmf_list[names(pmf_list) %in% s1_systems])

p_s1_pmf <- ggplot(s1_pmf, aes(x = RC, y = PMF, color = system)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = s1_colors) +
  labs(x = "HD-ART Distance (Å)", y = "PMF (kcal/mol)",
       title = "S1 PMF Landscapes with Conformational Barriers",
       color = NULL) +
  theme_7pt + theme(legend.position = "bottom", legend.key.size = unit(0.2, "cm"))

# ---- Barrier height comparison bar plot ----
p_barrier <- ggplot(barrier_summary, aes(x = reorder(system, max_barrier), 
                                          y = max_barrier, fill = architecture)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", max_barrier)), 
            hjust = -0.1, size = 2) +
  scale_fill_manual(values = c("S1" = "#2166AC", "S2" = "#B2182B")) +
  coord_flip() +
  labs(x = NULL, y = "Maximum Barrier Height (kcal/mol)",
       title = "Highest Conformational Barrier by System",
       fill = NULL) +
  theme_7pt + theme(legend.position = "bottom")

# ---- Barrier-to-depth ratio ----
p_ratio <- ggplot(barrier_summary, aes(x = reorder(system, max_barrier_frac),
                                         y = max_barrier_frac, fill = architecture)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", max_barrier_frac)),
            hjust = -0.1, size = 2) +
  scale_fill_manual(values = c("S1" = "#2166AC", "S2" = "#B2182B")) +
  coord_flip(ylim = c(0, max(barrier_summary$max_barrier_frac, na.rm = TRUE) * 1.15)) +
  labs(x = NULL, y = "Max Barrier / Well Depth",
       title = "Barrier-to-Depth Ratio",
       fill = NULL) +
  theme_7pt + theme(legend.position = "none")

# ---- Kramers rate ----
p_rate <- ggplot(barrier_summary, aes(x = reorder(system, -log10_rate),
                                        y = log10_rate, fill = architecture)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("S1" = "#2166AC", "S2" = "#B2182B")) +
  coord_flip() +
  labs(x = NULL, y = expression(log[10]~"k (s"^"-1"*")"),
       title = "Approximate Transition Rate (Kramers)",
       fill = NULL) +
  theme_7pt + theme(legend.position = "none")

# ---- Assemble figure ----
fig_barriers <- (p_s1_pmf / (p_barrier | p_ratio | p_rate)) +
  plot_annotation(
    title = "Conformational Transition Barrier Analysis",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  ) +
  plot_layout(heights = c(1.5, 1))

cairo_pdf(file.path(out_dir, "Fig_Transition_Barriers.pdf"),
          width = 190/25.4, height = 200/25.4, pointsize = 7)
print(fig_barriers)
dev.off()
cat("Saved: Fig_Transition_Barriers.pdf\n")

# ---- Key findings ----
cat("\n========== KEY BARRIER FINDINGS ==========\n")

# Compare Type II vs Type III in S1
type_assignment <- c(
  "APO" = "APO", "AZD5305" = "Type_II", "olaparib" = "Type_II",
  "veliparib" = "Type_III", "niraparib" = "Type_II",
  "rucaparib" = "Type_II", "talazoparib" = "Type_II"
)

barrier_summary$type <- type_assignment[barrier_summary$system]

s1_barriers <- barrier_summary[barrier_summary$architecture == "S1", ]
type2_s1 <- s1_barriers[s1_barriers$type == "Type_II", ]
type3_s1 <- s1_barriers[s1_barriers$type == "Type_III", ]

cat(sprintf("S1 Type II mean max barrier: %.1f kcal/mol (n=%d)\n",
            mean(type2_s1$max_barrier, na.rm = TRUE), nrow(type2_s1)))
cat(sprintf("S1 Type III mean max barrier: %.1f kcal/mol (n=%d)\n",
            mean(type3_s1$max_barrier, na.rm = TRUE), nrow(type3_s1)))

# AZD5305 specific
azd_s1 <- s1_barriers[s1_barriers$system == "AZD5305", ]
if (nrow(azd_s1) > 0) {
  cat(sprintf("AZD5305 S1: max barrier=%.1f, well depth=%.1f, ratio=%.2f\n",
              azd_s1$max_barrier, azd_s1$well_depth, azd_s1$max_barrier_frac))
}

message("===== Transition barrier analysis complete =====")
