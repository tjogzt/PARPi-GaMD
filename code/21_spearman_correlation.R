#!/usr/bin/env Rscript
# 21_spearman_correlation.R — two-state Spearman analysis, TABLE-CONVENTION values.
#
# Data sources (no in-code panel constants):
#   - data/01_curated/two_state_panel_values.csv — the Table S1/S5 C3 panel
#     values for all three states (neutral / niraparib-protonated /
#     rucaparib-protonated). AAI = S1 C3 / S2 CV2 C3 (Table 1 convention).
#   - data/analysis/s2_dbe_final.csv (from scripts/s08_consolidate_s2_dbe.py)
#     and data/analysis/s1_dbe_unified_wells.csv — the primary DBE pipeline
#     tables; the protonated panel columns are asserted against them on load.
#   - data/01_curated/trapping_potency.csv — trapping ranks + plot potencies.
# Neutral (original) and protonated niraparib panels. Tie-aware midranks,
# exact permutation test over all 5! = 120 rankings (n = 5).
# AZD5305 excluded (different scale), APO excluded (no trapping value).
#
# This script is ALSO the generator of association_input_table.csv (written at
# the end), the single machine-readable table for the association analysis.
#
# Purpose:  Two-state Spearman correlation analysis (n=5, tie-aware exact permutation), Table S1/S5 C3 convention.
# Created:  2026-09-15 (header standardised 2026-10-05)
# Inputs:   data/01_curated/two_state_panel_values.csv ; data/01_curated/trapping_potency.csv ; data/analysis/s2_dbe_final.csv ; data/analysis/s1_dbe_unified_wells.csv
# Outputs:  figures/pdf/Fig_Trapping_vs_Allostery.pdf ; data/analysis/association_input_table.csv
# Depends:  ggplot2, ggrepel, patchwork
# Run:      Rscript code/21_spearman_correlation.R   (from the repository root)
trapping <- read.csv("data/01_curated/trapping_potency.csv", stringsAsFactors = FALSE)
s2_dbe   <- read.csv("data/analysis/s2_dbe_final.csv", stringsAsFactors = FALSE)
ts <- read.csv("data/01_curated/two_state_panel_values.csv", stringsAsFactors = FALSE)
stopifnot(nrow(trapping) == 5, nrow(ts) == 5)
getc <- function(col) setNames(ts[[col]], ts$inhibitor)
s5 <- c("talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib")

# Panel vectors, all read from the curated table (Table S1/S5 C3 convention).
s2_cv2_c3   <- getc("s2cv2_neutral")      # S2 CV2 (HD–ART), neutral
s2_cv2_pro  <- getc("s2cv2_nira_prot")    # niraparib-protonated panel
s2_cv2_ruca <- getc("s2cv2_ruca_prot")    # rucaparib-protonated panel
s1_neu      <- getc("s1_neutral")         # S1 spans, neutral
s1_pro      <- getc("s1_nira_prot")       # niraparib-protonated panel
s1_ruca     <- getc("s1_ruca_prot")       # rucaparib-protonated panel
s2_cv1_neu  <- getc("s2cv1_neutral")      # S2 CV1 (protein–DNA), neutral
s2_cv1_c3   <- getc("s2cv1_nira_prot")    # niraparib-protonated panel
s2_cv1_ruca <- getc("s2cv1_ruca_prot")    # rucaparib-protonated panel
aai_ruca    <- (s1_ruca / s2_cv2_ruca)

# --- Assert the protonated panels against the primary DBE pipeline tables ----
chk1 <- setNames(round(s2_dbe$C3[s2_dbe$metric == "CV1"], 1), s2_dbe$ligand[s2_dbe$metric == "CV1"])
chk2 <- setNames(round(s2_dbe$C3[s2_dbe$metric == "CV2"], 1), s2_dbe$ligand[s2_dbe$metric == "CV2"])
s1u  <- read.csv("data/analysis/s1_dbe_unified_wells.csv", stringsAsFactors = FALSE)
s1u_c3 <- setNames(round(s1u$C3, 1), tolower(s1u$ligand))
stopifnot(all(s2_cv1_c3[s5] == chk1[s5]),
          all(s2_cv2_pro[s5] == chk2[s5]),
          all(s1_neu[s5]    == s1u_c3[s5]))

df <- data.frame(
  inhibitor = trapping$inhibitor,
  trap_rank = trapping$trapping_rank,
  s2_cv1    = s2_cv1_c3[trapping$inhibitor],      # protonated panel (Table S5)
  s2_neutral = s2_cv1_neu[trapping$inhibitor],    # original panel
  well_depth = s1_pro[trapping$inhibitor],        # protonated panel
  wd_neutral = s1_neu[trapping$inhibitor],        # original panel
  aai        = (s1_pro / s2_cv2_pro)[trapping$inhibitor],
  aai_neutral = (s1_neu / s2_cv2_c3)[trapping$inhibitor],
  stringsAsFactors = FALSE
)
df$trap_plot <- ts$trap_potency[match(df$inhibitor, ts$inhibitor)]
df$log_trap  <- log10(df$trap_plot)

cat("=== Table-convention values (n = 5) ===\n")
print(df[, c("inhibitor", "trap_rank", "s2_cv1", "s2_neutral",
             "well_depth", "wd_neutral", "aai", "aai_neutral")],
      row.names = FALSE, digits = 4)

# --- Exact permutation test ----------------------------------------------------
all_permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  prev <- all_permutations(n - 1L)
  out  <- matrix(0L, nrow = n * nrow(prev), ncol = n)
  idx  <- 0L
  for (i in seq_len(nrow(prev))) {
    for (pos in seq_len(n)) {
      idx <- idx + 1L
      out[idx, ] <- append(prev[i, ], n, after = pos - 1L)
    }
  }
  out
}

exact_permutation_p <- function(x, y) {
  x_rank  <- rank(x)
  y_rank  <- rank(y)
  rho_obs <- cor(x_rank, y_rank, method = "pearson")
  perms   <- all_permutations(length(x))
  rho_perm <- apply(perms, 1L, function(p) cor(x_rank[p], y_rank, method = "pearson"))
  mean(abs(rho_perm) >= abs(rho_obs) - 1e-12)
}

spearman_exact <- function(x, y) {
  list(rho = cor(x, y, method = "spearman"), p = exact_permutation_p(x, y))
}

t1 <- spearman_exact(df$s2_cv1, df$trap_rank)       # protonated panel
t2 <- spearman_exact(df$well_depth, df$trap_rank)
t3 <- spearman_exact(df$aai, df$trap_rank)
t1r <- spearman_exact(s2_cv1_ruca[df$inhibitor], df$trap_rank)  # rucaparib-prot panel
t2r <- spearman_exact(s1_ruca[df$inhibitor], df$trap_rank)
t3r <- spearman_exact(aai_ruca[df$inhibitor], df$trap_rank)
t1n <- spearman_exact(df$s2_neutral, df$trap_rank)  # neutral panel
t2n <- spearman_exact(df$wd_neutral, df$trap_rank)
t3n <- spearman_exact(df$aai_neutral, df$trap_rank)

cat(sprintf("\nProtonated: S2 rho=%.3f p=%.3f | S1 rho=%.3f p=%.3f | AAI rho=%+.3f p=%.3f\n",
            t1$rho, t1$p, t2$rho, t2$p, t3$rho, t3$p))
cat(sprintf("Neutral:    S2 rho=%.3f p=%.3f | S1 rho=%.3f p=%.3f | AAI rho=%+.3f p=%.3f\n",
            t1n$rho, t1n$p, t2n$rho, t2n$p, t3n$rho, t3n$p))
cat(sprintf("Rucap-prot: S2 rho=%.3f p=%.3f | S1 rho=%.3f p=%.3f | AAI rho=%+.3f p=%.3f\n",
            t1r$rho, t1r$p, t2r$rho, t2r$p, t3r$rho, t3r$p))

# Leave-one-out (protonated panel, primary metric)
loo <- sapply(seq_len(nrow(df)), function(i) {
  cor(df$s2_cv1[-i], df$trap_rank[-i], method = "spearman")
})
cat(sprintf("LOO (protonated): %+.3f to %+.3f; all negative: %s\n",
            min(loo), max(loo), all(loo < 0)))
print(data.frame(excluded = df$inhibitor, rho_loo = round(loo, 3)), row.names = FALSE)

# --- Figure --------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(patchwork)

args_h <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", args_h))) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
  source(file.path(script_dir, "..", "common", "helpers.R"))
} else {
  source("common/helpers.R")
}

df$label <- c(talazoparib = "Talazoparib", niraparib = "Niraparib\n(>1x, rank-only)",
              olaparib = "Olaparib", rucaparib = "Rucaparib",
              veliparib = "Veliparib")[df$inhibitor]
df$class <- ifelse(df$inhibitor %in% c("talazoparib", "olaparib"), "Type II", "Type III")

plot_df <- df
nira_neut <- df[df$inhibitor == "niraparib", ]; nira_neut$state <- "neutral"
nira_prot <- nira_neut; nira_prot$state <- "protonated"
nira_prot$s2_cv1 <- s2_cv1_c3["niraparib"]; nira_prot$well_depth <- s1_pro["niraparib"]
nira_prot$aai <- s1_pro["niraparib"] / s2_cv2_pro["niraparib"]
nira_neut$s2_cv1 <- s2_cv1_neu["niraparib"]; nira_neut$well_depth <- s1_neu["niraparib"]
nira_neut$aai <- s1_neu["niraparib"] / s2_cv2_c3["niraparib"]
others <- df[df$inhibitor != "niraparib", ]; others$state <- "neutral"
plot_df <- rbind(others, nira_neut, nira_prot)
# rucaparib dual state (protonated re-simulation)
ruca_neut <- df[df$inhibitor == "rucaparib", ]; ruca_neut$state <- "neutral"
ruca_prot <- ruca_neut; ruca_prot$state <- "protonated"
ruca_prot$s2_cv1 <- s2_cv1_ruca["rucaparib"]; ruca_prot$well_depth <- s1_ruca["rucaparib"]
ruca_prot$aai <- s1_ruca["rucaparib"] / s2_cv2_ruca["rucaparib"]
plot_df <- plot_df[plot_df$inhibitor != "rucaparib", ]
plot_df <- rbind(plot_df, ruca_neut, ruca_prot)
plot_df$is_nira <- plot_df$inhibitor == "niraparib"
plot_df$is_ruca <- plot_df$inhibitor == "rucaparib"
plot_df$state_f <- factor(plot_df$state, levels = c("neutral", "protonated"))

ann_a <- sprintf(
  "original:    rho = %.2f, p = %.3f\nnira +1:  rho = %.2f, p = %.3f\nruca +1:  rho = %.2f, p = %.3f\n(exact permutation, n = 5)",
  t1n$rho, t1n$p, t1$rho, t1$p, t1r$rho, t1r$p)
ann_b <- sprintf(
  "S1:   %.2f (p %.3f) / %.2f (p %.3f) / %.2f (p %.3f)\nAAI: %.2f / %.2f / %.2f (all n.s.)\n(neutral / nira +1 / ruca +1)",
  t2n$rho, t2n$p, t2$rho, t2$p, t2r$rho, t2r$p, t3n$rho, t3$rho, t3r$rho)

p_a <- ggplot(plot_df, aes(x = trap_plot, y = s2_cv1, color = class)) +
  geom_point(aes(shape = state_f), size = 3.0) +
  scale_shape_manual(values = c("neutral" = 16, "protonated" = 1), guide = "none") +
  geom_text_repel(aes(label = ifelse(is_nira & state == "protonated", "Niraparib-\nprotonated",
               ifelse(is_ruca & state == "protonated", "Rucaparib-\nprotonated", label))),
                  size = 2.9, max.overlaps = Inf, seed = 49,
                  min.segment.length = 0.3, box.padding = 0.35, force = 2, family = "Arial") +
  geom_segment(data = data.frame(x = 2, y1 = s2_cv1_neu["niraparib"], y2 = s2_cv1_c3["niraparib"]),
               aes(x = x, y = y1, xend = x, yend = y2),
               inherit.aes = FALSE, linetype = "dashed", linewidth = 0.3) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_color_manual(values = c("Type II" = "#E41A1C", "Type III" = "#377EB8")) +
  coord_cartesian(ylim = c(50, 66)) +
  labs(x = "Trapping Potency (x Olaparib)",
       y = "S2 CV1 Protein–DNA Span (kcal/mol)",
       title = "A  S2 Protein–DNA Span vs Trapping (two-state)",
       subtitle = "neutral: -0.82 (p 0.133)\nniraparib-protonated: -0.21 (p 0.767)\nrucaparib-protonated: -0.36 (p 0.633)",
       color = NULL) +
  theme_7pt + theme(legend.position = c(0.87, 0.87),
                    plot.subtitle = element_text(size = 8, hjust = 0)) +
  theme_open

p_b <- ggplot(plot_df, aes(x = aai, y = trap_plot, color = class)) +
  geom_point(aes(shape = state_f), size = 3.0) +
  scale_shape_manual(values = c("neutral" = 16, "protonated" = 1), guide = "none") +
  geom_text_repel(aes(label = ifelse(is_nira & state == "protonated", "Niraparib-\nprotonated",
               ifelse(is_ruca & state == "protonated", "Rucaparib-\nprotonated", label))),
                  size = 2.9, max.overlaps = Inf, seed = 49,
                  min.segment.length = 0.3, box.padding = 0.35, force = 2, family = "Arial") +
  geom_segment(data = data.frame(x1 = s1_neu["niraparib"] / s2_cv2_c3["niraparib"],
                                 x2 = s1_pro["niraparib"] / s2_cv2_pro["niraparib"], y = 2),
               aes(x = x1, y = y, xend = x2, yend = y),
               inherit.aes = FALSE, linetype = "dashed", linewidth = 0.3) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_color_manual(values = c("Type II" = "#E41A1C", "Type III" = "#377EB8")) +
  xlim(0.7, 1.3) +
  labs(x = "Allosteric Amplification Index (S1/S2)",
       y = "Trapping Potency (x Olaparib)",
       title = "B  AAI vs Trapping (sensitivity, two-state)",
       subtitle = "neutral: S1 -0.46 | AAI -0.36 (n.s.)\nniraparib-protonated: S1 -0.41 | AAI -0.05 (n.s.)\nrucaparib-protonated: S1 -0.82 | AAI +0.05 (n.s.)",
       color = NULL) +
  theme_7pt + theme(legend.position = "none",
                    plot.subtitle = element_text(size = 8, hjust = 0)) +
  theme_open

cairo_pdf("figures/pdf/Fig_Trapping_vs_Allostery.pdf",
          width = 157/25.4, height = 69.3/25.4, pointsize = 8)
print(p_a | p_b)
dev.off()
cat("Saved: figures/pdf/Fig_Trapping_vs_Allostery.pdf\n")

# Machine-readable input table (R1 remediation)
out <- data.frame(
  inhibitor   = df$inhibitor,
  trap_rank   = df$trap_rank,
  s2_cv1_prot = df$s2_cv1, s2_cv1_neutral = df$s2_neutral,
  s1_prot     = df$well_depth, s1_neutral = df$wd_neutral,
  aai_prot    = df$aai, aai_neutral = df$aai_neutral,
  stringsAsFactors = FALSE
)
write.csv(out, "data/analysis/association_input_table.csv", row.names = FALSE)
cat("Saved: data/analysis/association_input_table.csv\n")
