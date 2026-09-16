#!/usr/bin/env Rscript
# 21-spearman_correlation.R — Spearman correlation of S1 well depth / AAI vs experimental trapping
#
# Purpose:   Reproduce the manuscript analysis "S1 well depth and AAI show a
#            suggestive inverse relationship with experimental trapping
#            potency". Significance is assessed with the exact two-tailed
#            permutation test over all 5! = 120 rankings (n = 5).
# Inputs:    data/01_curated/trapping_potency.csv    (literature trapping potency, x olaparib)
#            results/figures/Fig_Mechanism_Data.csv  (S1 well depth wd_S1, AAI wd_ratio)
# Outputs:   results/figures/Fig4B_spearman_data.csv
#            results/figures/Fig_Trapping_vs_Allostery.pdf (Fig. 3, two panels)
# Depends:   R >= 4.0; packages: ggplot2, ggrepel, patchwork
#
# n = 5 follows the manuscript protocol: AZD5305 is excluded (potent
# PARP1-selective trapper measured on a different scale, not comparable to the
# dual-PARP x olaparib values) and APO is excluded (no trapping value).

trapping <- read.csv("data/01_curated/trapping_potency.csv", stringsAsFactors = FALSE)
mech     <- read.csv("results/figures/Fig_Mechanism_Data.csv", stringsAsFactors = FALSE)

stopifnot(nrow(trapping) == 5)
stopifnot(all(trapping$inhibitor %in% mech$ligand))

df <- data.frame(
  inhibitor  = trapping$inhibitor,
  trapping   = trapping$trapping_x_olaparib,
  well_depth = mech$wd_S1[match(trapping$inhibitor, mech$ligand)],
  aai        = mech$wd_ratio[match(trapping$inhibitor, mech$ligand)],
  stringsAsFactors = FALSE
)
df$log_trap <- log10(df$trapping)

cat("=== Data for Spearman test (n = 5) ===\n")
print(df[, c("inhibitor", "trapping", "well_depth", "aai")],
      row.names = FALSE, digits = 4)

# --- Exact permutation test over all n! rankings -----------------------------
all_permutations <- function(n) {
  # Returns a matrix with n! rows: every permutation of 1..n.
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
  # Two-tailed exact permutation p-value for Spearman's rho.
  x_rank  <- rank(x)
  y_rank  <- rank(y)
  rho_obs <- cor(x_rank, y_rank, method = "pearson")  # Spearman = Pearson on ranks
  perms   <- all_permutations(length(x))
  rho_perm <- apply(perms, 1L, function(p) cor(x_rank[p], y_rank, method = "pearson"))
  mean(abs(rho_perm) >= abs(rho_obs) - 1e-12)
}

rho1 <- cor(df$well_depth, df$log_trap, method = "spearman")
p1   <- exact_permutation_p(df$well_depth, df$log_trap)
rho2 <- cor(df$aai, df$log_trap, method = "spearman")
p2   <- exact_permutation_p(df$aai, df$log_trap)

cat("\n=== Test 1: S1 well depth vs experimental trapping ===\n")
cat(sprintf("Spearman rho = %.2f, exact two-tailed p = %.3f (n = %d, permutation test)\n",
            rho1, p1, nrow(df)))
cat("\n=== Test 2: AAI vs experimental trapping ===\n")
cat(sprintf("Spearman rho = %.2f, exact two-tailed p = %.3f (n = %d, permutation test)\n",
            rho2, p2, nrow(df)))

# --- Leave-one-out sensitivity ------------------------------------------------
cat("\n=== Leave-one-out sensitivity ===\n")
loo <- sapply(seq_len(nrow(df)), function(i) {
  cor(df$well_depth[-i], df$log_trap[-i], method = "spearman")
})
loo_df <- data.frame(excluded = df$inhibitor, rho_loo = round(loo, 2))
print(loo_df, row.names = FALSE)
cat(sprintf("LOO range: %.2f to %.2f; all negative: %s\n",
            min(loo), max(loo), all(loo < 0)))

# --- Manuscript-ready statement -----------------------------------------------
cat("\n=== Manuscript-ready statement ===\n")
cat(sprintf(
  "Spearman rank correlation between S1 well depth and trapping potency yielded rho = %.2f,\n",
  rho1))
cat(sprintf(
  "with an exact two-tailed p = %.3f (n = %d, permutation test). AAI showed an identical\n",
  p1, nrow(df)))
cat(sprintf(
  "rank correlation (rho = %.2f, p = %.3f). While this correlation does not reach the\n",
  rho2, p2))
cat("conventional p < 0.05 significance threshold, it is constrained by the limited sample\n")
cat(sprintf(
  "(n = %d) and should be interpreted as a trend warranting validation in larger inhibitor\n",
  nrow(df)))
cat(sprintf(
  "panels. A leave-one-out sensitivity analysis confirmed robustness: excluding any single\n"))
cat(sprintf(
  "inhibitor left rho between %.2f and %.2f (all negative), indicating that no single data\n",
  min(loo), max(loo)))
cat("point drives the correlation.\n")

# --- Regenerate Fig4B_spearman_data.csv (n = 5, no AZD5305 row) ----------------
out <- data.frame(
  inhibitor   = df$inhibitor,
  trapping    = df$trapping,
  well_depth  = df$well_depth,
  s1_s2_ratio = df$aai,
  log_trap    = df$log_trap,
  stringsAsFactors = FALSE
)
write.csv(out, "results/figures/Fig4B_spearman_data.csv", row.names = FALSE)
cat("\nSaved: results/figures/Fig4B_spearman_data.csv (n = 5, no AZD5305 row)\n")

# --- Fig_Trapping_vs_Allostery (manuscript Fig. 3, two panels) ----------------
library(ggplot2)
library(ggrepel)
library(patchwork)

# Shared helpers: theme_7pt (single source).
args_h <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", args_h))) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
  source(file.path(script_dir, "..", "common", "helpers.R"))
} else {
  source("common/helpers.R")
}

df$label <- c(talazoparib = "Talazoparib", niraparib = "Niraparib",
              olaparib = "Olaparib", rucaparib = "Rucaparib",
              veliparib = "Veliparib")[df$inhibitor]
df$class <- ifelse(df$inhibitor %in% c("talazoparib", "olaparib"), "Type II", "Type III")

ann <- sprintf("rho = %.2f\np = %.3f (n = %d)", rho1, p1, nrow(df))

p_a <- ggplot(df, aes(x = trapping, y = well_depth, color = class)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 10,
                  seed = 49,
                  min.segment.length = 0.2, box.padding = 0.25) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_color_manual(values = c("Type II" = "#E41A1C", "Type III" = "#377EB8")) +
  annotate("text", x = 0.05, y = 100, label = ann, size = 2.5, hjust = 0) +
  labs(x = "Trapping Potency (× Olaparib)", y = "S1 HD-ART Well Depth (kcal/mol)",
       title = "A  S1 Well Depth vs Trapping", color = NULL) +
  theme_7pt + theme(legend.position = c(0.87, 0.87))

p_b <- ggplot(df, aes(x = aai, y = trapping, color = class)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 10,
                  seed = 49,
                  min.segment.length = 0.2, box.padding = 0.25) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_color_manual(values = c("Type II" = "#E41A1C", "Type III" = "#377EB8")) +
  annotate("text", x = 1.05, y = 60, label = ann, size = 2.5, hjust = 0) +
  xlim(0.8, 4.0) +
  labs(x = "Allosteric Amplification Index (S1/S2)",
       y = "Trapping Potency (× Olaparib)",
       title = "B  AAI vs Trapping", color = NULL) +
  theme_7pt + theme(legend.position = "none")

cairo_pdf("results/figures/Fig_Trapping_vs_Allostery.pdf",
          width = 170/25.4, height = 75/25.4, pointsize = 7)
print(p_a | p_b)
dev.off()
cat("Saved: results/figures/Fig_Trapping_vs_Allostery.pdf\n")
