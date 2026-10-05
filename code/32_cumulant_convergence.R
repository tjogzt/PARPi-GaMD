#!/usr/bin/env Rscript
# 32-cumulant_convergence.R — regenerates data/cumulant_convergence.csv
#
# The original producer of this CSV was lost; this script is the canonical
# reconstruction (2026-09-16 four-team audit, revised 2026-09-29 code audit).
# It is a PURE AGGREGATOR: all well-depth and C2-C3-correlation numbers are
# read from data/analysis/s1_table1_rows.csv (written by
# scripts/collect_s1_table1.py, the canonical S1 producer whose C1/C2/C3
# values match Table S1); only rc_min_C3 is taken from the mean-curve xvg
# files (data/analysis/sys1_*_pmf_c3.xvg), because s1_table1_rows.csv does
# not carry that column.
#
# Columns: wd_Ci = C1-C3 span, wd_mean/wd_sd/wd_range across orders,
# rc_min_C3, pmf_cor_C2C3 (= the collect_s1 "r" column, Table S1's r>0.80 claim).
# Note: the data/analysis sys1_*_pmf_c2/c3.xvg files were historically
# regenerated with a different reweighting protocol than Table S1 and must NOT
# be used to recompute wd or r (they produced r~0.2-0.8 for some systems and
# broke the Table S1 r>0.80 claim); the collect_s1 CSV is the single source.
#
# Purpose:  Regenerate data/cumulant_convergence.csv (per-window cumulant well depths used by the diagnostic figures).
# Author:   Tao Zhu (tjogzt@gmail.com)
# Created:  2026-09-16 (header standardised 2026-10-05)
# Inputs:   data/analysis
# Depends:  data.table
# Run:      Rscript code/32_cumulant_convergence.R   (from the repository root)
suppressMessages(library(data.table))

# Shared helpers: read_pmf (single source).
args_h <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", args_h))) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
  source(file.path(script_dir, "..", "common", "helpers.R"))
} else {
  source("common/helpers.R")
}

data_dir <- "data/analysis"
sys <- c("APO", "AZD5305", "niraparib", "olaparib", "rucaparib", "talazoparib", "veliparib")

# Canonical S1 producer output (Table S1 values; 3-seed means for tala/veli).
t1 <- fread(file.path(data_dir, "s1_table1_rows.csv"))
stopifnot(all(sys %in% t1$ligand))

out <- list()
for (lig in sys) {
  row <- t1[ligand == lig]
  p3 <- read_pmf(file.path(data_dir, sprintf("sys1_%s_pmf_c3.xvg", lig)))
  wd1 <- row$C1; wd2 <- row$C2; wd3 <- row$C3
  out[[lig]] <- data.frame(
    ligand = lig, wd_C1 = wd1, wd_C2 = wd2, wd_C3 = wd3,
    wd_mean = mean(c(wd1, wd2, wd3)), wd_sd = sd(c(wd1, wd2, wd3)),
    wd_range = max(c(wd1, wd2, wd3)) - min(c(wd1, wd2, wd3)),
    rc_min_C3 = p3$RC[which.min(p3$PMF)], pmf_cor_C2C3 = row$r)
}
res <- do.call(rbind, out)

# ---- Asserts against the legacy CSV (columns that must not drift) ----
legacy <- fread("data/cumulant_convergence.csv")
stopifnot(nrow(legacy) == 7)
for (lig in sys) {
  L <- legacy[ligand == lig]; N <- res[res$ligand == lig, ]
  for (col in c("wd_C1", "wd_C2", "wd_C3", "wd_mean", "wd_sd", "wd_range",
                "rc_min_C3", "pmf_cor_C2C3")) {
    # tolerance 0.1: legacy was written at 0.1 print precision; wd_sd/wd_range
    # propagate the ~0.05 rounding of wd_C1-C3, rc_min_C3 and r carry ~0.01
    # rounding differences from the current pipeline.
    stopifnot(abs(L[[col]] - N[[col]]) < 0.1)
  }
}

write.csv(res, "data/cumulant_convergence.csv", row.names = FALSE)
cat("Saved: data/cumulant_convergence.csv (all asserts passed)\n")
