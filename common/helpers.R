# common/helpers.R — shared R helpers (single source of truth).
#
# Sourced by code/*.R scripts to keep per-script implementations from drifting
# (four-team code audit, 2026-09-16). Contents:
#   theme_7pt        — project-wide 7-pt figure theme (authoritative version,
#                      previously duplicated with variants in 12+ scripts)
#   extract_features — PMF descriptor extraction (canonical; was duplicated in
#                      12-descriptive_analysis.R and 28-pca_features.R)
#   read_pmf         — comment-aware xvg reader (header line counts vary across
#                      files; a fixed skip silently drops data rows)
#   source_common_helpers() — script-relative source helper for the above.

# Project-wide 7-pt figure theme (authoritative; formerly 13-mechanism_figure.R).
# Fully namespace-qualified so scripts without library(ggplot2) can still source us.
theme_7pt <- ggplot2::theme_bw(base_size = 7) +
  ggplot2::theme(
    panel.grid.minor  = ggplot2::element_blank(),
    panel.grid.major  = ggplot2::element_line(color = "grey92", linewidth = 0.2),
    plot.title        = ggplot2::element_text(size = 7, face = "bold"),
    axis.title        = ggplot2::element_text(size = 7),
    axis.text         = ggplot2::element_text(size = 6),
    legend.text       = ggplot2::element_text(size = 6),
    legend.title      = ggplot2::element_text(size = 7),
    legend.key.size   = grid::unit(0.3, "cm"),
    strip.text        = ggplot2::element_text(size = 7, face = "bold"),
    strip.background  = ggplot2::element_rect(fill = "grey95")
  )

# PMF descriptor extraction. Input: data.frame with columns RC and PMF_norm
# (min-normalized PMF). Returns a named numeric vector with 8 descriptors.
# n_states note: the "pmf <= pmf + 0.5" filter is tautological and never fires,
# so n_states counts ALL local minima (see data_manifest.md, n_states semantics).
extract_features <- function(pmf_df) {
  rc <- pmf_df$RC; pmf <- pmf_df$PMF_norm
  min_idx <- which.min(pmf)
  barrier_left  <- if (min_idx > 1) max(pmf[1:min_idx]) - pmf[min_idx] else NA
  barrier_right <- if (min_idx < length(rc)) max(pmf[min_idx:length(rc)]) - pmf[min_idx] else NA
  barrier <- max(barrier_left, barrier_right, na.rm = TRUE)
  hm <- max(pmf) / 2; idx <- which(pmf <= hm)
  fwhm <- if (length(idx) > 1) rc[tail(idx, 1)] - rc[idx[1]] else NA
  d1 <- diff(pmf); local_mins <- which(diff(sign(d1)) == 2) + 1
  n_states <- sum(pmf[local_mins] <= (pmf[local_mins] + 0.5))
  c(rc_min = rc[min_idx], well_depth = max(pmf) - min(pmf),
    rc_range = max(rc) - min(rc), barrier_left = barrier_left,
    barrier_right = barrier_right, barrier = barrier, fwhm = fwhm, n_states = n_states)
}

# Comment-aware xvg reader: drops lines starting with # or @ (header counts vary
# across files), then parses the two-column data. Returns data.frame(RC, PMF, PMF_norm).
read_pmf <- function(f) {
  lines <- readLines(f)
  lines <- lines[!grepl("^\\s*[#@]", lines)]
  d <- data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE)
  d <- as.data.frame(d); names(d) <- c("RC", "PMF")
  d$PMF_norm <- d$PMF - min(d$PMF)
  d
}

# Source this file from a script (handles both Rscript --file and interactive runs).
source_common_helpers <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args)
  if (length(m)) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", args[m[1]])))
    source(file.path(script_dir, "..", "common", "helpers.R"))
  } else {
    source("common/helpers.R")
  }
}
