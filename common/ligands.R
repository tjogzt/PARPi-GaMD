# common/ligands.R — ligand panel metadata (single source: common/ligands.csv)
#
# Columns: ligand, label, class, color, lty, shape, trap_potency.
# trap_potency is the literature trapping potency relative to olaparib
# (empty = not determined; see data/01_curated/trapping_potency.csv).

load_ligands <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args)
  if (length(m)) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", args[m[1]])))
    p <- file.path(script_dir, "..", "common", "ligands.csv")
  } else {
    p <- "common/ligands.csv"
  }
  if (!file.exists(p)) {
    stop("cannot find common/ligands.csv (run scripts from the repository root)")
  }
  lig <- read.csv(p, stringsAsFactors = FALSE, na.strings = "")
  lig
}
