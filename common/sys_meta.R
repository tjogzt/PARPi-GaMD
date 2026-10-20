# common/sys_meta.R — canonical system metadata for the PMF / panel figures.
#
# Single source of truth for per-system labels, classes, colours and linetypes
# (code/10_s2_pmf_figures.R and code/11_s1_s2_combined_pmf.R previously carried
# independent copies that could fork silently).
#
# sys_meta           : per-ligand attributes; row order = display order of the
#                      combined S1/S2 facet figure (code/11).
# sys_display_order  : panel order of the per-system S2 figure (code/10).
#
# Purpose:  Canonical label/class/colour/linetype metadata for the seven S2 systems.
# Created:  2026-10-06 (extracted from code/10 and code/11 in the P1 refactor)
.sm_csv <- if (file.exists("common/ligands.csv")) "common/ligands.csv" else {
  stop("cannot find common/ligands.csv (run scripts from the repository root)")
}
.li <- read.csv(.sm_csv, stringsAsFactors = FALSE, na.strings = "")
sys_meta <- .li[, c("ligand", "label", "class", "color", "lty")]
sys_meta$class <- gsub("_", " ", sys_meta$class)   # "Type II" / "Type III"
sys_display_order <- c("APO", "AZD5305", "olaparib", "talazoparib",
                       "niraparib", "rucaparib", "veliparib")
