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
sys_meta <- data.frame(
  ligand = c("APO", "AZD5305", "olaparib", "talazoparib", "veliparib", "niraparib", "rucaparib"),
  label  = c("APO", "AZD5305", "Olaparib", "Talazoparib", "Veliparib", "Niraparib", "Rucaparib"),
  class  = c("APO", "Unknown", "Type_II", "Type_II", "Type_III", "Type_III", "Type_III"),
  color  = c("grey40", "darkorange", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3"),
  lty    = c("dotted", "dashed", "solid", "solid", "solid", "solid", "solid"),
  stringsAsFactors = FALSE
)

sys_display_order <- c("APO", "AZD5305", "olaparib", "talazoparib",
                       "niraparib", "rucaparib", "veliparib")
