#!/bin/bash
# rebuild_all.sh — one-command rebuild + validation of the manuscript tables and figures.
#
# Chain (canonical order; dependencies flow data -> analysis tables -> figures):
#   [1] s02 seed dataset        [4] figures (13, 21, 12, 16)
#   [2] s03 trajectory metadata [5] s10 manifest checksum verification
#   [3] s08 consolidated S2 DBE + AAI table (golden check vs the frozen table)
#   [6] s57 manuscript-number audit (optional; needs the manuscript tex tree,
#       skipped automatically in checkouts that do not carry it)
#
# History: the pre-R6 steps (s04 dihedral energies, s05 textbook reweighting,
# s06 reconstruction, s07 xvg regeneration, s09 2D DBE grids) are superseded by
# the R6 corrected-span reweighting; their archived products are validated in
# place by s08 (golden check) and s10. The scripts remain in scripts/ for
# provenance — see docs/standardisation_report_v2.md.
#
# Requires: R (>=4.0; ggplot2/ggrepel/patchwork/data.table), python (numpy).
#           The trajectory data layer is read from $DATA_ROOT (archived volume).
# Usage:    DATA_ROOT=/path/to/PARPi_data bash scripts/rebuild_all.sh
#
# Purpose:  One-command rebuild + validation of the main tables and figures.
# Created:  2026-09-17 (header standardised 2026-10-05; chain-order rewrite 2026-10-06;
#           R6-supersession rewrite 2026-10-11)
# Run:      bash scripts/rebuild_all.sh   (from the repository root)
set -euo pipefail
cd "$(dirname "$0")/.."
: "${DATA_ROOT:?Set DATA_ROOT to the trajectory data directory}"
mkdir -p results results/analysis results/figures figures/pdf
export PYTHONHASHSEED=0
PY="${PYTHON:-}"
if [ -z "$PY" ]; then PY="$(command -v python3 || true)"; fi
if [ -z "$PY" ]; then echo "ERROR: python3 not found; set PYTHON=/abs/path/to/python3" >&2; exit 1; fi
export RETICULATE_PYTHON="$PY"

echo "== [1/6] trapping seed dataset =="
"$PY" scripts/s02_build_seed_dataset.py

echo "== [2/6] trajectory metadata (all archival systems) =="
"$PY" scripts/s03_build_trajectory_metadata.py

echo "== [3/6] consolidated S2 DBE + AAI table (golden check vs frozen) =="
"$PY" scripts/s08_consolidate_s2_dbe.py

echo "== [4/6] figures =="
Rscript code/13_mechanism_figure.R
Rscript code/21_spearman_correlation.R
Rscript code/12_descriptive_analysis.R
Rscript code/16_2d_landscape.R

echo "== [5/6] manifest verification =="
"$PY" scripts/s10_verify_manifest.py

echo "== [6/6] manuscript-number audit (optional) =="
TEX_DIR="${TEX_DIR:-../manuscripts/jcim}"
if [ -d "$TEX_DIR" ]; then
  "$PY" scripts/s57_manuscript_number_audit.py --tex-dir "$TEX_DIR"
else
  echo "s57 skipped: $TEX_DIR not present in this checkout;"
  echo "             number audit runs wherever the manuscript tex tree is available."
fi

echo "REBUILD COMPLETE. figures/pdf/*.pdf regenerated; manifest verified."
