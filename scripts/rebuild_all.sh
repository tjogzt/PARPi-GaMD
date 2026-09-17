#!/bin/bash
# rebuild_all.sh — one-command rebuild of main tables and figures (reviewer Section 8).
#
# Chain (canonical order, from data_manifest.md generating scripts):
#   data   -> analysis CSVs (python) -> figure CSVs -> figures (R)
# Usage:  DATA_ROOT=/Volumes/tjogzt4T/PARPi_data bash scripts/rebuild_all.sh
# Requires: R (>=4.0, ggplot2/ggrepel/patchwork/data.table), python (MDAnalysis,
#           scipy), OpenMM for the dihedral-energy scripts (rebuild only).
set -euo pipefail
cd "$(dirname "$0")/.."
export DATA_ROOT="${DATA_ROOT:-/Volumes/tjogzt4T/PARPi_data}"
export RETICULATE_PYTHON="/opt/anaconda3/bin/python3"

echo "== [1/8] trajectory metadata =="
/opt/anaconda3/bin/python3 scripts/build_trajectory_metadata.py

echo "== [2/8] trapping seed dataset =="
/opt/anaconda3/bin/python3 scripts/build_seed_dataset.py

echo "== [3/8] S2 DBE weights (exact logs + trajectory reconstruction) =="
/opt/anaconda3/bin/python3 scripts/reweight_s2_dbe.py
/opt/anaconda3/bin/python3 scripts/reconstruct_s2_dbe.py
/opt/anaconda3/bin/python3 scripts/consolidate_s2_dbe.py

echo "== [4/8] PMF xvg regeneration (20 S2 sets) =="
/opt/anaconda3/bin/python3 scripts/regenerate_s2_pmf_xvgs.py

echo "== [5/8] 2D landscapes (7 systems) =="
/opt/anaconda3/bin/python3 scripts/compute_2d_dbe.py

echo "== [6/8] mechanism + trapping-correlation figures =="
Rscript code/13-mechanism_figure.R
Rscript code/21-spearman_correlation.R

echo "== [7/8] descriptive + 2D figure + landscape R figures =="
Rscript code/12-descriptive_analysis.R
Rscript code/16-2d_landscape.R

echo "== [8/8] verification =="
/opt/anaconda3/bin/python3 scripts/verify_manifest.py

echo "REBUILD COMPLETE. Check: results/figures/*.pdf regenerated; manifest verified."
