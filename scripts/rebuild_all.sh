#!/bin/bash
# rebuild_all.sh — one-command rebuild of main tables and figures (reviewer Section 8).
#
# Chain (canonical order, from data_manifest.md generating scripts):
#   data   -> analysis CSVs (python) -> figure CSVs -> figures (R)
# Usage:  DATA_ROOT=/path/to/PARPi_data bash scripts/rebuild_all.sh
# Requires: R (>=4.0, ggplot2/ggrepel/patchwork/data.table), python (MDAnalysis,
#           scipy), OpenMM for the dihedral-energy scripts (rebuild only).
set -euo pipefail
cd "$(dirname "$0")/.."
: "${DATA_ROOT:?Set DATA_ROOT to the trajectory data directory}"
PY="${PYTHON:-python3}"
export RETICULATE_PYTHON="$PY"

echo "== [1/8] trajectory metadata =="
"$PY" scripts/build_trajectory_metadata.py

echo "== [2/8] trapping seed dataset =="
"$PY" scripts/build_seed_dataset.py

echo "== [3/8] S2 DBE weights (exact logs + trajectory reconstruction) =="
"$PY" scripts/reweight_s2_dbe.py
# Reconstruction for niraparib/olaparib/rucaparib (no archived boost energies):
#   1. per-frame dihedral-group energy (OpenMM, stride 10) -> /tmp/dihed_<lig>.csv
#   2. regenerate_s2_pmf_xvgs.py recon_weights() (E1 = Vmax - 157.0 kcal/mol)
"$PY" scripts/dihedral_group_energy.py \
    "$DATA_ROOT/dbe_rebuild/niraparib_solute.prmtop" \
    "$DATA_ROOT/dbe_rebuild/niraparib_solute.dcd" /tmp/dihed_niraparib.csv 10
"$PY" scripts/dihedral_group_energy.py \
    "$DATA_ROOT/dbe_rebuild/olaparib_solute.prmtop" \
    "$DATA_ROOT/dbe_rebuild/olaparib_solute.dcd" /tmp/dihed_olaparib.csv 10
"$PY" scripts/dihedral_group_energy.py \
    "$DATA_ROOT/dbe_rebuild/rucaparib_solute.prmtop" \
    "$DATA_ROOT/dbe_rebuild/rucaparib_solute.dcd" /tmp/dihed_rucaparib.csv 10
"$PY" scripts/consolidate_s2_dbe.py

echo "== [4/8] PMF xvg regeneration (20 S2 sets) =="
"$PY" scripts/regenerate_s2_pmf_xvgs.py

echo "== [5/8] 2D landscapes (7 systems) =="
"$PY" scripts/compute_2d_dbe.py

echo "== [6/8] mechanism + trapping-correlation figures =="
Rscript code/13-mechanism_figure.R
Rscript code/21-spearman_correlation.R

echo "== [7/8] descriptive + 2D figure + landscape R figures =="
Rscript code/12-descriptive_analysis.R
Rscript code/16-2d_landscape.R

echo "== [8/8] verification =="
"$PY" scripts/verify_manifest.py

echo "REBUILD COMPLETE. Check: results/figures/*.pdf regenerated; manifest verified."
