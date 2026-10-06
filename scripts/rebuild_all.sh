#!/bin/bash
# rebuild_all.sh — one-command rebuild of the main tables and figures (reviewer Section 8).
#
# Chain (canonical order; dependencies flow data -> analysis CSVs -> figures):
#   [1] s02 seed dataset        [5] s08 consolidated S2 DBE + AAI table
#   [2] s03 trajectory metadata [6] s07 S2 PMF xvg regeneration (skip-guarded)
#   [3] s05 S2 DBE well depths  [7] s09 2D DBE landscapes
#   [4] s04 dihedral energies   [8] figures (13, 21, 12, 16) + s10 manifest check
#
# Requires: R (>=4.0; ggplot2/ggrepel/patchwork/data.table), python (MDAnalysis,
#           numpy, scipy), OpenMM for the dihedral-energy step. The trajectory
#           data layer is read from $DATA_ROOT (archived volume).
# Usage:    DATA_ROOT=/path/to/PARPi_data bash scripts/rebuild_all.sh
#
# Purpose:  One-command rebuild of the main tables and figures from the archived data layer.
# Created:  2026-09-17 (header standardised 2026-10-05; chain-order rewrite 2026-10-06)
# Run:      bash scripts/rebuild_all.sh   (from the repository root)
set -euo pipefail
cd "$(dirname "$0")/.."
: "${DATA_ROOT:?Set DATA_ROOT to the trajectory data directory}"
mkdir -p results results/analysis results/figures figures/pdf
# Resolve an absolute interpreter path: reticulate requires an absolute
# RETICULATE_PYTHON (a bare "python3" is rejected as non-existent).
PY="${PYTHON:-}"
if [ -z "$PY" ]; then
  PY="$(command -v python3 || true)"
fi
if [ -z "$PY" ]; then
  echo "ERROR: python3 not found on PATH; set PYTHON=/abs/path/to/python3" >&2
  exit 1
fi
export RETICULATE_PYTHON="$PY"
# Deterministic hashing (defensive; no helper in this chain depends on set order).
export PYTHONHASHSEED=0

echo "== [1/8] trapping seed dataset =="
"$PY" scripts/s02_build_seed_dataset.py

echo "== [2/8] trajectory metadata (all archival systems) =="
"$PY" scripts/s03_build_trajectory_metadata.py

echo "== [3/8] S2 DBE well depths (exact archived boost logs) =="
"$PY" scripts/s05_reweight_s2_dbe.py

echo "== [4/8] dihedral-group energies (reconstruction support; OpenMM) =="
# positional-arg CLI; invoked here with explicit arguments (mechanical-run exemption)
"$PY" scripts/s04_dihedral_group_energy.py \
    "$DATA_ROOT/dbe_rebuild/niraparib_solute.prmtop" \
    "$DATA_ROOT/dbe_rebuild/niraparib_solute.dcd" results/dihed_niraparib.csv 10
"$PY" scripts/s04_dihedral_group_energy.py \
    "$DATA_ROOT/dbe_rebuild/olaparib_solute.prmtop" \
    "$DATA_ROOT/dbe_rebuild/olaparib_solute.dcd" results/dihed_olaparib.csv 10
"$PY" scripts/s04_dihedral_group_energy.py \
    "$DATA_ROOT/dbe_rebuild/rucaparib_solute.prmtop" \
    "$DATA_ROOT/dbe_rebuild/rucaparib_solute.dcd" results/dihed_rucaparib.csv 10

echo "== [5/8] consolidated S2 DBE + AAI table =="
"$PY" scripts/s08_consolidate_s2_dbe.py

echo "== [6/8] S2 PMF xvg regeneration (skip-guarded; keeps archived C3 curves) =="
"$PY" scripts/s07_regenerate_s2_pmf_xvgs.py

echo "== [7/8] 2D DBE landscapes =="
"$PY" scripts/s09_compute_2d_dbe.py

echo "== [8/8] figures + manifest verification =="
Rscript code/13_mechanism_figure.R
Rscript code/21_spearman_correlation.R
Rscript code/12_descriptive_analysis.R
Rscript code/16_2d_landscape.R
"$PY" scripts/s10_verify_manifest.py

echo "REBUILD COMPLETE. Check: figures/pdf/*.pdf regenerated; manifest verified."
