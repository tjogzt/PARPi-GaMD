#!/bin/bash
# download_replicates.sh — Pull replicate analysis data from the simulation host
# Usage: HOST=user@host PORT=22 REMOTE_ROOT=/path/to/runs bash download_replicates.sh
# No machine-specific credentials or paths are stored in this script; provide
# them via environment variables (required for reproducibility elsewhere).

#
# Purpose:  Pull the replicate analysis data from the simulation host via environment-variable parameters.
# Author:   Tao Zhu (tjogzt@gmail.com)
# Created:  2026-09-16 (header standardised 2026-10-05)
# Run:      bash code/48_download_replicates.sh   (from the repository root)
: "${HOST:?Set HOST (e.g. user@host) to the simulation host}"
: "${PORT:=22}"
: "${REMOTE_ROOT:?Set REMOTE_ROOT to the remote runs/ parent directory}"

# Download target: the repository's runs/ directory (repo-root-relative; no
# machine-specific absolute paths, so this script is portable).
LOCAL="$(cd "$(dirname "$0")/.." && pwd)/runs"

mkdir -p "$LOCAL"

LIGANDS=("rep_tala" "rep_veli" "rep_azd" "rep_eb47")
REPS=(1 2)

for lig in "${LIGANDS[@]}"; do
    for rep in "${REPS[@]}"; do
        TAG="${lig}_${rep}"
        REMOTE="${REMOTE_ROOT}/${TAG}/analysis/"
        LOCAL_DIR="${LOCAL}/${TAG}/"

        echo "=== Downloading ${TAG} ==="
        mkdir -p "$LOCAL_DIR"

        # PMF xvg files
        scp -P "$PORT" "${HOST}:${REMOTE}*pmf*.xvg" "$LOCAL_DIR/" 2>/dev/null && echo "  PMF: OK" || echo "  PMF: not found"

        # CV data npy
        scp -P "$PORT" "${HOST}:${REMOTE}*cv*.npy" "$LOCAL_DIR/" 2>/dev/null && echo "  CV: OK"

        # Weights
        scp -P "$PORT" "${HOST}:${REMOTE}*weights*" "$LOCAL_DIR/" 2>/dev/null

        # Done tag
        scp -P "$PORT" "${HOST}:${REMOTE_ROOT}/${TAG}/DONE" "$LOCAL_DIR/" 2>/dev/null && echo "  DONE ✓" || echo "  not finished"

        echo ""
    done
done

echo "All downloads complete. Run: Rscript code/18-replicate_analysis.R"
