#!/bin/bash
# download_replicates.sh — Pull replicate analysis data from AutoDL
# Usage: bash download_replicates.sh

AUTODL="root@connect.westc.seetacloud.com"
PORT="22250"
LOCAL="$HOME/clacky_workspace/PARPi_design/runs"

mkdir -p "$LOCAL"

LIGANDS=("rep_tala" "rep_veli" "rep_azd" "rep_eb47")
REPS=(1 2)

for lig in "${LIGANDS[@]}"; do
    for rep in "${REPS[@]}"; do
        TAG="${lig}_${rep}"
        REMOTE="/root/autodl-tmp/PARPi_design/runs/${TAG}/analysis/"
        LOCAL_DIR="${LOCAL}/${TAG}/"
        
        echo "=== Downloading ${TAG} ==="
        mkdir -p "$LOCAL_DIR"
        
        # PMF xvg files
        scp -P "$PORT" "${AUTODL}:${REMOTE}*pmf*.xvg" "$LOCAL_DIR/" 2>/dev/null && echo "  PMF: OK" || echo "  PMF: not found"
        
        # CV data npy
        scp -P "$PORT" "${AUTODL}:${REMOTE}*cv*.npy" "$LOCAL_DIR/" 2>/dev/null && echo "  CV: OK"
        
        # Weights
        scp -P "$PORT" "${AUTODL}:${REMOTE}*weights*" "$LOCAL_DIR/" 2>/dev/null
        
        # Done tag
        scp -P "$PORT" "${AUTODL}:/root/autodl-tmp/PARPi_design/runs/${TAG}/DONE" "$LOCAL_DIR/" 2>/dev/null && echo "  DONE ✓" || echo "  not finished"
        
        echo ""
    done
done

echo "All downloads complete. Run: Rscript code/18-replicate_analysis.R"
