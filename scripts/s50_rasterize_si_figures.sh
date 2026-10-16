#!/bin/bash
# rasterize_si_figures.sh — regenerate the 300-dpi raster copies of the SI figures.
# Keeps the combined SI PDF under the submission size limit while preserving
# byte-identical appearance. Run after any change to the source figures.
# Usage: bash scripts/rasterize_si_figures.sh
set -euo pipefail
cd "$(dirname "$0")/.."

FIGURES=("Fig_DCCM_Allostery")

TMPD=$(mktemp -d /tmp/si_raster.XXXXXX)
mkdir -p "$TMPD"
for fig in "${FIGURES[@]}"; do
    src="results/figures/${fig}.pdf"
    dst="manuscripts/jcim/si_figures/${fig}.pdf"
    [ -f "$src" ] || { echo "missing $src"; exit 1; }
    mkdir -p "$(dirname "$dst")"
    rm -rf "$TMPD"/"$fig"
    pdftoppm -png -r 300 "$src" "$TMPD"/"$fig"
    pngs=("$TMPD"/"$fig"-*.png)
    if [ "${#pngs[@]}" -eq 1 ]; then
        sips -s format pdf "${pngs[0]}" --out "$dst" >/dev/null
        echo "rasterized $fig -> $dst"
    else
        echo "multi-page raster for $fig; only single-page supported" >&2
        exit 1
    fi
done
echo "SI raster copies regenerated."
