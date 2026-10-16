#!/usr/bin/env bash
# Download the GaMD starting / reference PDB structures (mmCIF + PDB)
set -euo pipefail
DEST="data/00_raw/pdb"
mkdir -p "$DEST"
PDBS=(4DQY 6VKK 6VKO 6VKQ 6NTU 6BHV)
for p in "${PDBS[@]}"; do
  lc=$(echo "$p" | tr 'A-Z' 'a-z')
  for ext in cif pdb; do
    out="$DEST/${p}.${ext}"
    if [[ -s "$out" ]]; then echo "skip $out (exists)"; continue; fi
    url="https://files.rcsb.org/download/${p}.${ext}"
    if curl -fsSL --retry 3 --retry-delay 2 -o "$out" "$url"; then
      echo "OK   $out ($(wc -c < "$out") bytes)"
    else
      echo "FAIL $url"; rm -f "$out"
    fi
  done
done
echo "--- listing ---"
ls -lh "$DEST"
