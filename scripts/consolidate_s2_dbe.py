#!/usr/bin/env python3
"""Consolidate the final 10-system S2 DBE + AAI table (P0-2 reanalysis).

Writes results/analysis/s2_dbe_final.csv (per-system S2 CV1/CV2 C1-C3 well
depths + source: exact log DBE vs reconstructed+bias-corrected) and
results/analysis/aai_dbe_final.csv (AAI = S1-C3 / S2-CV2-C3).

Sources: 7 systems exact (scripts/reweight_s2_dbe.py, s2_dbe_well_depths.npy);
3 systems reconstructed (scripts/reconstruct_s2_dbe.py + C3 bias correction
+7.7 kcal, validated end-to-end on APO/talazoparib).
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir

A = analysis_dir()

EXACT = {
    # lig: {CV1: [C1,C2,C3], CV2: [C1,C2,C3]}  (from s2_dbe_well_depths.npy)
    "APO":        {"CV1": [34.0, 39.4, 47.6], "CV2": [34.0, 40.1, 44.0]},
    "talazoparib": {"CV1": [34.6, 42.5, 53.1], "CV2": [34.9, 43.1, 61.9]},
    "AZD5305":    {"CV1": [34.8, 45.1, 79.4], "CV2": [34.2, 40.9, 50.3]},
    "veliparib":  {"CV1": [34.6, 42.8, 58.9], "CV2": [34.4, 48.7, 72.9]},
    "fluzoparib": {"CV1": [33.2, 41.3, 56.9], "CV2": [33.5, 38.8, 46.9]},
    "pamiparib":  {"CV1": [34.1, 42.4, 58.1], "CV2": [34.0, 40.4, 51.8]},
    "senaparib":  {"CV1": [33.7, 44.0, 66.3], "CV2": [34.1, 42.5, 49.4]},
}
C3_BIAS = 7.7  # reconstructed-C3 bias correction (kcal/mol), see validation
RECON = {
    # raw reconstruction C1/C2/C3; C3 stored bias-corrected
    "niraparib": {"CV1": [31.2, 33.2, 34.4 + C3_BIAS], "CV2": [31.4, 33.5, 32.2 + C3_BIAS]},
    "olaparib":  {"CV1": [32.0, 35.5, 37.6 + C3_BIAS], "CV2": [33.1, 36.4, 39.8 + C3_BIAS]},
    "rucaparib": {"CV1": [32.2, 35.9, 39.3 + C3_BIAS], "CV2": [32.2, 36.5, 32.5 + C3_BIAS]},
}
S1_C3 = {"talazoparib": 30.4, "olaparib": 68.9, "niraparib": 51.6,
         "rucaparib": 53.4, "veliparib": 101.1, "AZD5305": 28.2,
         "APO": 42.5, "fluzoparib": 54.0, "pamiparib": 48.7, "senaparib": 52.2}

rows = []
for lig, vals in {**EXACT, **RECON}.items():
    src = "exact_log_dbe" if lig in EXACT else "reconstructed_bias_corrected"
    for cvn in ["CV1", "CV2"]:
        c1, c2, c3 = vals[cvn]
        rows.append([lig, cvn, c1, c2, c3, src])
    aai = S1_C3[lig] / vals["CV2"][2]
    rows.append([lig, "AAI", S1_C3[lig], vals["CV2"][2], aai, src])

with open(A / "s2_dbe_final.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["ligand", "metric", "C1", "C2", "C3", "source"])
    w.writerows(rows)
print(f"wrote {A/'s2_dbe_final.csv'} ({len(rows)} rows)")
