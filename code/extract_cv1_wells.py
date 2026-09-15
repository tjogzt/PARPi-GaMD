#!/usr/bin/env python3
"""
extract_cv1_wells.py — S2 CV1/CV2 well depths from the C3 cumulant pipeline

Purpose:   Extract the S2 CV1 (protein-DNA) and CV2 (HD-ART) PMF well depths
           from the pmf-c3 reweighted xvg files, so the C3-pipeline CV1/CV2
           values used in the manuscript (Table S5 and the AAI denominators)
           have CSV provenance.
Inputs:    results/analysis/pmf-c3-sys2_<system>_CV{1,2}_cv.dat.xvg
Outputs:   results/analysis/cumulant_wells_C3.csv
Depends:   numpy. Run from the repository root.
"""
import csv
import glob
import os

import numpy as np

ANALYSIS_DIR = "results/analysis"
OUT_CSV = os.path.join(ANALYSIS_DIR, "cumulant_wells_C3.csv")

rows = []
for cv in ["CV1", "CV2"]:
    print(f"--- {cv} (C3 cumulant pipeline) ---")
    for f in sorted(glob.glob(os.path.join(ANALYSIS_DIR, f"pmf-c3-sys2_*_{cv}_cv.dat.xvg"))):
        base = os.path.basename(f)
        name = base.replace("pmf-c3-sys2_", "").replace(f"_{cv}_cv.dat.xvg", "")
        d = np.loadtxt(f, comments=["#", "@"])
        if d.ndim == 1:
            d = d.reshape(-1, 2)
        rc, pmf = d[:, 0], d[:, 1]
        imin = int(pmf.argmin())
        imax = int(pmf.argmax())
        well = float(pmf.max() - pmf.min())
        print(f"{name:<12} rc range {rc.min():.1f}-{rc.max():.1f}  "
              f"well depth {well:.2f}  (min@{rc[imin]:.1f}, max@{rc[imax]:.1f})")
        rows.append({
            "system": name,
            "cv": cv,
            "n_points": int(d.shape[0]),
            "rc_min": round(float(rc.min()), 4),
            "rc_max": round(float(rc.max()), 4),
            "well_depth": round(well, 4),
            "min_at_rc": round(float(rc[imin]), 4),
            "max_at_rc": round(float(rc[imax]), 4),
        })

with open(OUT_CSV, "w", newline="") as fh:
    writer = csv.DictWriter(
        fh, fieldnames=["system", "cv", "n_points", "rc_min", "rc_max",
                        "well_depth", "min_at_rc", "max_at_rc"])
    writer.writeheader()
    writer.writerows(rows)

print(f"\nSaved: {OUT_CSV} ({len(rows)} rows)")
