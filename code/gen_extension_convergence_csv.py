#!/usr/bin/env python3
"""
gen_extension_convergence_csv.py — extension S1 convergence-curve data

Purpose:   Generate the cumulative well-depth convergence table (20 ns blocks)
           for the extension inhibitors as a CSV for the R figure script.
Inputs:    results/analysis/new_drugs/sys1_<drug>_{cv,weights}.dat
Outputs:   results/analysis/extension_s1_convergence.csv
Depends:   numpy; common (pmf, paths)
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir
from common.pmf import run_pyrew

DATA = analysis_dir() / "new_drugs"

rows = []
for d in ["fluzoparib", "pamiparib", "senaparib"]:
    cv = np.loadtxt(DATA / f"sys1_{d}_cv.dat")
    w = np.loadtxt(DATA / f"sys1_{d}_weights.dat")
    for t in range(20, 201, 20):
        n = t * 20  # 50 ps/frame
        r = run_pyrew(cv[:n], w[:n])
        rows.append({"ligand": d, "time_ns": t, "well_depth": r.c3})

out = analysis_dir() / "extension_s1_convergence.csv"
with open(out, "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=["ligand", "time_ns", "well_depth"])
    wr.writeheader()
    wr.writerows(rows)
print(f"wrote {out} ({len(rows)} rows)")
