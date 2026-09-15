#!/usr/bin/env python3
"""
s1_convergence_blocks.py — S1 cumulative-block convergence

Purpose:   Cumulative well depth vs simulation length (20 ns blocks) for the
           extension inhibitors; verifies the manuscript's convergence claim.
Inputs:    results/analysis/new_drugs/sys1_<drug>_{cv,weights}.dat
Outputs:   console table
Depends:   numpy; common (pmf, paths)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir
from common.pmf import run_pyrew

DATA = analysis_dir() / "new_drugs"

drugs = ["fluzoparib", "pamiparib", "senaparib"]
print("=== Cumulative well depth (C3) ===")
print(f'{"drug":<12}' + "".join([f"{t:>9}" for t in
                                 ["20", "40", "60", "80", "100", "120", "140",
                                  "160", "180", "200ns"]]))
for d in drugs:
    cv = np.loadtxt(DATA / f"sys1_{d}_cv.dat")
    w = np.loadtxt(DATA / f"sys1_{d}_weights.dat")
    row = []
    for t in range(20, 201, 20):
        n = t * 20  # 50 ps/frame -> 20 ns = 400 frames
        r = run_pyrew(cv[:n], w[:n])
        row.append(f"{r.c3:>9.1f}")
    print(f"{d:<12}" + "".join(row))
