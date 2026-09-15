#!/usr/bin/env python3
"""
s1_length_sensitivity.py — S1 sampling-length sensitivity analysis

Purpose:   Truncate the extension S1 trajectories (200 ns) to the legacy-panel
           production lengths and recompute well depths, checking the
           manuscript's matched-length robustness claim.
           Legacy production lengths (frames x 1 ps): APO 21.5 / AZD5305 25 /
           olaparib 19 / niraparib 31.2 / rucaparib 31.2 / veliparib 30.8 ns.
           Extension data: cv.dat + weights.dat (4001 frames x 50 ps = 200 ns).
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
# Truncation frames (50 ps/frame): 19 ns = 380, 22 ns = 440, 31 ns = 620
cutoffs = {"19ns": 380, "22ns": 440, "31ns": 620}


def c3_min_of(stdout):
    """Parse the PyReweighting pmf_min-c3 value from its stdout."""
    for line in stdout.splitlines():
        if "pmf_min-c3" in line:
            return float(line.split("=")[1])
    return None


print("=== Extension S1 well depth: length-truncation sensitivity ===")
print(f'{"drug":<12}{"200ns(full)":>10}' + "".join([f"{k:>10}" for k in cutoffs]))
for d in drugs:
    cv = np.loadtxt(DATA / f"sys1_{d}_cv.dat")
    w = np.loadtxt(DATA / f"sys1_{d}_weights.dat")
    full = run_pyrew(cv, w, tmp=str(DATA / d))  # full 200 ns (re-check of cached run)
    results = [f"{full.c3:>10.1f}"]
    for k, n in cutoffs.items():
        r = run_pyrew(cv[:n], w[:n])
        results.append(f"{r.c3:>10.1f}")
    print(f"{d:<12}" + "".join(results))
print()
print("Reference (legacy C3): veliparib 101.1 / olaparib 68.9 / rucaparib 53.4 / "
      "niraparib 51.6 / APO 42.5 / AZD5305 28.2")
print("(legacy lengths 19-31 ns; truncated extension runs are protocol-comparable)")
