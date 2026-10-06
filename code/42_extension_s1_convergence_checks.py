#!/usr/bin/env python3
"""
42_extension_s1_convergence_checks.py — extension S1 convergence diagnostics

Purpose:   Two console diagnostics on the 200 ns extension-inhibitor S1
           trajectories (50 ps/frame), supporting the SI convergence and
           matched-length robustness claims:
             (1) cumulative C3 well depth in 20 ns blocks (convergence);
             (2) truncation of the 200 ns series to the legacy-panel production
                 lengths (19 / 22 / 31 ns) with recomputed well depths
                 (matched-length sensitivity).
           Merged from the former s1_convergence_blocks.py and
           s1_length_sensitivity.py (same input, single entry point).
Created:   2026-09-15 (original scripts); merged 2026-10-05 (repo standardisation)
Inputs:    data/new_drugs_s1/sys1_<lig>_{cv,weights}.dat
Outputs:  console tables (no files written)
Depends:   numpy; common.pmf (run_pyrew)
Run:       python3 code/42_extension_s1_convergence_checks.py   (repo root)
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from common.pmf import run_pyrew  # noqa: E402

DATA = ROOT / "data" / "new_drugs_s1"
(ROOT / "results").mkdir(exist_ok=True)
DRUGS = ["fluzoparib", "pamiparib", "senaparib"]

# ---- (1) cumulative well depth in 20 ns blocks ------------------------------
print("=== Cumulative well depth (C3) ===")
print(f'{"drug":<12}' + "".join(f"{t:>9}" for t in
      ["20", "40", "60", "80", "100", "120", "140", "160", "180", "200ns"]))
for d in DRUGS:
    cv = np.loadtxt(DATA / f"sys1_{d}_cv.dat")
    w = np.loadtxt(DATA / f"sys1_{d}_weights.dat")
    row = []
    for t in range(20, 201, 20):
        n = t * 20  # 50 ps/frame -> 20 ns = 400 frames
        r = run_pyrew(cv[:n], w[:n])
        row.append(f"{r.c3:>9.1f}")
    print(f"{d:<12}" + "".join(row))

# ---- (2) matched-length truncation sensitivity ------------------------------
# Truncation frames (50 ps/frame): 19 ns = 380, 22 ns = 440, 31 ns = 620
CUTOFFS = {"19ns": 380, "22ns": 440, "31ns": 620}

print("\n=== Extension S1 well depth: length-truncation sensitivity ===")
print(f'{"drug":<12}{"200ns(full)":>10}' + "".join(f"{k:>10}" for k in CUTOFFS))
for d in DRUGS:
    cv = np.loadtxt(DATA / f"sys1_{d}_cv.dat")
    w = np.loadtxt(DATA / f"sys1_{d}_weights.dat")
    full = run_pyrew(cv, w, tmp=str(ROOT / "results"))  # full 200 ns re-check
    out = [f"{full.c3:>10.1f}"]
    for _, n in CUTOFFS.items():
        r = run_pyrew(cv[:n], w[:n])
        out.append(f"{r.c3:>10.1f}")
    print(f"{d:<12}" + "".join(out))
print()
print("Reference (legacy C3): veliparib 101.1 / olaparib 68.9 / rucaparib 53.4 / "
      "niraparib 51.6 / APO 42.5 / AZD5305 28.2")
print("(legacy lengths 19-31 ns; truncated extension runs are protocol-comparable)")
