#!/usr/bin/env python3
"""
s2_weight_protocol_test.py — decisive test of the two S2 weight schemes

Purpose:   Reproduce the manuscript C3 value (talazoparib S2 CV2 = 30.6) from the
           legacy S2 data under two reweighting schemes:
           (a) legacy:  weight = exp(beta * DFW)  [old weight files]
           (b) standard: weight = exp(beta * DBE)  [gamd.log column 8]
Inputs:    results/analysis/sys2_talazoparib_CV2_cv.dat
           $DATA_ROOT/sys2_talazoparib/gamd.log
Outputs:   console table
Depends:   numpy; common (pmf, paths); DATA_ROOT environment variable
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir, data_root
from common.pmf import run_pyrew

CV = analysis_dir() / "sys2_talazoparib_CV2_cv.dat"
LOG = data_root() / "sys2_talazoparib" / "gamd.log"
beta = 1.677571

cv = np.loadtxt(CV)
log = np.loadtxt(LOG)
n = min(len(cv), len(log))
cv, log = cv[:n], log[:n]
prod = log[:, 0] == 1
cv_p, dfw, dbe = cv[prod], log[prod, 5], log[prod, 7]
print(f"production frames: {prod.sum()}, "
      f"DFW {dfw.mean():.3f}+-{dfw.std():.3f}, DBE {dbe.mean():.2f}+-{dbe.std():.2f}")

w3a = np.column_stack([beta * dfw, np.zeros(len(cv_p)), dfw])   # legacy: col1 = beta*DFW
w3b = np.column_stack([beta * dbe, np.zeros(len(cv_p)), dbe])   # standard: col1 = beta*DBE
wa = run_pyrew(cv_p, w3a, wfmt="%.6f")
wb = run_pyrew(cv_p, w3b, wfmt="%.6f")
print(f"(a) DFW scheme: C1={wa.c1:.1f} C2={wa.c2:.1f} C3={wa.c3:.1f}")
print(f"(b) DBE scheme: C1={wb.c1:.1f} C2={wb.c2:.1f} C3={wb.c3:.1f}")
print("manuscript talazoparib S2 CV2 C3 = 30.6")
