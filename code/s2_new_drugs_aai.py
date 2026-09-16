#!/usr/bin/env python3
"""LEGACY (2026-09-17): extension AAI console values computed with the old DFW
S2 denominators. Superseded by the textbook-DBE reanalysis
(scripts/reweight_s2_dbe.py, scripts/reconstruct_s2_dbe.py,
scripts/consolidate_s2_dbe.py -> results/analysis/s2_dbe_final.csv), whose AAI
values are the ones cited in the revised manuscript. Kept for provenance of the
original extension-panel protocol only; do not cite its printed AAI values.
"""
"""
s2_new_drugs_aai.py — extension AAI = S1 well depth (C3, DBE) / S2 CV2 (C3, DFW)

Purpose:   Compute the Allosteric Amplification Index of the extension
           inhibitors from the reweighted C3 well depths.
Inputs:    results/analysis/new_drugs/sys1_<drug>_{cv,weights}.dat
           results/analysis/new_drugs_s2/<drug>/analysis_CV2.dat + analysis_weights_dfw.dat
Outputs:   console table
Depends:   numpy; common (pmf, paths)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir
from common.pmf import run_pyrew

S1DIR = analysis_dir() / "new_drugs"
S2DIR = analysis_dir() / "new_drugs_s2"


def s1_cumulants(drug):
    """S1 C1/C2/C3 (DBE weights, full 200 ns)."""
    cv = np.loadtxt(S1DIR / f"sys1_{drug}_cv.dat")
    w = np.loadtxt(S1DIR / f"sys1_{drug}_weights.dat")
    r = run_pyrew(cv, w)
    return [r.c1, r.c2, r.c3]


def s2_c3(drug):
    """S2 CV2 C3 (DFW weights)."""
    cv = np.loadtxt(S2DIR / drug / "analysis_CV2.dat")
    w = np.loadtxt(S2DIR / drug / "analysis_weights_dfw.dat")
    r = run_pyrew(cv, w, wfmt="%.6f")
    return r.c3


print(f'{"drug":<12}{"S1 C1":>8}{"S1 C2":>8}{"S1 C3":>8}{"S1 C3+-SD":>12}'
      f'{"S2 CV2":>8}{"AAI":>8}{"AAI+-SD":>10}')
for d in ["fluzoparib", "pamiparib", "senaparib"]:
    c1, c2, c3 = s1_cumulants(d)
    s2 = s2_c3(d)
    sd = np.std([c1, c2, c3])
    aai = c3 / s2              # manuscript convention: C3-based central value
    aai_sd = sd / s2           # uncertainty propagated from the C1-C3 spread
    print(f"{d:<12}{c1:>8.1f}{c2:>8.1f}{c3:>8.1f}{c3:>8.1f}+-{sd:.1f}"
          f"{s2:>8.1f}{aai:>8.2f}{aai_sd:>8.2f}")

print("\nManuscript AAI (original panel): tala 0.99 / AZD5305 0.99 / nira 1.68 / "
      "ruca 1.75 / ola 2.27 / veli 3.60 / APO 1.39")
