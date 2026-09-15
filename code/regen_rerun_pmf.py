#!/usr/bin/env python3
"""
regen_rerun_pmf.py — regenerate AZD5305/veliparib S2 pmf-c3 files (DFW weights)

Purpose:   Rebuild the S2 pmf-c3 xvg files of the two rebuilt systems from the
           rerun data (DFW weights), replacing the 2.6 ns legacy values.
Inputs:    results/analysis/new_drugs_s2/<drug>/{analysis_CV1,analysis_CV2}.dat,
           analysis_weights_dfw.dat
Outputs:   results/analysis/pmf-c3-sys2_<drug>_CV{1,2}_cv.dat.xvg
Depends:   numpy; common (pmf, paths)
"""
import sys
import shutil
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir
from common.pmf import run_pyrew

A = analysis_dir()
S2 = A / "new_drugs_s2"

for d in ["AZD5305", "veliparib"]:
    w = np.loadtxt(S2 / d / "analysis_weights_dfw.dat")
    for cvn in ["CV1", "CV2"]:
        cv = np.loadtxt(S2 / d / f"analysis_{cvn}.dat")
        r = run_pyrew(cv, w, wfmt="%.6f")
        dst = A / f"pmf-c3-sys2_{d}_{cvn}_cv.dat.xvg"
        shutil.copy(r.tmpdir / "pmf-c3-cv.dat.xvg", dst)
        print(f"{d} {cvn} -> {dst} (replaces the 2.6 ns legacy value)")
print("DONE")
