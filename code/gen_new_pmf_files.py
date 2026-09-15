#!/usr/bin/env python3
"""
gen_new_pmf_files.py — generate extension pmf-c3 files for figure scripts

Purpose:   Write the reweighted C3 PMF curves of the extension inhibitors to the
           locations expected by the figure scripts:
           - S1: results/analysis/sys1_<lig>_pmf_c3.xvg (DBE weights, 200 ns)
           - S2: results/analysis/pmf-c3-sys2_<lig>_CV1/CV2_cv.dat.xvg (DFW weights)
Inputs:    results/analysis/new_drugs/sys1_<lig>_{cv,weights}.dat
           results/analysis/new_drugs_s2/<lig>/analysis_{CV1,CV2}.dat + analysis_weights_dfw.dat
Outputs:   xvg files under results/analysis/
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
S1 = A / "new_drugs"
S2 = A / "new_drugs_s2"

for d in ["fluzoparib", "pamiparib", "senaparib"]:
    # S1 (DBE weights)
    cv = np.loadtxt(S1 / f"sys1_{d}_cv.dat")
    w = np.loadtxt(S1 / f"sys1_{d}_weights.dat")
    r = run_pyrew(cv, w, wfmt="%.6f")
    shutil.copy(r.tmpdir / "pmf-c3-cv.dat.xvg", A / f"sys1_{d}_pmf_c3.xvg")
    # S2 CV1 + CV2 (DFW weights)
    w2 = np.loadtxt(S2 / d / "analysis_weights_dfw.dat")
    for cvn in ["CV1", "CV2"]:
        cv2 = np.loadtxt(S2 / d / f"analysis_{cvn}.dat")
        r2 = run_pyrew(cv2, w2, wfmt="%.6f")
        shutil.copy(r2.tmpdir / "pmf-c3-cv.dat.xvg",
                    A / f"pmf-c3-sys2_{d}_{cvn}_cv.dat.xvg")
    print(f"{d}: S1 + S2 CV1/CV2 pmf-c3 generated")
print("DONE")
