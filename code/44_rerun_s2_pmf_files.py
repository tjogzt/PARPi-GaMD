#!/usr/bin/env python3
"""
44_rerun_s2_pmf_files.py — regenerate AZD5305/veliparib S2 pmf-c3 files (DFW weights)

Purpose:   Rebuild the S2 pmf-c3 xvg files of the two rebuilt systems from the
           rerun data (DFW weights), replacing the 2.6 ns legacy values.
           Skip-guard: existing canonical pmf-c3 xvg files are never
           overwritten unless --force is passed.
Inputs:    data/new_drugs_s2/<drug>/{analysis_CV1,analysis_CV2}.dat,
           analysis_weights_dfw.dat
Outputs:   data/analysis/pmf-c3-sys2_<drug>_CV{1,2}_cv.dat.xvg
Depends:   numpy; common (pmf, paths)

Created:  2026-09-15 (header standardised 2026-10-05)
Run:      python3 code/44_rerun_s2_pmf_files.py   (from the repository root)
"""
import sys
import shutil
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir, REPO_ROOT
from common.pmf import run_pyrew

A = analysis_dir()
S2 = REPO_ROOT / "data" / "new_drugs_s2"

FORCE = "--force" in sys.argv[1:]

for d in ["AZD5305", "veliparib"]:
    w = np.loadtxt(S2 / d / "analysis_weights_dfw.dat")
    for cvn in ["CV1", "CV2"]:
        dst = A / f"pmf-c3-sys2_{d}_{cvn}_cv.dat.xvg"
        if dst.exists() and not FORCE:
            # canonical archived C3 curve present: never overwrite silently
            print(f"  [skip] {dst.name} exists (archived C3 curve kept; --force to overwrite)")
            continue
        cv = np.loadtxt(S2 / d / f"analysis_{cvn}.dat")
        r = run_pyrew(cv, w, wfmt="%.6f")
        shutil.copy(r.tmpdir / "pmf-c3-cv.dat.xvg", dst)
        print(f"{d} {cvn} -> {dst} (replaces the 2.6 ns legacy value)")
print("DONE")
