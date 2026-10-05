#!/usr/bin/env python3
"""
41_extension_panel_data_files.py — derived data files for the extension panel

Purpose:   Produce the derived data files for the three extension inhibitors
           (fluzoparib / pamiparib / senaparib) that the downstream figure
           scripts consume:
             (1) data/extension_s1_convergence.csv — cumulative S1 C3 well
                 depth in 20 ns blocks (200 ns trajectories, 50 ps/frame);
             (2) data/analysis/sys1_<lig>_pmf_c3.xvg — S1 HD-ART C3 PMF curves
                 (dihedral-boost-energy weights, 200 ns);
             (3) data/analysis/pmf-c3-sys2_<lig>_CV{1,2}_cv.dat.xvg — S2
                 CV1/CV2 C3 PMF curves (DFW weights, 22 ns).
           Merged from the former gen_extension_convergence_csv.py and
           gen_new_pmf_files.py (identical logic, single entry point).
Author:    Tao Zhu (tjogzt@gmail.com)
Created:   2026-09-16 (original scripts); merged 2026-10-05 (repo standardisation)
Inputs:    data/new_drugs_s1/sys1_<lig>_{cv,weights}.dat
           data/new_drugs_s2/<lig>/analysis_{CV1,CV2}.dat
           data/new_drugs_s2/<lig>/analysis_weights_dfw.dat
Outputs:   data/extension_s1_convergence.csv
           data/analysis/sys1_<lig>_pmf_c3.xvg
           data/analysis/pmf-c3-sys2_<lig>_CV{1,2}_cv.dat.xvg
Depends:   numpy; common.pmf (run_pyrew)
Run:       python3 code/41_extension_panel_data_files.py   (from the repo root)
"""
import csv
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from common.pmf import run_pyrew  # noqa: E402

S1_DIR = ROOT / "data" / "new_drugs_s1"
S2_DIR = ROOT / "data" / "new_drugs_s2"
OUT_DIR = ROOT / "data" / "analysis"
LIGANDS = ["fluzoparib", "pamiparib", "senaparib"]

# ---- (1) cumulative S1 convergence table (20 ns blocks) --------------------
rows = []
for d in LIGANDS:
    cv = np.loadtxt(S1_DIR / f"sys1_{d}_cv.dat")
    w = np.loadtxt(S1_DIR / f"sys1_{d}_weights.dat")
    for t in range(20, 201, 20):
        n = t * 20  # 50 ps per frame -> t ns = t*20 frames
        r = run_pyrew(cv[:n], w[:n])
        rows.append({"ligand": d, "time_ns": t, "well_depth": r.c3})

conv_csv = ROOT / "data" / "extension_s1_convergence.csv"
with open(conv_csv, "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=["ligand", "time_ns", "well_depth"])
    wr.writeheader()
    wr.writerows(rows)
print(f"wrote {conv_csv} ({len(rows)} rows)")

# ---- (2)(3) PMF curve files for the figure scripts --------------------------
for d in LIGANDS:
    # S1: dihedral-boost-energy weights, full 200 ns
    cv = np.loadtxt(S1_DIR / f"sys1_{d}_cv.dat")
    w = np.loadtxt(S1_DIR / f"sys1_{d}_weights.dat")
    r = run_pyrew(cv, w, wfmt="%.6f")
    shutil.copy(r.tmpdir / "pmf-c3-cv.dat.xvg", OUT_DIR / f"sys1_{d}_pmf_c3.xvg")
    # S2: DFW weights, CV1 + CV2
    w2 = np.loadtxt(S2_DIR / d / "analysis_weights_dfw.dat")
    for cvn in ["CV1", "CV2"]:
        cv2 = np.loadtxt(S2_DIR / d / f"analysis_{cvn}.dat")
        r2 = run_pyrew(cv2, w2, wfmt="%.6f")
        shutil.copy(r2.tmpdir / "pmf-c3-cv.dat.xvg",
                    OUT_DIR / f"pmf-c3-sys2_{d}_{cvn}_cv.dat.xvg")
    print(f"{d}: S1 + S2 CV1/CV2 pmf-c3 files generated")
print("DONE")
