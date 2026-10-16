#!/usr/bin/env python3
#
# Purpose:  Rebuild the S1 niraparib PROTONATED C3 PMF xvg from the archived per-frame data (sys1_nira_prot_frames.csv.gz), same protocol as collect_s1_table1.py. Writes results/analysis/sys1_niraparib_pmf_c3_prot.xvg (new file, ...
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   <DATA_ROOT> archived per-frame data (sys1_nira_prot_frames.csv.gz)
# Outputs:  Writes results/analysis/sys1_niraparib_pmf_c3_prot.xvg (new file, so the
# Depends:  numpy, common
# Run:      python3 scripts/s23_rebuild_s1_nira_prot.py

"""Rebuild the S1 niraparib PROTONATED C3 PMF xvg from the archived per-frame
data (sys1_nira_prot_frames.csv.gz), same protocol as collect_s1_table1.py.
Writes results/analysis/sys1_niraparib_pmf_c3_prot.xvg (new file, so the
neutral file sys1_niraparib_pmf_c3.xvg is no longer overwritten by either
producer). Expected C3 span ~59.3 kcal/mol (Table S1, protonated).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
A = Path(__file__).resolve().parents[1] / "results" / "analysis"

frames = np.loadtxt(A.parent / "archive" / "per_frame" / "sys1_nira_prot_frames.csv.gz",
                    skiprows=1)  # header row is comma-separated, data rows space-separated
# columns: frame, step_ps, dV_D, dV_T, cv_hd_art_com
cv = frames[:, 4]
dv = frames[:, 2]
n = min(len(cv), len(dv))
w = np.column_stack([BETA * dv[:n], np.zeros(n), dv[:n]])
r = run_pyrew(cv[:n], w)
print(f"S1 nira-protonated from archive: C1={r.c1:.1f} C2={r.c2:.1f} C3={r.c3:.1f}")

# save the C3 profile (run_pyrew keeps tmpdir artifacts; rebuild from pmf text)
import shutil
src = r.tmpdir / "pmf-c3-cv.dat.xvg"
dst = A / "sys1_niraparib_pmf_c3_prot.xvg"
shutil.copyfile(src, dst)
print(f"WROTE {dst}")
