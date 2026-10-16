#!/usr/bin/env python3
#
# Purpose:  Niraparib protonated re-run: generate the PMF xvg files needed for the consistency matrix (S2 CV1/CV2 + S1).
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   (see docstring)
# Depends:  numpy, common, MDAnalysis
# Run:      python3 scripts/s27_nira_prot_xvgs.py

"""Niraparib protonated re-run: generate the PMF xvg files needed for the consistency matrix (S2 CV1/CV2 + S1)."""
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
BASE = data_root() / "rerun_202609/runs_nira_prot"
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis"
PROD_STEP = 6.0e6


def main():
    import MDAnalysis as mda
    # ---- S2 ----
    s2 = BASE / "sys2"
    log = np.loadtxt(s2 / "gamd.log", comments="#")
    prod = log[:, 1] >= PROD_STEP
    dv = log[:, 7][prod]
    u = mda.Universe(str(s2 / "niraparib_best_prot.prmtop"), str(s2 / "output.dcd"))
    hd = u.select_atoms("resid 662:787 and name CA")
    art = u.select_atoms("resid 788:1011 and name CA")
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    cv1, cv2 = [], []
    for ts in u.trajectory:
        cv1.append(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
        cv2.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    w = np.column_stack([BETA * dv, np.zeros(len(dv)), dv])
    for cvn, cv in [("CV1", np.array(cv1)[prod]), ("CV2", np.array(cv2)[prod])]:
        n = min(len(cv), len(w))
        r = run_pyrew(cv[:n], w[:n])
        dst = OUT / f"pmf-c3-sys2_niraparib_{cvn}_cv.dat.xvg"
        shutil.copyfile(r.tmpdir / "pmf-c3-cv.dat.xvg", dst)
        print(f"S2 {cvn}: C1={r.c1:.1f} C2={r.c2:.1f} C3={r.c3:.1f} -> {dst.name}")
    # ---- S1 ----
    s1 = BASE / "sys1"
    log = np.loadtxt(s1 / "gamd.log", comments="#")
    prod = log[:, 1] >= PROD_STEP
    dv = log[:, 7][prod]
    u1 = mda.Universe(str(s1 / "sys1_niraparib_prot.prmtop"), str(s1 / "output.dcd"))
    hd1 = u1.select_atoms("resid 1:119 and name CA")
    art1 = u1.select_atoms("resid 120:345 and name CA")
    cv = []
    for ts in u1.trajectory:
        cv.append(np.linalg.norm(hd1.center_of_mass() - art1.center_of_mass()))
    cv = np.array(cv)[prod]
    n = min(len(cv), len(dv))
    w1 = np.column_stack([BETA * dv[:n], np.zeros(n), dv[:n]])
    r = run_pyrew(cv[:n], w1)
    dst = OUT / "sys1_niraparib_pmf_c3.xvg"
    shutil.copyfile(r.tmpdir / "pmf-c3-cv.dat.xvg", dst)
    print(f"S1: C1={r.c1:.1f} C2={r.c2:.1f} C3={r.c3:.1f} -> {dst.name}")


if __name__ == "__main__":
    main()
