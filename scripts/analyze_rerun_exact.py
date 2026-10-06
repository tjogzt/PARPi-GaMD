#!/usr/bin/env python3
"""Exact reweighting of the six rerun trajectories (2026-09-26).

S2 (protein+DNA, 26 ns production, DBE):  CV1 = protein CA COM - DNA heavy
  atom COM; CV2 = HD(662-787) - ART(788-1011) CA COM.  Weights = canonical DBE
  e^{beta*dV_D} (log col 8; dV_T is identically 0 in these runs).
S1 (CAT-only, 24 ns production, dual boost):  CV = HD(1:119) - ART(120:345)
  CA COM (internal numbering, same as Table S1).  Weights = full lower-dual
  boost e^{beta*(dV_T + dV_D)} (log cols 7+8; the 2026-09-21 audit showed
  dV_T active in S1, so both columns enter the weight).

Production frames: total_nstep >= 6.0e6 (0.1M prep + 1M cMD + 0.1M eq-prep
+ 4.8M GaMD eq; production-start-step.txt = 6000000).  PyReweighting: amdweight_CE, T=300, disc 0.1, Emax 20,
cutoff 2 (identical to reweight_s2_dbe.py used for the manuscript values).

Output: results/analysis/rerun_cumulant_wells.csv

Purpose:  Exact DBE reweighting of the six re-run trajectories; produces rerun_cumulant_wells.csv (olaparib/rucaparib finals).
Created:  packaged 2026-10-06
Depends:  common.paths, common.pmf, numpy
Run:      DATA_ROOT=<DATA_ROOT> python3 scripts/analyze_rerun_exact.py   (from the repository root)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
RERUN = data_root() / "rerun_202609"
DATA = data_root()
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis"
PROD_STEP = 6.0e6
LIGS = ["niraparib", "olaparib", "rucaparib"]


def load_log(tag):
    log = np.loadtxt(RERUN / "rerun_logs" / f"{tag}_gamd.log", comments="#")
    return log, log[:, 1] >= PROD_STEP


def s1_cv_mda(tag):
    import MDAnalysis as mda
    u = mda.Universe(str(RERUN / "sys_files" / f"{tag}.prmtop"),
                     str(RERUN / "runs" / tag / "gamd_out" / "output.dcd"))
    hd = u.select_atoms("resid 1:119 and name CA")
    art = u.select_atoms("resid 120:345 and name CA")
    cv = []
    for ts in u.trajectory:
        cv.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    return np.array(cv)


def s2_cv_mda(tag):
    import MDAnalysis as mda
    lig = tag.split("_")[1]
    prmtop = DATA / f"sys2_{lig}" / f"sys2_{lig}.prmtop"
    u = mda.Universe(str(prmtop),
                     str(RERUN / "runs" / tag / "gamd_out" / "output.dcd"))
    hd = u.select_atoms("resid 662:787 and name CA")
    art = u.select_atoms("resid 788:1011 and name CA")
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    assert len(dna) > 0, "no DNA heavy atoms selected"
    cv1, cv2 = [], []
    for ts in u.trajectory:
        cv1.append(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
        cv2.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    return np.array(cv1), np.array(cv2)


def main():
    rows = []
    for lig in LIGS:
        # ---- S2 ----
        tag = f"sys2_{lig}"
        log, prod = load_log(tag)
        dv_d = log[:, 7][prod]
        cv1, cv2 = s2_cv_mda(tag)
        print(f"{tag}: dV_D mean {dv_d.mean():.2f} "
              f"CV1 {cv1[prod][:len(dv_d)].mean():.1f} "
              f"CV2 {cv2[prod][:len(dv_d)].mean():.1f} nprod={len(dv_d)}",
              flush=True)
        w = np.column_stack([BETA * dv_d, np.zeros(len(dv_d)), dv_d])
        for cvn, cv in [("CV1", cv1[prod]), ("CV2", cv2[prod])]:
            n = min(len(cv), len(w))
            r = run_pyrew(cv[:n], w[:n])
            rows.append([tag, cvn, round(r.c1, 1), round(r.c2, 1), round(r.c3, 1),
                         "rerun_exact_dbe"])
            print(f"    {cvn}: C1={r.c1:.1f} C2={r.c2:.1f} C3={r.c3:.1f}", flush=True)
        # ---- S1 ----
        tag = f"sys1_{lig}"
        log, prod = load_log(tag)
        dv = (log[:, 6] + log[:, 7])[prod]
        cv = s1_cv_mda(tag)
        print(f"{tag}: dV_T+dV_D mean {dv.mean():.2f} "
              f"CV {cv[prod][:len(dv)].mean():.1f} nprod={len(dv)}", flush=True)
        w = np.column_stack([BETA * dv, np.zeros(len(dv)), dv])
        n = min(len(cv[prod]), len(w))
        r = run_pyrew(cv[prod][:n], w[:n])
        rows.append([tag, "CV", round(r.c1, 1), round(r.c2, 1), round(r.c3, 1),
                     "rerun_exact_dualboost"])
        print(f"    CV: C1={r.c1:.1f} C2={r.c2:.1f} C3={r.c3:.1f}", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "rerun_cumulant_wells.csv", "w") as f:
        f.write("system,cv,C1,C2,C3,source\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    print("WROTE results/analysis/rerun_cumulant_wells.csv")


if __name__ == "__main__":
    main()
