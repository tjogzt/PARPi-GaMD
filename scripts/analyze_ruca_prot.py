#!/usr/bin/env python3
"""Rucaparib protonated re-run: exact DBE reweighting (S1 + S2).
Data: <DATA_ROOT>/rerun_202609/runs_ruca_prot/{sys1,sys2}/
Convention: identical to analyze_nira_prot.py -- dV_D = log col8 (index 7),
    production segment = production-start-step.txt (6.0e6), weight = exp(beta*dV_D).
Output: results/analysis/ruca_prot_dbe.csv

Purpose:  Protonated rucaparib re-run engine: exact DBE reweighting (S1 + S2); produces the rucaparib protonated-panel values.
Created:  packaged 2026-10-06
Depends:  common.paths, common.pmf, numpy
Run:      DATA_ROOT=<DATA_ROOT> python3 scripts/analyze_ruca_prot.py   (from the repository root)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
BASE = data_root() / "rerun_202609/runs_ruca_prot"
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis"

PROD_STEP = 6.0e6  # == production-start-step.txt


def load_dv(log_path):
    log = np.loadtxt(log_path, comments="#")
    prod = log[:, 1] >= PROD_STEP
    return log[:, 7][prod], prod


def s1_cv(prmtop, dcd):
    import MDAnalysis as mda
    u = mda.Universe(str(prmtop), str(dcd))
    hd = u.select_atoms("resid 1:119 and name CA")
    art = u.select_atoms("resid 120:345 and name CA")
    assert len(hd) > 0 and len(art) > 0, f"S1 selection failed: {len(hd)},{len(art)}"
    cv = []
    for ts in u.trajectory:
        cv.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    return np.array(cv)


def s2_cv(prmtop, dcd):
    import MDAnalysis as mda
    u = mda.Universe(str(prmtop), str(dcd))
    hd = u.select_atoms("resid 662:787 and name CA")
    art = u.select_atoms("resid 788:1011 and name CA")
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    print(f"  S2 selections: HD {len(hd)} ART {len(art)} prot {len(prot)} dna {len(dna)}",
          flush=True)
    assert len(hd) > 0 and len(art) > 0, "S2 HD/ART selection failed"
    assert len(dna) > 0, "no DNA heavy atoms selected"
    cv1, cv2 = [], []
    for ts in u.trajectory:
        cv1.append(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
        cv2.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    return np.array(cv1), np.array(cv2)


def main():
    rows = []
    # ---- S2 ----
    s2dir = BASE / "sys2"
    if (s2dir / "output.dcd").exists():
        dv, prod = load_dv(s2dir / "gamd.log")
        cv1, cv2 = s2_cv(s2dir / "rucaparib_best_prot.prmtop", s2dir / "output.dcd")
        cv1, cv2 = cv1[prod], cv2[prod]
        n = min(len(cv1), len(dv))
        print(f"S2: dV_D mean {dv[:n].mean():.2f} "
              f"CV1 {cv1[:n].mean():.1f} CV2 {cv2[:n].mean():.1f} nprod={n}", flush=True)
        w = np.column_stack([BETA * dv[:n], np.zeros(n), dv[:n]])
        for cvn, cv in [("CV1", cv1), ("CV2", cv2)]:
            r = run_pyrew(cv[:n], w)
            rows.append([f"sys2_rucaparib_prot", cvn, round(r.c1, 1), round(r.c2, 1),
                         round(r.c3, 1), "ruca_prot_dbe"])
            print(f"    S2 {cvn}: C1={r.c1:.1f} C2={r.c2:.1f} C3={r.c3:.1f}", flush=True)
    else:
        print("S2 data missing — skipped")
    # ---- S1 ----
    s1dir = BASE / "sys1"
    if (s1dir / "output.dcd").exists():
        dv, prod = load_dv(s1dir / "gamd.log")
        cv = s1_cv(s1dir / "sys1_rucaparib_prot.prmtop", s1dir / "output.dcd")
        cv = cv[prod]
        n = min(len(cv), len(dv))
        print(f"S1: dV_D mean {dv[:n].mean():.2f} CV mean {cv[:n].mean():.1f} nprod={n}",
              flush=True)
        w = np.column_stack([BETA * dv[:n], np.zeros(n), dv[:n]])
        r = run_pyrew(cv[:n], w)
        rows.append(["sys1_rucaparib_prot", "CV", round(r.c1, 1), round(r.c2, 1),
                     round(r.c3, 1), "ruca_prot_dbe"])
        print(f"    S1 CV: C1={r.c1:.1f} C2={r.c2:.1f} C3={r.c3:.1f}", flush=True)
    else:
        print("S1 data missing — skipped")
    out = OUT / "ruca_prot_dbe.csv"
    with open(out, "w") as f:
        f.write("system,cv,c1,c2,c3,method\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    print(f"WROTE {out}")


if __name__ == "__main__":
    main()
