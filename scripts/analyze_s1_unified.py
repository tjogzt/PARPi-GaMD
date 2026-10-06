#!/usr/bin/env python3
"""Unified S1 panel reweighting (24 ns production, dual-boost weights).

S1 CV = HD(1:119)-ART(120:345) CA COM (internal numbering).
Weights = e^{beta*(dV_T + dV_D)} (log cols 7+8; dV_T is active in S1, see
the 2026-09-21 audit).  All systems truncated/aligned to 24 ns production:
  - APO, AZD5305: original 24 ns trajectories (log @1 ps, 24000 frames)
  - talazoparib, veliparib: 2026-09-21 200 ns replicates (s7/s49/s123),
    first 480 production frames (= 24 ns @50 ps) per seed, reported as
    mean +/- SD over seeds
  - niraparib/olaparib/rucaparib: 2026-09-26 rerun (24 ns, exact dual boost)
PyReweighting: amdweight_CE, T=300, disc 0.1, Emax 20, cutoff 2.

Output: results/analysis/s1_unified_wells.csv

Purpose:  Unified S1 panel reweighting; produces the unified S1 well depths consumed by the AAI numerators.
Created:  packaged 2026-10-06
Depends:  common.paths, common.pmf, numpy
Run:      DATA_ROOT=<DATA_ROOT> python3 scripts/analyze_s1_unified.py   (from the repository root)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
DATA = data_root()
REPL = Path(__file__).resolve().parents[1] / "results" / "analysis" / "replicates"
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis"
RERUN = data_root() / "rerun_202609"


def pyrew(cv, w):
    r = run_pyrew(cv, w)
    return round(r.c1, 1), round(r.c2, 1), round(r.c3, 1)


def main():
    rows = []

    # ---- APO / AZD5305: original 24 ns, exact dual boost ----
    for tag, prm in [("APO", "sys1_APO.prmtop"),
                     ("AZD5305", "system.prmtop")]:
        import MDAnalysis as mda
        d = DATA / f"sys1_{tag}"
        u = mda.Universe(str(d / prm), str(d / "output.dcd"))
        hd = u.select_atoms("resid 1:119 and name CA")
        art = u.select_atoms("resid 120:345 and name CA")
        assert len(hd) and len(art), f"{tag}: HD/ART selection empty"
        cv = []
        for ts in u.trajectory:
            cv.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
        cv = np.array(cv)
        log = np.loadtxt(d / "gamd.log", comments="#")
        dv = log[:, 6] + log[:, 7]
        n = min(len(cv), len(dv))
        w = np.column_stack([BETA * dv[:n], np.zeros(n), dv[:n]])
        c1, c2, c3 = pyrew(cv[:n], w[:n])
        rows.append([tag, c1, c2, c3, "orig24_exact"])
        print(f"S1 {tag}: C1={c1} C2={c2} C3={c3} (orig 24ns)", flush=True)

    # ---- tala / veli: 3 seeds x first 480 frames ----
    for tag in ["talazoparib", "veliparib"]:
        c1s, c2s, c3s = [], [], []
        for s in ["s7", "s49", "s123"]:
            d = REPL / f"sys1_{tag}_{s}"
            cv = np.loadtxt(d / "cv.dat")[:480]
            w = np.loadtxt(d / "weights.dat")[:480]
            if w.ndim == 1:
                w = np.column_stack([w, np.zeros(len(w)), w])
            c1, c2, c3 = pyrew(cv, w)
            c1s.append(c1); c2s.append(c2); c3s.append(c3)
            print(f"S1 {tag}_{s} (24ns trunc): C1={c1} C2={c2} C3={c3}", flush=True)
        rows.append([tag,
                     f"{np.mean(c1s):.1f}±{np.std(c1s):.1f}",
                     f"{np.mean(c2s):.1f}±{np.std(c2s):.1f}",
                     f"{np.mean(c3s):.1f}±{np.std(c3s):.1f}",
                     "rep3x24_exact"])

    # ---- nira/ola/ruca: rerun 24 ns (already computed) ----
    rr = np.loadtxt(OUT / "rerun_cumulant_wells.csv", delimiter=",",
                    dtype=str, skiprows=1)
    for lig in ["niraparib", "olaparib", "rucaparib"]:
        row = rr[(rr[:, 0] == f"sys1_{lig}")][0]
        rows.append([lig.capitalize(), row[2], row[3], row[4], "rerun24_exact"])

    with open(OUT / "s1_unified_wells.csv", "w") as f:
        f.write("ligand,C1,C2,C3,source\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    print("WROTE results/analysis/s1_unified_wells.csv")


if __name__ == "__main__":
    main()
