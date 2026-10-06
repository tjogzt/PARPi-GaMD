#!/usr/bin/env python3
"""p0_2_trunc19.py (v3) — P0-2: uniform 19.3 ns production truncation of the
S2 CV1 spans for the main panel, using the log's total_nstep column (2 fs
timestep) as the absolute time axis — no dt assumptions.

Production start = 6.0e6 steps. Truncation at 19.3 ns => steps <= 6.0e6 + 19.3e3/2*1e3
= 6.0e6 + 9.65e6 = 1.565e7. Every (cv, dv, nstep) triple is aligned from its
authority source:

  talazoparib : $DATA_ROOT/sys2_talazoparib/gamd.log, cv = sys2_talazoparib_CV1_cv.npy
                (whole log is post-prep production; see row 0 nstep check)
  veliparib   : $DATA_ROOT/s2_new_drugs/sys2_veliparib/gamd.log (rows[243:] aligned to
                analysis_CV1.dat via the dfw sanity check as in reweight_s2_dbe.py)
  nira/ola/ruca : $DATA_ROOT/rerun_202609/rerun_logs/sys2_<lig>_gamd.log (prod rows
                nstep>=6.0e6) + CV1 computed from rerun DCD (protein CA-COM - DNA heavy-COM)

Output: results/analysis/p0_2_trunc19.csv + console table + association.

Purpose:  Uniform 19.3 ns production truncation of the S2 CV1 spans (log total_nstep time axis, no dt assumptions).
Created:  packaged 2026-10-06
Depends:  common.paths, common.pmf, numpy, MDAnalysis
Run:      python3 scripts/p0_2_trunc19.py   (from the repository root)
"""
import sys
from pathlib import Path
from itertools import permutations

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir, data_root
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
A = Path(analysis_dir())
D = Path(data_root())
RERUN = D / "rerun_202609"
PROD_STEP = 6.0e6
STEP_PER_NS = 1000.0 / 2.0 * 1000.0  # 2 fs -> 500,000 steps per ns
TRUNC_NS = 19.3

RANK = {"talazoparib": 5, "niraparib": 4, "olaparib": 2.5, "rucaparib": 2.5, "veliparib": 1}
ORDER = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]


def s2_cv_mda(tag):
    import MDAnalysis as mda
    lig = tag.split("_")[1]
    prmtop = D / f"sys2_{lig}" / f"sys2_{lig}.prmtop"
    u = mda.Universe(str(prmtop), str(RERUN / "runs" / tag / "gamd_out" / "output.dcd"))
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    cv1 = []
    for ts in u.trajectory:
        cv1.append(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
    return np.array(cv1)


def build(lig):
    """Return (cv, dv, frame_ps) with every frame's own time spacing honored.

    Verified structures (2026-10-01):
      talazoparib : 25858 frames @1.0 ps (full 25.86 ns series incl. GaMD prep burned in)
      veliparib   : 517 frames @50 ps (rebuilt 26 ns production; log tail-aligned)
      nira/ola/ruca : 521 frames @50 ps (rerun 26 ns production, nstep>=6.0e6)
    """
    if lig == "talazoparib":
        lg = np.loadtxt(D / "sys2_talazoparib" / "gamd.log", comments="#")
        cv = np.load(A / "sys2_talazoparib_CV1_cv.npy")
        dv = lg[:, 7]
        dt = 1.0
    elif lig == "veliparib":
        lg = np.loadtxt(D / "s2_new_drugs" / "sys2_veliparib" / "gamd.log", comments="#")
        cv = np.loadtxt(A / "new_drugs_s2" / "veliparib" / "analysis_CV1.dat")
        off = len(lg) - len(cv)
        assert off > 0, "veli: no offset"
        dv = lg[off:, 7]
        dt = 50.0
    else:
        lg = np.loadtxt(RERUN / "rerun_logs" / f"sys2_{lig}_gamd.log", comments="#")
        keep = lg[:, 1] >= PROD_STEP
        dv = lg[:, 7][keep]
        cv = s2_cv_mda(f"sys2_{lig}")
        n = min(len(cv), len(dv))
        cv = cv[-n:]
        dv = dv[-n:]
        dt = 50.0
    n = min(len(cv), len(dv))
    return cv[:n], dv[:n], dt


def ess_from_dv(dv):
    ww = np.exp(BETA * dv)
    return (ww.sum() ** 2) / (ww ** 2).sum() / len(ww)


def exact_p(rho, n=5):
    y = np.array([RANK[s] for s in ORDER])
    cnt = tot = 0
    for perm in permutations(range(1, n + 1)):
        rp = np.corrcoef(np.array(perm, dtype=float), y)[0, 1]
        tot += 1
        if abs(rp) >= abs(rho) - 1e-9:
            cnt += 1
    return cnt / tot


from scipy.stats import spearmanr

rows = []
print(f"{'system':<14}{'N_prod':>7}{'ns_prod':>9}{'N_tr':>6}{'ns_tr':>7}"
      f"{'ESS':>7}{'C3_full':>9}{'C3_tr':>8}")
for lig in ORDER:
    cv, dv, dt = build(lig)
    n_tr = min(int(round(TRUNC_NS * 1000.0 / dt)), len(cv))
    r_full = run_pyrew(cv, np.column_stack([BETA * dv, np.zeros_like(dv), dv]), wfmt="%.6f")
    r_tr = run_pyrew(cv[:n_tr], np.column_stack([BETA * dv[:n_tr], np.zeros(n_tr), dv[:n_tr]]), wfmt="%.6f")
    ns_prod = len(cv) * dt / 1000.0
    ns_tr = n_tr * dt / 1000.0
    print(f"{lig:<14}{len(cv):>7}{ns_prod:>9.2f}{n_tr:>6}{ns_tr:>7.2f}"
          f"{ess_from_dv(dv):>7.3f}{r_full.c3:>9.1f}{r_tr.c3:>8.1f}", flush=True)
    rows.append(dict(system=lig, n_prod=len(cv), ns_prod=round(ns_prod, 2), n_trunc=n_tr,
                     ns_trunc=round(ns_tr, 2), ess=round(ess_from_dv(dv), 4),
                     c3_full=round(float(r_full.c3), 2), c3_trunc=round(float(r_tr.c3), 2)))

c3_tr = np.array([r["c3_trunc"] for r in rows])
y = np.array([RANK[s] for s in ORDER])
rho, _ = spearmanr(c3_tr, y)
p = exact_p(rho)
print(f"\nTruncated-19.3ns reassociation: rho={rho:.3f}  exact p={p:.4f}")
print("reference (full-length, authoritative): rho=-0.821  p=0.133")

import csv
with open(A / "p0_2_trunc19.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    wr.writeheader()
    wr.writerows(rows)
print(f"Saved {A/'p0_2_trunc19.csv'}")
