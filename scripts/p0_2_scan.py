#!/usr/bin/env python3
"""p0_2_scan.py — truncation-length scan for the S2 CV1 association.

For each truncation length T in {8, 10, 14, 19.3, 23, 26} ns, truncate every
main-panel system to its first T ns of available data and recompute the n=5
Spearman association of the C3 span vs trapping rank (exact permutation p).

Frame intervals (verified): tala 1 ps; nira/ola/ruca 50 ps; veli 50 ps.

Purpose:  Truncation-length scan for the S2 CV1 association: recomputes the n=5 Spearman association at truncation lengths {8, 10, 14, 19.3, 23, 26} ns.
Created:  packaged 2026-10-06
Depends:  common.paths, common.pmf, numpy, MDAnalysis
Run:      python3 scripts/p0_2_scan.py   (from the repository root)
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

RANK = {"talazoparib": 5, "niraparib": 4, "olaparib": 2.5, "rucaparib": 2.5, "veliparib": 1}
ORDER = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
SCANS = [8, 10, 14, 19.3, 23, 26]


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
    if lig == "talazoparib":
        lg = np.loadtxt(D / "sys2_talazoparib" / "gamd.log", comments="#")
        cv = np.load(A / "sys2_talazoparib_CV1_cv.npy")
        return cv, lg[:, 7], 1.0
    if lig == "veliparib":
        lg = np.loadtxt(D / "s2_new_drugs" / "sys2_veliparib" / "gamd.log", comments="#")
        cv = np.loadtxt(A / "new_drugs_s2" / "veliparib" / "analysis_CV1.dat")
        off = len(lg) - len(cv)
        return cv, lg[off:, 7], 50.0
    lg = np.loadtxt(RERUN / "rerun_logs" / f"sys2_{lig}_gamd.log", comments="#")
    keep = lg[:, 1] >= PROD_STEP
    dv = lg[:, 7][keep]
    cv = s2_cv_mda(f"sys2_{lig}")
    n = min(len(cv), len(dv))
    return cv[-n:], dv[-n:], 50.0


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

pairs = {lig: build(lig) for lig in ORDER}

print(f"{'T_ns':>5}  " + "  ".join(f"{lig[:6]:>7}" for lig in ORDER) + f"  {'rho':>7}{'p':>8}")
import csv
rows = []
for T in SCANS:
    c3 = []
    for lig in ORDER:
        cv, dv, dt = pairs[lig]
        n_tr = min(int(round(T * 1000.0 / dt)), len(cv))
        if n_tr < 50:
            c3.append(np.nan)
            continue
        r = run_pyrew(cv[:n_tr], np.column_stack([BETA * dv[:n_tr], np.zeros(n_tr), dv[:n_tr]]),
                      wfmt="%.6f")
        c3.append(round(float(r.c3), 1))
    ok = [i for i, v in enumerate(c3) if not np.isnan(v)]
    if len(ok) == 5:
        rho, _ = spearmanr(np.array(c3), np.array([RANK[s] for s in ORDER]))
        p = exact_p(rho)
    else:
        rho, p = np.nan, np.nan
    print(f"{T:>5}  " + "  ".join(f"{v:>7.1f}" if not np.isnan(v) else f"{'--':>7}" for v in c3)
          + f"  {rho:>7.3f}{p:>8.4f}", flush=True)
    rows.append(dict(T_ns=T, **{lig: c3[i] for i, lig in enumerate(ORDER)},
                     rho=round(rho, 3) if not np.isnan(rho) else None,
                     p=round(p, 4) if not np.isnan(p) else None))

with open(A / "p0_2_trunc_scan.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    wr.writeheader()
    wr.writerows(rows)
print(f"Saved {A / 'p0_2_trunc_scan.csv'}")
