#!/usr/bin/env python3
#
# Purpose:  EB-47 S1 well depth under the primary (cumulant/DBE) protocol.
# Created:  2026-10 (header standardised 2026-10-11)
# Inputs:   (see docstring)
# Outputs:  Output: results/analysis/eb47_primary_wells.csv
# Depends:  csv, numpy, common
# Run:      python3 scripts/s20_eb47_primary_dbe.py

"""analyze_eb47_primary.py — EB-47 S1 well depth under the primary (cumulant/DBE) protocol.

CV: HD-ART CA COM from the archived per-frame distances (dist_raw.npy).
Weights: per-trajectory effective applied bias e^{beta*(dV_T + dV_D)}:
  - seed 7:  production dV_T == 0  -> reduces to e^{beta*dV_D}
  - seed 49: dV_T active           -> full e^{beta*(dV_T + dV_D)}
PyReweighting amdweight_CE, T=300, disc=0.1, Emax=20, cutoff=2
(identical to scripts/s19_s1_unified_dbe.py).

Output: results/analysis/eb47_primary_wells.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "analysis"

rows = []
for seed in ["s7", "s49"]:
    d = ROOT / "results" / "replicates" / f"rep_eb47_{seed}"
    cv = np.load(d / "dist_raw.npy").astype(float)
    log = np.loadtxt(d / "gamd.log", comments="#")
    dvT, dvD = log[:, 6], log[:, 7]
    n = min(len(cv), len(log))
    dv = dvT[:n] + dvD[:n]           # effective applied bias (total + dihedral)
    w = np.column_stack([BETA * dv, np.zeros(n), dv])
    res = run_pyrew(cv[:n], w)
    c1, c2, c3 = res.c1, res.c2, res.c3
    print(f"EB-47 {seed}: n={n}  <dV_T>={dvT[:n].mean():.2f}  <dV_D>={dvD[:n].mean():.2f}  "
          f"C1={c1:.1f}  C2={c2:.1f}  C3={c3:.1f}", flush=True)
    rows.append({"seed": seed, "n": n, "dV_T_mean": round(float(dvT[:n].mean()), 3),
                 "dV_D_mean": round(float(dvD[:n].mean()), 3),
                 "C1": c1, "C2": c2, "C3": c3})

c3s = [r["C3"] for r in rows]
mean_c3 = float(np.mean(c3s))
sd_c3 = float(np.std(c3s, ddof=1))
print(f"\nEB-47 mean C3 = {mean_c3:.1f} ± {sd_c3:.1f} kcal/mol (N=2)")

with open(OUT / "eb47_primary_wells.csv", "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["seed", "n", "dV_T_mean", "dV_D_mean", "C1", "C2", "C3"])
    writer.writeheader()
    writer.writerows(rows)
print(f"WROTE {OUT / 'eb47_primary_wells.csv'}")
