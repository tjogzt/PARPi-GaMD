#!/usr/bin/env python3
"""s54_state_combinations.py — complete eight-combination chemical-state lattice.

Purpose:  Extends the chemical-state panel table (Table S8) from the four reported
          panels (neutral, each single correction, all-three) to all 2^3 = 8
          combinations of the three audit-corrected ligands' neutral/protonated
          states (niraparib, rucaparib, veliparib). Descriptors use the frozen
          table convention (C3 cumulant, primary pipeline): S2 CV1 span,
          S2 CV2 span, and S1 span master values (identical to
          scripts/s38_estimator_matrix.py), with the span ratio computed as
          S1 / S2 CV2 per system. Statistics: tie-aware Spearman against the
          literature trapping ranks y = (5, 4, 2.5, 2.5, 1) (talazoparib >
          niraparib > olaparib = rucaparib > veliparib), exact two-sided
          permutation over all 120 data orderings (60 distinct tied vectors;
          minimum attainable p = 0.033). Asserts that the four panels published
          in Table S8 reproduce exactly.
Inputs:   none (frozen master values embedded; provenance documented)
Outputs:  data/analysis/r6/r9_state_lattice.csv
Depends:  numpy, itertools, csv
Run:      python3 scripts/s54_state_combinations.py   (from the repository root)
"""
import csv
import itertools

import numpy as np

# frozen master values (table convention; identical to scripts/s38_estimator_matrix.py)
SPAN = {
    ("S2CV1", "talazoparib"): 22.9, ("S2CV1", "niraparib"): 26.9,
    ("S2CV1", "niraparib_prot"): 34.1, ("S2CV1", "olaparib"): 28.8,
    ("S2CV1", "rucaparib"): 28.7, ("S2CV1", "rucaparib_prot"): 48.7,
    ("S2CV1", "veliparib"): 28.1, ("S2CV1", "veliparib_prot"): 25.5,
    ("S2CV2", "talazoparib"): 36.3, ("S2CV2", "niraparib"): 21.9,
    ("S2CV2", "niraparib_prot"): 19.2, ("S2CV2", "olaparib"): 16.6,
    ("S2CV2", "rucaparib"): 20.1, ("S2CV2", "rucaparib_prot"): 45.1,
    ("S2CV2", "veliparib"): 41.2, ("S2CV2", "veliparib_prot"): 16.8,
    ("S1", "talazoparib"): 42.2, ("S1", "niraparib"): 50.8,
    ("S1", "niraparib_prot"): 59.2, ("S1", "olaparib"): 55.4,
    ("S1", "rucaparib"): 27.6, ("S1", "rucaparib_prot"): 27.7,
    ("S1", "veliparib"): 40.3, ("S1", "veliparib_prot"): 23.7,
}

DRUGS = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
CORRECTABLE = ["niraparib", "rucaparib", "veliparib"]
Y = [5, 4, 2.5, 2.5, 1]  # trapping ranks (weakest-to-strongest coding)


def midrank(v):
    v = np.asarray(v, float)
    order = np.argsort(v, kind="stable")
    r = np.empty(len(v))
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        r[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return r


def rho_p(xs, ys=Y):
    rx = midrank(xs)
    ry = midrank(ys)
    obs = float(np.corrcoef(rx, ry)[0, 1])
    cnt = 0
    for perm in itertools.permutations(list(xs)):
        rx2 = midrank(list(perm))
        if rx2.std() == 0:
            continue
        if abs(float(np.corrcoef(rx2, ry)[0, 1])) >= abs(obs) - 1e-12:
            cnt += 1
    return obs, cnt / 120.0


def system_key(drug, prot_set):
    return drug + "_prot" if drug in prot_set else drug


PANELS = [
    ("neutral", frozenset()),
    ("niraparib+1", frozenset({"niraparib"})),
    ("rucaparib+1", frozenset({"rucaparib"})),
    ("veliparib+1", frozenset({"veliparib"})),
    ("niraparib+1,rucaparib+1", frozenset({"niraparib", "rucaparib"})),
    ("niraparib+1,veliparib+1", frozenset({"niraparib", "veliparib"})),
    ("rucaparib+1,veliparib+1", frozenset({"rucaparib", "veliparib"})),
    ("audit-corrected", frozenset(CORRECTABLE)),
]

# published Table S8 values (the four panels the manuscript reports); used as assertions
PUBLISHED = {
    ("neutral", "S1"): (0.308, 0.667), ("neutral", "S2CV1"): (-0.667, 0.267),
    ("neutral", "S2CV2"): (-0.051, 1.000), ("neutral", "ratio"): (0.205, 0.767),
    ("niraparib+1", "S1"): (0.462, 0.433), ("niraparib+1", "S2CV1"): (-0.051, 1.000),
    ("niraparib+1", "S2CV2"): (-0.205, 0.767), ("niraparib+1", "ratio"): (0.205, 0.767),
    ("rucaparib+1", "S1"): (0.308, 0.667), ("rucaparib+1", "S2CV1"): (-0.667, 0.267),
    ("rucaparib+1", "S2CV2"): (-0.308, 0.667), ("rucaparib+1", "ratio"): (0.308, 0.667),
    ("audit-corrected", "S1"): (0.616, 0.300), ("audit-corrected", "S2CV1"): (-0.205, 0.767),
    ("audit-corrected", "S2CV2"): (0.410, 0.500), ("audit-corrected", "ratio"): (-0.103, 0.900),
}

rows = []
for panel, protset in PANELS:
    vals = {}
    for metric in ("S2CV1", "S2CV2", "S1"):
        vals[metric] = [SPAN[(metric, system_key(d, protset))] for d in DRUGS]
    vals["ratio"] = [vals["S1"][i] / vals["S2CV2"][i] for i in range(5)]
    for metric in ("S1", "S2CV1", "S2CV2", "ratio"):
        rho, p = rho_p(vals[metric])
        rows.append(dict(panel=panel, metric=metric,
                         nira_state="protonated" if "niraparib" in protset else "neutral",
                         ruca_state="protonated" if "rucaparib" in protset else "neutral",
                         veli_state="protonated" if "veliparib" in protset else "neutral",
                         rho=round(rho, 3), p=round(p, 3)))

# consistency assertions against the published Table S8
bad = []
for r in rows:
    key = (r["panel"], r["metric"])
    if key in PUBLISHED:
        rho0, p0 = PUBLISHED[key]
        if abs(r["rho"] - rho0) > 1e-9 or abs(r["p"] - p0) > 1e-9:
            bad.append((key, (r["rho"], r["p"]), (rho0, p0)))
if bad:
    raise SystemExit(f"ASSERTION FAILED vs Table S8: {bad}")
print("4 published panels x 4 metrics reproduce exactly (16/16)")

print("\npanel                          | S1            | S2CV1        | S2CV2        | ratio")
for panel, _ in PANELS:
    cells = []
    for metric in ("S1", "S2CV1", "S2CV2", "ratio"):
        r = [x for x in rows if x["panel"] == panel and x["metric"] == metric][0]
        cells.append(f"{r['rho']:+.3f} ({r['p']:.3f})")
    print(f"{panel:30s} | " + " | ".join(cells))

s2 = [r for r in rows if r["metric"] == "S2CV1"]
print("\nS2 CV1 across the 8 panels: rho range "
      f"[{min(r['rho'] for r in s2):+.3f}, {max(r['rho'] for r in s2):+.3f}]; "
      f"max p = {max(r['p'] for r in s2):.3f}")
allp = [r["p"] for r in rows]
allrho = [r["rho"] for r in rows]
print(f"all 32 cells: p in [{min(allp):.3f}, {max(allp):.3f}], "
      f"rho in [{min(allrho):+.3f}, {max(allrho):+.3f}], "
      f"cells with p < 0.05: {sum(1 for p in allp if p < 0.05)}")

out = "data/analysis/r6/r9_state_lattice.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["panel", "metric", "nira_state", "ruca_state",
                                      "veli_state", "rho", "p"])
    w.writeheader()
    w.writerows(rows)
print("\nWROTE", out, len(rows), "cells")
