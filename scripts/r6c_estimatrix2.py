#!/usr/bin/env python3
"""r6c_estimatrix2.py — final estimator x metric matrix on the corrected pipeline.
mm row  = frozen table-convention valid spans (master values, hardcoded).
q90/t95/eff = computed from the sampled-support curves (hist >= 1, empty bins excluded);
              S1 seeds use r6s1c_* dumps (canonical [240:720] window; per-seed values averaged).
"""
import numpy as np, itertools, csv
from pathlib import Path

AUD = Path("/Users/taozhu/.hermes/cache/scratch/audit")
OUT = Path("/Users/taozhu/clacky_workspace/PARPi_design/results/analysis/r6/r6c_estimatrix.csv")

def load_support(f):
    d = np.load(AUD / f, allow_pickle=True)
    h = np.asarray(d["hist"]); e = np.asarray(d["binsX"], float)
    rc = (e[:-1] + e[1:]) / 2; pmf = np.asarray(d["pmf_c3"], float)
    m = (h >= 1) & np.isfinite(pmf)
    return rc[m], pmf[m]

def variants(rc, pmf):
    q90 = np.percentile(pmf, 95) - np.percentile(pmf, 5)
    t95 = np.percentile(pmf, 97.5) - np.percentile(pmf, 2.5)
    sel = pmf <= pmf.min() + 3.0
    eff = float(rc[sel].max() - rc[sel].min()) if sel.sum() > 1 else np.nan
    return {"q90": q90, "t95": t95, "eff": eff}

MM = {
    ("S2CV1", "talazoparib"): 22.9, ("S2CV1", "niraparib"): 26.9, ("S2CV1", "olaparib"): 28.8,
    ("S2CV1", "rucaparib"): 28.7, ("S2CV1", "veliparib"): 28.1,
    ("S2CV1", "niraparib_prot"): 34.1, ("S2CV1", "rucaparib_prot"): 48.7,
    ("S2CV1", "veliparib_prot"): 25.5,
    ("S2CV2", "talazoparib"): 36.3, ("S2CV2", "niraparib"): 21.9, ("S2CV2", "olaparib"): 16.6,
    ("S2CV2", "rucaparib"): 20.1, ("S2CV2", "veliparib"): 41.2,
    ("S2CV2", "niraparib_prot"): 19.2, ("S2CV2", "rucaparib_prot"): 45.1,
    ("S2CV2", "veliparib_prot"): 16.8,
    ("S1", "talazoparib"): 42.2, ("S1", "niraparib"): 50.8, ("S1", "olaparib"): 55.4,
    ("S1", "rucaparib"): 27.6, ("S1", "veliparib"): 40.3,
    ("S1", "niraparib_prot"): 59.2, ("S1", "rucaparib_prot"): 27.7,
    ("S1", "veliparib_prot"): 23.7,
}

VAR = {}
for k in ["talazoparib", "veliparib", "niraparib", "niraparib_prot", "olaparib", "rucaparib", "rucaparib_prot"]:
    VAR[("S2CV1", k)] = variants(*load_support(f"r6c2_{k}_d0.1_c2_e20.npz"))
    VAR[("S2CV2", k)] = variants(*load_support(f"r6d_cv2_{k}_d0.1_c2_e20.npz"))
for k in ["talazoparib", "veliparib"]:
    sv = [variants(*load_support(f"r6s1c_{k}_{s}.npz")) for s in ("s7", "s49", "s123")]
    VAR[("S1", k)] = {e: float(np.nanmean([x[e] for x in sv])) for e in ("q90", "t95", "eff")}
for k in ["niraparib", "niraparib_prot", "olaparib", "rucaparib", "rucaparib_prot"]:
    VAR[("S1", k)] = variants(*load_support(f"r6s1_{k}.npz"))
# R9: veliparib(+1) dumps (protonation-state re-simulation, 2026-10-11)
VAR[("S2CV1", "veliparib_prot")] = variants(*load_support("veli_prot_s2_prot_CV1.npz"))
VAR[("S2CV2", "veliparib_prot")] = variants(*load_support("veli_prot_s2_prot_CV2.npz"))
VAR[("S1", "veliparib_prot")] = variants(*load_support("veli_prot_prot_s1.npz"))

def midrank(v):
    v = np.asarray(v, float); order = np.argsort(v, kind="stable"); r = np.empty(len(v)); i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        r[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return r

def rho_p(xs, ys):
    if any(not np.isfinite(x) for x in xs):
        return None
    rx = midrank(xs); ry = midrank(ys)
    if rx.std() == 0 or ry.std() == 0:
        return None
    obs = float(np.corrcoef(rx, ry)[0, 1])
    cnt = 0
    for perm in itertools.permutations(list(xs)):
        rx2 = midrank(list(perm))
        if rx2.std() == 0:
            continue
        if abs(float(np.corrcoef(rx2, ry)[0, 1])) >= abs(obs) - 1e-12:
            cnt += 1
    return obs, cnt / 120.0

names = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
Y = [5, 4, 2.5, 2.5, 1]
subs = {"neutral": {}, "nira_prot": {"niraparib": "niraparib_prot"},
        "ruca_prot": {"rucaparib": "rucaparib_prot"},
        "combined": {"niraparib": "niraparib_prot", "rucaparib": "rucaparib_prot",
                     "veliparib": "veliparib_prot"}}

rows = []
for pname, sub in subs.items():
    for metric in ("S2CV1", "S2CV2", "S1", "AAI"):
        for e in ("mm", "q90", "t95", "eff"):
            xs = []
            for k in names:
                kk = sub.get(k, k)
                if metric == "AAI":
                    num = MM[("S1", kk)] if e == "mm" else VAR[("S1", kk)][e]
                    den = MM[("S2CV2", kk)] if e == "mm" else VAR[("S2CV2", kk)][e]
                    xs.append(num / den if np.isfinite(den) and den != 0 else np.nan)
                else:
                    xs.append(MM[(metric, kk)] if e == "mm" else VAR[(metric, kk)][e])
            res = rho_p(xs, Y)
            if res is None:
                rows.append([pname, metric, e, "n.d.", "n.d."])
            else:
                rows.append([pname, metric, e, round(res[0], 3), round(res[1], 3)])

def show(pname):
    print(f"\n--- {pname} ---")
    for metric in ("S2CV1", "S2CV2", "S1", "AAI"):
        line = f"{metric:7s}"
        for e in ("mm", "q90", "t95", "eff"):
            r = [x for x in rows if x[0] == pname and x[1] == metric and x[2] == e][0]
            if r[3] == "n.d.":
                line += f"  {e}: n.d."
            else:
                line += f"  {e}: {r[3]:+.3f}({r[4]:.3f})"
        print(line)

for p in subs:
    show(p)

chk = [x for x in rows if x[0] == "neutral" and x[1] == "S2CV1" and x[2] == "mm"][0]
assert chk[3] == -0.667, chk
chk2 = [x for x in rows if x[0] == "neutral" and x[1] == "S1" and x[2] == "mm"][0]
assert chk2[3] == 0.308, f"S1 mm neutral: {chk2}"
chk3 = [x for x in rows if x[0] == "neutral" and x[1] == "AAI" and x[2] == "mm"][0]
assert chk3[3] == 0.205, f"AAI mm neutral: {chk3}"
print("\nkey mm rows verified: S2CV1 -0.667 / S1 +0.308 / AAI +0.205 ✓")
chk4 = [x for x in rows if x[0] == "combined" and x[1] == "S2CV1" and x[2] == "mm"][0]
assert chk4[3] == -0.205, f"combined S2CV1 mm: {chk4}"
print("combined (R9 triple) S2CV1 mm =", chk4[3], chk4[4], "✓")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["panel", "metric", "estimator", "rho", "p"]); w.writerows(rows)
print("WROTE", OUT)
