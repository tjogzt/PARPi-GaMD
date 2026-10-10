#!/usr/bin/env python3
"""R6d: rebuilt ordinal sensitivity model (fully explicit; replaces MCMC + across-cumulant-SD inputs).

Model: latent span theta_i = x_i + eps_i, eps_i ~ N(0, s^2) i.i.d.
  x_i = audited C3 span/AAI (Table convention); s = ASSUMED uncertainty band
  (no replicate-based per-system error exists; s is an explicit assumption grid).
Statistic: rho = tie-aware Spearman(ranks(theta), trapping midranks), n = 5.
Posterior: fully explicit -> direct i.i.d. draws (1e6), MC standard errors reported.
Prior predictive (flat): exact enumeration over all 120 orderings.
Arcsine tilt: weight w ~ 1/sqrt(1-rho^2) on the draws (sensitivity only).
"""
import itertools
import numpy as np
from pathlib import Path

OUT = Path("/Users/taozhu/clacky_workspace/PARPi_design/results/analysis/r6")
OUT.mkdir(parents=True, exist_ok=True)

KS = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
Y = {"talazoparib": 5.0, "niraparib": 4.0, "olaparib": 2.5, "rucaparib": 2.5, "veliparib": 1.0}

BASE = {
    "S2CV1": {"talazoparib": 22.9, "niraparib": 26.9, "olaparib": 28.8, "rucaparib": 28.7, "veliparib": 28.1},
    "S1":    {"talazoparib": 42.2, "niraparib": 50.8, "olaparib": 55.4, "rucaparib": 27.6, "veliparib": 40.3},
    "AAI":   {"talazoparib": 1.16, "niraparib": 2.32, "olaparib": 3.34, "rucaparib": 1.37, "veliparib": 0.98},
}
PANELS = {
    "neutral":   {"S2CV1": {}, "S1": {}, "AAI": {}},
    "nira_prot": {"S2CV1": {"niraparib": 34.1}, "S1": {"niraparib": 59.2}, "AAI": {"niraparib": 3.08}},
    "ruca_prot": {"S2CV1": {"rucaparib": 48.7}, "S1": {"rucaparib": 27.7}, "AAI": {"rucaparib": 0.61}},
}
# uncertainty grids: spans absolute (kcal/mol); AAI relative fraction of value
S_GRID = {"S2CV1": [6.0, 12.0, 18.0], "S1": [6.0, 12.0, 18.0], "AAI": [0.15, 0.30, 0.45]}
MID = {"S2CV1": 12.0, "S1": 12.0, "AAI": 0.30}
N_DRAW = 1_000_000


def midrank(v):
    v = np.asarray(v, dtype=float)
    order = np.argsort(v, kind="mergesort")
    r = np.empty(len(v))
    r[order] = np.arange(1, len(v) + 1, dtype=float)
    # average tied ranks
    out = r.copy()
    sv = v[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        if j > i:
            out[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    return out


def spearman(a, b):
    ra, rb = midrank(a), midrank(b)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    return float((ra * rb).sum() / np.sqrt((ra * ra).sum() * (rb * rb).sum()))


def panel_vals(metric, panel):
    v = dict(BASE[metric])
    v.update(PANELS[panel][metric])
    return v


rows = []
print("── deterministic (s=0) checks vs r6c grid ──")
EXPECT = {("neutral", "S2CV1"): -0.667, ("nira_prot", "S2CV1"): -0.051, ("ruca_prot", "S2CV1"): -0.667,
          ("neutral", "S1"): +0.308, ("nira_prot", "S1"): +0.462, ("ruca_prot", "S1"): +0.308,
          ("neutral", "AAI"): +0.205, ("nira_prot", "AAI"): +0.205, ("ruca_prot", "AAI"): +0.308}
for panel in PANELS:
    for metric in ("S2CV1", "S1", "AAI"):
        v = panel_vals(metric, panel)
        rho0 = spearman([v[k] for k in KS], [Y[k] for k in KS])
        exp = EXPECT[(panel, metric)]
        ok = abs(rho0 - exp) < 0.002
        print(f"  [{panel:9s} {metric:6s}] rho = {rho0:+.4f}  (grid {exp:+.3f}) {'OK' if ok else '*** MISMATCH ***'}")
        assert ok, f"deterministic mismatch {panel}/{metric}"

print("\n── MC cells (direct draws, 1e6) ──")
for panel in PANELS:
    for metric in ("S2CV1", "S1", "AAI"):
        v = panel_vals(metric, panel)
        x = np.array([v[k] for k in KS])
        y = np.array([Y[k] for k in KS])
        for s in [0.0] + S_GRID[metric]:
            # per-system sd
            if metric == "AAI":
                sd = s * x
            else:
                sd = np.full(5, s)
            if s == 0.0:
                rho_mc = np.array([spearman(x, y)])
                p_neg, med, lo, hi = (rho_mc < 0).mean(), rho_mc[0], rho_mc[0], rho_mc[0]
                mcse = 0.0
                p_arc = float(rho_mc[0] < 0)
            else:
                rng = np.random.default_rng(49)
                draws = x[None, :] + rng.normal(0.0, sd[None, :], size=(N_DRAW, 5))
                # vectorized midrank-corr for the fixed y vector
                ry = midrank(y)
                ry = ry - ry.mean()
                rho_mc = np.empty(N_DRAW)
                yd = np.sqrt((ry * ry).sum())
                for i0 in range(0, N_DRAW, 200_000):
                    chunk = draws[i0:i0 + 200_000]
                    rx = np.argsort(np.argsort(chunk, axis=1, kind="mergesort"), axis=1).astype(float)
                    # ties inside chunk are measure-zero under continuous noise
                    rx -= rx.mean(axis=1, keepdims=True)
                    rho_mc[i0:i0 + 200_000] = (rx * ry).sum(axis=1) / (np.sqrt((rx * rx).sum(axis=1)) * yd)
                p_neg = float((rho_mc < 0).mean())
                med, lo, hi = [float(q) for q in np.percentile(rho_mc, [50, 2.5, 97.5])]
                mcse = float(np.sqrt(p_neg * (1 - p_neg) / N_DRAW))
                w = 1.0 / np.sqrt(np.maximum(1e-12, 1.0 - rho_mc ** 2))
                p_arc = float((w * (rho_mc < 0)).sum() / w.sum())
            rows.append((panel, metric, s, p_neg, med, lo, hi, mcse, p_arc))
            print(f"  [{panel:9s} {metric:6s} s={s:<5.2g}] P(rho<0)={p_neg:.3f} med={med:+.3f} CI[{lo:+.2f},{hi:+.2f}] MCSE={mcse:.4f} | arcsine {p_arc:.3f}")

print("\n── MC stability across 4 seeds (mid-grid cells) ──")
for panel in PANELS:
    for metric in ("S2CV1", "S1", "AAI"):
        v = panel_vals(metric, panel)
        x = np.array([v[k] for k in KS])
        s = MID[metric]
        sd = s * x if metric == "AAI" else np.full(5, s)
        ps = []
        for seed in (1, 2, 3, 4):
            rng = np.random.default_rng(seed)
            draws = x[None, :] + rng.normal(0.0, sd[None, :], size=(200_000, 5))
            ry = midrank(np.array([Y[k] for k in KS])); ry = ry - ry.mean()
            rx = np.argsort(np.argsort(draws, axis=1, kind="mergesort"), axis=1).astype(float)
            rx -= rx.mean(axis=1, keepdims=True)
            rho_mc = (rx * ry).sum(axis=1) / (np.sqrt((rx * rx).sum(axis=1)) * np.sqrt((ry * ry).sum()))
            ps.append((rho_mc < 0).mean())
        spread = max(ps) - min(ps)
        print(f"  [{panel:9s} {metric:6s}] seed P(rho<0): " + " ".join(f"{p:.3f}" for p in ps) + f"  |Δ|max={spread:.4f}")

print("\n── prior predictive (flat over latent orderings): exact enumeration ──")
y = np.array([Y[k] for k in KS])
vals = sorted({spearman(np.array(p, dtype=float), y) for p in itertools.permutations(range(1, 6))})
from collections import Counter
c = Counter(round(spearman(np.array(p, dtype=float), y), 6) for p in itertools.permutations(range(1, 6)))
tot = sum(c.values())
print(f"  distinct rho values: {len(c)}")
for rho_v, cnt in sorted(c.items()):
    print(f"    rho={rho_v:+.3f}  P={cnt}/{tot}={cnt/tot:.4f}")
p_neg_pp = sum(cnt for rho_v, cnt in c.items() if rho_v < 0) / tot
print(f"  flat-prior P(rho<0) = {p_neg_pp:.4f} (by symmetry); support includes +/-1 at {c.get(1.0,0)}/{tot} each")

# write CSV
with open(OUT / "r6d_ordinal_sensitivity.csv", "w") as f:
    f.write("panel,metric,s,P_rho_neg,median,ci_lo,ci_hi,mcse,P_rho_neg_arcsine_tilt\n")
    for r in rows:
        f.write(f"{r[0]},{r[1]},{r[2]},{r[3]:.4f},{r[4]:+.4f},{r[5]:+.4f},{r[6]:+.4f},{r[7]:.5f},{r[8]:.4f}\n")
print(f"\nWROTE {OUT/'r6d_ordinal_sensitivity.csv'}")
