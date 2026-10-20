#!/usr/bin/env python3
"""s60_rank_encoding_seed_checks.py — reviewer statistical batch, items 11/13/20.

A. Archive-consistency checks: production-only 50-ps spans recomputed for
   APO / rucaparib-neutral against their published Table S5 values, and 2-ns
   window SDs for the systems added to the nearest-neighbour comparison
   (veliparib; talazoparib comes from s59).
B. Phenotype-rank encoding audit: every chemically defensible trapping-rank
   encoding crossed with the eight chemical-state lattice cells; exact
   two-tailed permutation p for each (metric: S2 CV1 C3).
C. Seed-summary sensitivity for niraparib: replace the single production value
   with the three-seed mean or median in the neutral and protonated panels.

Writes data/analysis/r6/s60_rank_encodings.csv and s60_seed_sensitivity.csv.
"""
import contextlib
import gzip
import importlib.util
import io
import itertools
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "pyrew1d", ROOT / "tools/PyReweighting/PyReweighting-1D.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

BETA = 1.0 / (0.001987 * 300.0)
DISC, EMAX, CUT, T = 0.1, 20.0, 2, 300.0


def spans(cv, dv):
    data, dV = np.asarray(cv, float), np.asarray(dv, float)
    mx = DISC * (int(np.amax(data) / DISC) + 1)
    mn = DISC * (int(np.amin(data) / DISC) - 1)
    binsX = mod.assignbins([mn, mx], DISC)
    with contextlib.redirect_stdout(io.StringIO()):
        hist, _e, c1, c2, c3 = mod.reweight_CE(data, CUT, binsX, DISC, dV, T, False)
    c1, c2, c3 = -c1 / BETA, -c2 / BETA, -c3 / BETA
    pmf = mod.prephist(hist, T, EMAX)
    p3 = mod.normalize(pmf + c1 + c2 + c3, EMAX)
    valid = hist >= CUT
    return float(p3[valid].max() - p3[valid].min())


def load(name):
    p = ROOT / "data/per_frame" / f"{name}_frames.csv.gz"
    with gzip.open(p, "rt") as f:
        cols = f.readline().strip().split(",")
        arr = np.loadtxt(f)
    return arr[:, cols.index("cv1_prot_dna_com")], arr[:, cols.index("dV_D")]


def midrank(v):
    v = list(v)
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for t in range(i, j + 1):
            r[order[t]] = (i + j + 2) / 2.0
        i = j + 1
    return r


def rho_p(xs, ys):
    n = len(xs)
    rx, ry = midrank(xs), midrank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    sxx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    syy = math.sqrt(sum((b - my) ** 2 for b in ry))

    def r_of(seq):
        m2 = sum(seq) / n
        s2 = math.sqrt(sum((a - m2) ** 2 for a in seq))
        return sum((a - m2) * (b - my) for a, b in zip(seq, ry)) / (s2 * syy)

    obs = r_of(rx)
    cnt = 0
    for perm in itertools.permutations(xs):
        if max(perm) == min(perm):
            continue
        if abs(r_of(midrank(perm))) >= abs(obs) - 1e-12:
            cnt += 1
    return round(obs, 3), cnt


# ── A. archive-consistency checks ───────────────────────────────────────
print("── A. archive consistency (production 50-ps) ──")
A = [("orig_s2_APO", 19.5), ("sys2_rucaparib", 28.7)]
for name, pub in A:
    cv, dv = load(name)
    s = spans(cv, dv)
    print(f"{name:18s} n={len(cv):4d} span={s:6.2f} (published {pub:4.1f}) Δ={s-pub:+.2f}")
cv, dv = load("orig_s2_veliparib")
ws = [spans(cv[s:s + 40], dv[s:s + 40]) for s in range(0, len(cv) - 39, 40)]
print(f"veliparib 2-ns window SD (50-ps series, {len(ws)} windows): {np.std(ws):.2f}")

# ── B. phenotype-rank encoding audit over the lattice ───────────────────
print("\n── B. rank-encoding × lattice audit (S2 CV1 C3) ──")
# S2 CV1 C3 values (Table S5; production series)
V = {"tala": 22.9, "veli": 28.1,
     "nira": {0: 26.9, 1: 34.1}, "ruca": {0: 28.7, 1: 48.7}, "veliP": {0: 28.1, 1: 25.5}}
V = {"tala": 22.9, "ola": 28.8, "nira": {0: 26.9, 1: 34.1},
     "ruca": {0: 28.7, 1: 48.7}, "veli": {0: 28.1, 1: 25.5}}

ENCODINGS = {
    "main (tala>nira>ola~ruca>veli)": [5, 4, 2.5, 2.5, 1],
    "Pommier (nira~ola~ruca tied)": [5, 3, 3, 3, 1],
    "ola>ruca split": [5, 4, 3, 2, 1],
    "ruca>ola split": [5, 4, 2, 3, 1],
    "nira~ola tied": [5, 3.5, 3.5, 2, 1],
}
rows = []
print(f"{'encoding':38s} {'state':4s} {'rho':>7s} {'p':>7s}")
for ename, Y in ENCODINGS.items():
    for s_nira, s_ruca, s_veli in itertools.product((0, 1), repeat=3):
        xs = [V["tala"], V["nira"][s_nira], V["ola"], V["ruca"][s_ruca], V["veli"][s_veli]]
        r, cnt = rho_p(xs, Y)
        state = f"{s_nira}{s_ruca}{s_veli}"
        rows.append((ename, state, r, cnt / 120, cnt))
        print(f"{ename:38s} {state:4s} {r:+7.3f} {cnt/120:7.3f}")
r_all = [x[2] for x in rows]
p_all = [x[3] for x in rows]
print(f"\nρ across 40 encoding×lattice combinations: [{min(r_all):+.2f}, {max(r_all):+.2f}]; "
      f"min exact p = {min(p_all):.3f}")

# ── C. niraparib seed-summary sensitivity ───────────────────────────────
print("\n── C. niraparib seed summary sensitivity (main encoding) ──")
Y = [5, 4, 2.5, 2.5, 1]
variants = {
    "neutral: seed42 26.9 (published)": (26.9, 28.7, 28.1),
    "neutral: 3-seed mean 33.0": (33.0, 28.7, 28.1),
    "neutral: 3-seed median 26.9": (26.9, 28.7, 28.1),
    "audit: seed42 34.1 (published)": (34.1, 48.7, 25.5),
    "audit: 3-seed mean 55.2": (55.2, 48.7, 25.5),
    "audit: 3-seed median 49.3": (49.3, 48.7, 25.5),
}
srows = []
for name, (nira, ruca, veli) in variants.items():
    xs = [22.9, nira, 28.8, ruca, veli]
    r, cnt = rho_p(xs, Y)
    srows.append((name, r, cnt, cnt / 120))
    print(f"{name:38s} rho={r:+.3f}  p={cnt}/120={cnt/120:.3f}")

with open(ROOT / "data/analysis/r6/s60_rank_encodings.csv", "w") as f:
    f.write("encoding,state_bitmask,rho,p_exact,cnt120\n")
    for r in rows:
        f.write(",".join(str(x) for x in r) + "\n")
with open(ROOT / "data/analysis/r6/s60_seed_sensitivity.csv", "w") as f:
    f.write("variant,rho,p_exact,cnt120\n")
    for r in srows:
        f.write(",".join(str(x) for x in r) + "\n")
print("\nwrote s60 CSVs")
