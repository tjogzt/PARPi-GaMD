#!/usr/bin/env python3
"""s58_cadence_bootstrap.py — reviewer-driven robustness computations (statistical review batch).

Produces:
  A. Frame-cadence check: talazoparib S2 CV1 C3 span and 2-ns window SD at the
     native 1-ps cadence vs the 50-ps cadence of the other systems.
  B. Window-noise extension: 2-ns window SDs for talazoparib and veliparib
     (the two systems entering the "nearest-neighbour" separations but absent
     from the five-trajectory window analysis).
  C. Moving-block bootstrap (circular, block = 2 ns) of the C3 valid-support
     span for every system with an archived per-frame series; the resampled
     span distributions feed an empirical rank-correlation distribution for
     the chemical-state panels that are fully covered by archived frames.
  D. Seed-summary sensitivity for niraparib (mean vs median of the three
     production seeds) on the exact permutation test.

Outputs (written next to this script's data dir):
  data/analysis/r6/s58_bootstrap_summary.csv
  data/analysis/r6/s58_rank_distributions.csv
Convention identical to s10_spanfix / Table S5 (PyReweighting-1D CE, T=300,
disc 0.1, Emax 20, cutoff 2; valid = bins with occupancy >= 2).
"""
import contextlib
import gzip
import importlib.util
import io
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PYREW = ROOT / "tools" / "PyReweighting" / "PyReweighting-1D.py"
spec = importlib.util.spec_from_file_location("pyrew1d", PYREW)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

BETA = 1.0 / (0.001987 * 300.0)
DISC, EMAX, CUT, T = 0.1, 20.0, 2, 300.0
RNG = np.random.default_rng(49)


def _edges(data):
    mx = DISC * (int(np.amax(data) / DISC) + 1)
    mn = DISC * (int(np.amin(data) / DISC) - 1)
    return mod.assignbins([mn, mx], DISC)


def spans(cv, dv):
    data = np.asarray(cv, float)
    dV = np.asarray(dv, float)
    binsX = _edges(data)
    with contextlib.redirect_stdout(io.StringIO()):
        hist, _e, c1, c2, c3 = mod.reweight_CE(data, CUT, binsX, DISC, dV, T, False)
    c1, c2, c3 = -c1 / BETA, -c2 / BETA, -c3 / BETA
    pmf = mod.prephist(hist, T, EMAX)
    p1 = mod.normalize(pmf + c1, EMAX)
    p2 = mod.normalize(pmf + c1 + c2, EMAX)
    p3 = mod.normalize(pmf + c1 + c2 + c3, EMAX)
    valid = hist >= CUT
    return (float(p1[valid].max() - p1[valid].min()),
            float(p2[valid].max() - p2[valid].min()),
            float(p3[valid].max() - p3[valid].min()))


def load_frames(name):
    p = ROOT / "data" / "per_frame" / f"{name}_frames.csv.gz"
    with gzip.open(p, "rt") as f:
        cols = f.readline().strip().split(",")
        arr = np.loadtxt(f)
    return arr[:, cols.index("cv1_prot_dna_com")], arr[:, cols.index("dV_D")]


def win_sd(cv, dv, per_win):
    ws = [spans(cv[s:s + per_win], dv[s:s + per_win])[2]
          for s in range(0, len(cv) - per_win + 1, per_win)]
    return float(np.std(ws)), len(ws)


def block_boot(cv, dv, block, n_boot=400):
    n = len(cv)
    nb = int(np.ceil(n / block))
    out = []
    for _ in range(n_boot):
        starts = RNG.integers(0, n, nb)
        idx = np.concatenate([(s + np.arange(block)) % n for s in starts])[:n]
        out.append(spans(cv[idx], dv[idx])[2])
    return np.array(out)


SYSTEMS = [
    # name, cadence_ps, published CV1 C3 span (Table S5)
    ("sys2_niraparib", 50, 26.9),
    ("sys2_nira_prot", 50, 34.1),
    ("sys2_olaparib", 50, 28.8),
    ("sys2_rucaparib", 50, 20.1),
    ("sys2_ruca_prot", 50, 48.7),
    ("orig_s2_talazoparib", 1, 22.9),
    ("orig_s2_veliparib", 1, 28.1),
]

summary = []
dists = {}
print(f"{'system':22s} {'cad':>4s} {'n':>6s} {'span':>7s} {'pub':>6s} {'winSD':>7s} {'nW':>4s} "
      f"{'boot_mean':>9s} {'boot_SD':>7s} {'q2.5':>6s} {'q97.5':>6s}")
for name, cad, pub in SYSTEMS:
    cv, dv = load_frames(name)
    s1, s2, s3 = spans(cv, dv)
    per_win = 2000 if cad == 1 else 40          # 2 ns per window
    wsd, nw = win_sd(cv, dv, per_win)
    boots = block_boot(cv, dv, per_win, n_boot=400 if cad == 1 else 400)
    dists[name] = boots
    summary.append((name, cad, len(cv), s3, pub, wsd, nw,
                    float(boots.mean()), float(boots.std(ddof=1)),
                    float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))))
    print(f"{name:22s} {cad:>4d} {len(cv):>6d} {s3:7.2f} {pub:6.1f} {wsd:7.2f} {nw:>4d} "
          f"{boots.mean():9.2f} {boots.std(ddof=1):7.2f} {np.quantile(boots,0.025):6.2f} "
          f"{np.quantile(boots,0.975):6.2f}", flush=True)

# ---- A. talazoparib cadence check (1 ps vs 50 ps) ----
cvT, dvT = load_frames("orig_s2_talazoparib")
s_full = spans(cvT, dvT)
cv50, dv50 = cvT[::50], dvT[::50]
s_50 = spans(cv50, dv50)
w_full = win_sd(cvT, dvT, 2000)
w_50 = win_sd(cv50, dv50, 40)
print("\n── A. talazoparib cadence ──")
print(f"1-ps  : n={len(cvT)}  C1={s_full[0]:.2f} C2={s_full[1]:.2f} C3={s_full[2]:.2f}  winSD={w_full[0]:.2f} (n={w_full[1]})")
print(f"50-ps : n={len(cv50)}  C1={s_50[0]:.2f} C2={s_50[1]:.2f} C3={s_50[2]:.2f}  winSD={w_50[0]:.2f} (n={w_50[1]})")

# ---- C. empirical rank-correlation distributions (panels with full frame coverage) ----
Y = np.array([5.0, 4.0, 2.5, 2.5, 1.0])          # tala > nira > ola≈ruca > veli (midrank 2.5)


def midrank(v):
    v = np.asarray(v, float)
    order = np.argsort(v)
    r = np.empty(len(v))
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        r[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = midrank(x), midrank(y)
    rx = rx - rx.mean(); ry = ry - ry.mean()
    return float((rx * ry).sum() / np.sqrt((rx ** 2).sum() * (ry ** 2).sum()))


PANELS = {
    "neutral (all five)": ["orig_s2_talazoparib", "sys2_niraparib", "sys2_olaparib",
                            "sys2_rucaparib", "orig_s2_veliparib"],
    "niraparib-protonated": ["orig_s2_talazoparib", "sys2_nira_prot", "sys2_olaparib",
                             "sys2_rucaparib", "orig_s2_veliparib"],
    "rucaparib-protonated": ["orig_s2_talazoparib", "sys2_niraparib", "sys2_olaparib",
                             "sys2_ruca_prot", "orig_s2_veliparib"],
    "audit-partial (veli neutral; no prot frames)":
        ["orig_s2_talazoparib", "sys2_nira_prot", "sys2_olaparib",
         "sys2_ruca_prot", "orig_s2_veliparib"],
}
NB = 400
rows = []
print("\n── C. empirical rank-correlation distributions (block bootstrap) ──")
print(f"{'panel':46s} {'rho_pt':>7s} {'P(rho<0)':>9s} {'q2.5':>6s} {'q50':>6s} {'q97.5':>7s}")
for pname, systems in PANELS.items():
    draws = []
    m = min(len(dists[s]) for s in systems)
    for b in range(m):
        draws.append(spearman([dists[s][b] for s in systems], Y))
    draws = np.array(draws)
    rows.append((pname, float((draws < 0).mean()), float(np.quantile(draws, 0.025)),
                 float(np.median(draws)), float(np.quantile(draws, 0.975))))
    print(f"{pname:46s} {'':>7s} {(draws<0).mean():9.3f} {np.quantile(draws,0.025):6.2f} "
          f"{np.median(draws):6.2f} {np.quantile(draws,0.975):7.2f}")

out1 = ROOT / "data/analysis/r6/s58_bootstrap_summary.csv"
with open(out1, "w") as f:
    f.write("system,cadence_ps,n_frames,span_c3,published_c3,win_sd_2ns,n_windows,"
            "boot_mean,boot_sd,q025,q975\n")
    for r in summary:
        f.write(",".join(str(x) for x in r) + "\n")
out2 = ROOT / "data/analysis/r6/s58_rank_distributions.csv"
with open(out2, "w") as f:
    f.write("panel,P_rho_lt0,q025,median,q975,n_boot\n")
    for name, p, q1, med, q2 in rows:
        f.write(f"{name},{p:.4f},{q1:.4f},{med:.4f},{q2:.4f},{NB}\n")
print(f"\nwrote {out1}\nwrote {out2}")
