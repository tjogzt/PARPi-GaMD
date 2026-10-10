#!/usr/bin/env python3
"""C7 toy benchmark — window-noise behavior and span-estimator tail sensitivity on
synthetic processes with known ground truth.

(A) Window-SD vs window length for OU processes: fast (tau = 0.5 ns) vs slow (tau = 50 ns),
    vs the naive L^-1/2 reference. 20 realizations x 400 ns @ 50 ps.
(B) Span-estimator variability on a rare-excursion process (OU tau = 2 ns + sparse
    1-ns +3 sigma excursions, rate 0.02/ns): min-max span vs 5-95% quantile span,
    200 replicates x 26 ns.

Outputs: results/analysis/r6/c7_toy_benchmark.csv (curves) + c7_estimator_stats.csv
"""
import numpy as np
from pathlib import Path

OUT = Path("/Users/taozhu/clacky_workspace/PARPi_design/results/analysis/r6")
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(49)

dt = 0.05          # ns per frame (50 ps)
T = 400.0          # ns per realization
n = int(T / dt)
L_LIST = [2, 4, 6, 12, 24, 48, 96]   # window lengths (ns)
N_REAL = 20

def ou_series(tau, sigma, n, dt, rng):
    x = np.empty(n)
    x[0] = rng.normal(0, sigma)
    a = dt / tau
    s = sigma * np.sqrt(2 * a)
    noise = rng.normal(0, 1, n)
    for i in range(1, n):
        x[i] = x[i - 1] - a * x[i - 1] + s * noise[i]
    return x

rows = []
for name, tau in [("fast_tau0.5ns", 0.5), ("medium_tau8ns", 8.0), ("slow_tau50ns", 50.0)]:
    print(f"== {name} ==", flush=True)
    for L in L_LIST:
        blk = int(L / dt)
        sds = []
        for r in range(N_REAL):
            x = ou_series(tau, 1.0, n, dt, rng)
            ns = n // blk
            spans = np.array([x[i * blk:(i + 1) * blk].max() - x[i * blk:(i + 1) * blk].min()
                              for i in range(ns)])
            sds.append(spans.std(ddof=1))
        m, s = np.mean(sds), np.std(sds, ddof=1)
        rows.append({"process": name, "L_ns": L, "window_sd_mean": round(m, 4), "window_sd_se": round(s, 4)})
        print(f"  L={L:3d} ns: window-SD = {m:.3f} +/- {s:.3f} (sigma units)", flush=True)

# naive 1/sqrt(L) reference normalized at L=2
ref = np.array([np.sqrt(2.0 / L) for L in L_LIST])

import csv
with open(OUT / "c7_toy_benchmark.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=["process", "L_ns", "window_sd_mean", "window_sd_se"])
    wr.writeheader()
    wr.writerows(rows)
    for L, v in zip(L_LIST, ref):
        wr.writerow({"process": "naive_L^-1/2_ref", "L_ns": L, "window_sd_mean": round(float(v), 4), "window_sd_se": 0.0})

# ---- Part B: estimator tail sensitivity (inflation of span by rare short spikes) ----
print("\n== rare-spike estimator inflation ==", flush=True)
T2 = 26.0
n2 = int(T2 / dt)
N_REP = 200
mm_base, q90_base, mm_spk, q90_spk = [], [], [], []
for r in range(N_REP):
    x = ou_series(2.0, 1.0, n2, dt, rng)
    xb = x.copy()
    nev = rng.poisson(0.02 * T2)
    for _ in range(nev):
        t0 = rng.uniform(0, T2 - 0.5)
        i0 = int(t0 / dt); i1 = min(int((t0 + 0.4) / dt), n2)
        x[i0:i1] += 3.0   # short spikes (0.4 ns <= 1.6% of frames: below the 2.5% trim)
    mm_base.append(xb.max() - xb.min()); q90_base.append(np.percentile(xb, 95) - np.percentile(xb, 5))
    mm_spk.append(x.max() - x.min());    q90_spk.append(np.percentile(x, 95) - np.percentile(x, 5))
mm_base, q90_base = np.array(mm_base), np.array(q90_base)
mm_spk, q90_spk = np.array(mm_spk), np.array(q90_spk)

stats = []
for nm, b, s in [("minmax_span", mm_base, mm_spk), ("q90_span", q90_base, q90_spk)]:
    infl = s.mean() / b.mean() - 1.0
    stats.append({"estimator": nm, "mean_base": round(float(b.mean()), 4),
                  "mean_spiked": round(float(s.mean()), 4), "inflation": round(float(infl), 4),
                  "cv_spiked": round(float(s.std(ddof=1) / s.mean()), 4)})
    print(f"  {nm}: base={b.mean():.3f} spiked={s.mean():.3f} inflation={100*infl:+.1f}% CV={s.std(ddof=1)/s.mean():.3f}", flush=True)

with open(OUT / "c7_estimator_stats.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=["estimator", "mean_base", "mean_spiked", "inflation", "cv_spiked"])
    wr.writeheader()
    wr.writerows(stats)

ir = stats[0]["inflation"] / max(stats[1]["inflation"], 1e-9)
print(f"\ninflation ratio (min-max / q90) = {ir:.1f}x")
print("WROTE c7_toy_benchmark.csv + c7_estimator_stats.csv")
