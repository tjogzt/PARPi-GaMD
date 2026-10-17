#!/usr/bin/env python3
"""s52_positive_control.py — positive control for the six-diagnostic audit.

Purpose:  Planted-association recovery benchmark for the audit's own test,
          with the exact tie structure of the analyzed panel. Synthetic n = 5
          panels use the reference trapping rank vector
          y = (5, 4, 2.5, 2.5, 1)  (talazoparib > niraparib >
          olaparib = rucaparib > veliparib; midrank 2.5 for the
          olaparib--rucaparib tie, weakest-to-strongest coding), latent spans
          y * delta, and per-system Gaussian noise sigma. Detection = tie-aware
          Spearman rho between the noisy observations and y, with the exact
          two-sided permutation test over all 5! = 120 labelings of y (60
          distinct tied rank vectors; minimum attainable p = 4/120 = 0.0333,
          reached only by an essentially perfect monotone ordering) at the
          p < 0.05 threshold. SIGMA grid = {2, 4.6, 8, 15, 18.4} kcal/mol:
          4.6-18.4 spans the measured 2-ns window variability of the five
          re-simulated S2 trajectories; 2 is deliberately below that range.
          DELTA grid = 0-48 kcal/mol in steps of 2; the delta = 0 rows
          calibrate the nominal false-positive level of the test (no planted
          signal). Seed 49; 1,500 worlds per cell; vectorized.
Inputs:   none (synthetic)
Outputs:  data/analysis/r6/c7g_positive_control.csv
Depends:  numpy
Run:      python3 scripts/s52_positive_control.py   (from the repository root)
"""
import csv
import itertools
import math

import numpy as np

rng = np.random.default_rng(49)

# reference trapping ranks: talazoparib, niraparib, olaparib, rucaparib, veliparib
Y = np.array([5.0, 4.0, 2.5, 2.5, 1.0])
DELTAS = np.arange(0.0, 50.0, 2.0)
SIGMAS = [2.0, 4.6, 8.0, 15.0, 18.4]
N_WORLDS = 1500
N_PERM = math.factorial(5)


def midranks(v):
    """Midranks (average ranks) of a 1-D array."""
    order = np.argsort(v, kind="mergesort")
    sv = v[order]
    r = np.empty(len(v))
    j = 0
    m = len(v)
    while j < m:
        k = j
        while k + 1 < m and sv[k + 1] == sv[j]:
            k += 1
        r[order[j:k + 1]] = (j + k) / 2.0 + 1.0
        j = k + 1
    return r


def midranks_matrix(X):
    """Row-wise midranks for an (n_worlds, 5) array."""
    n, _ = X.shape
    R = np.empty_like(X)
    for i in range(n):
        R[i] = midranks(X[i])
    return R


# exact permutation null: all 5! labelings of the reference rank vector y
# (the duplicated 2.5 gives 60 distinct tied rank vectors over 120 labelings)
Y_PERMS = np.array([Y[list(p)] for p in itertools.permutations(range(5))], float)
RY_PERMS = midranks_matrix(Y_PERMS)
RY_C = RY_PERMS - RY_PERMS.mean(axis=1, keepdims=True)
RY_NORM = np.sqrt((RY_C ** 2).sum(axis=1))

ry_obs = midranks(Y)
ry_obs_c = ry_obs - ry_obs.mean()
ry_obs_norm = np.sqrt((ry_obs_c ** 2).sum())


def detect_rate(delta, sigma):
    latent = Y * float(delta)
    obs = latent[None, :] + rng.normal(0.0, sigma, (N_WORLDS, 5))
    rx = midranks_matrix(obs)
    rx_c = rx - rx.mean(axis=1, keepdims=True)
    rx_norm = np.sqrt((rx_c ** 2).sum(axis=1))
    rho_perm = (rx_c @ RY_C.T) / (rx_norm[:, None] * RY_NORM[None, :])  # (W, 120)
    rho_obs = (rx_c @ ry_obs_c) / (rx_norm * ry_obs_norm)               # (W,)
    cnt = (np.abs(rho_perm) >= np.abs(rho_obs)[:, None] - 1e-12).sum(axis=1)
    return float((cnt / N_PERM < 0.05).mean())


rows = []
for sigma in SIGMAS:
    for delta in DELTAS:
        p = detect_rate(delta, sigma)
        mcse = math.sqrt(max(p * (1.0 - p), 0.0) / N_WORLDS)
        rows.append(dict(delta_kcal=float(delta), sigma_kcal=float(sigma),
                         detect_rate=p, n_worlds=N_WORLDS, mcse=mcse))
        print(f"  sigma={sigma:>4g} delta={delta:>4g}  rate={p:.3f}  mcse={mcse:.3f}", flush=True)

out = "data/analysis/r6/c7g_positive_control.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["delta_kcal", "sigma_kcal", "detect_rate", "n_worlds", "mcse"])
    w.writeheader()
    w.writerows(rows)
print("saved", out, len(rows), "cells (", len(SIGMAS), "sigmas x", len(DELTAS), "deltas )")
