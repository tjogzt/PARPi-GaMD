#!/usr/bin/env python3
"""s52_positive_control.py — positive control for the six-diagnostic audit.

Purpose:  Planted-association recovery benchmark: synthetic n = 5 panels with a
          known rank association (latent spans 1..5 x delta_kcal) plus Gaussian
          trajectory noise sigma; detection = tie-aware Spearman with exact
          two-tailed permutation p over all 120 orderings (mirroring the audit's
          own test), significance threshold p < 0.05. Sigma levels bracket the
          measured 2-ns window noise (4.6-18.4 kcal/mol); delta sweeps
          2-48 kcal/mol. Seed 49; 1,500 worlds per cell; vectorized.
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

N_RANK = 5
DELTAS = np.arange(2.0, 50.0, 2.0)
SIGMAS = [2.0, 8.0, 15.0]
N_WORLDS = 1500
N_PERM = math.factorial(N_RANK)

perms = np.array(list(itertools.permutations(range(1, N_RANK + 1))), float)  # (120,5)
P_c = perms - perms.mean(axis=1, keepdims=True)
P_norm = np.sqrt((P_c ** 2).sum(axis=1))  # all equal sqrt(10)


def midranks_matrix(X):
    """Row-wise midranks for (n_worlds, 5) array."""
    n, m = X.shape
    R = np.empty_like(X)
    for i in range(n):
        v = X[i]
        order = np.argsort(v, kind="mergesort")
        sv = v[order]
        r = np.empty(m)
        j = 0
        while j < m:
            k = j
            while k + 1 < m and sv[k + 1] == sv[j]:
                k += 1
            r[order[j:k + 1]] = (j + k) / 2.0 + 1.0
            j = k + 1
        R[i] = r
    return R


def detect_rate(latent, sigma):
    obs = latent[None, :] + rng.normal(0.0, sigma, (N_WORLDS, N_RANK))
    Rx = midranks_matrix(obs)
    Rx_c = Rx - Rx.mean(axis=1, keepdims=True)
    Rx_norm = np.sqrt((Rx_c ** 2).sum(axis=1))
    rho_all = (Rx_c @ P_c.T) / (Rx_norm[:, None] * P_norm[None, :])  # (W,120)
    # observed rho: compared against each relabeling of the truth (identity included)
    rho_obs = rho_all[:, 0]  # identity permutation is the first row of perms (1,2,3,4,5)
    cnt = (np.abs(rho_all) >= np.abs(rho_obs)[:, None] - 1e-12).sum(axis=1)
    return float((cnt < 0.05 * N_PERM).mean())


rows = []
for sigma in SIGMAS:
    for delta in DELTAS:
        latent = np.arange(1, N_RANK + 1, dtype=float) * float(delta)
        rows.append(dict(delta_kcal=float(delta), sigma_kcal=float(sigma),
                         detect_rate=detect_rate(latent, sigma), n_worlds=N_WORLDS))
        print(f"  sigma={sigma:>4g} delta={delta:>4g}  rate={rows[-1]['detect_rate']:.3f}", flush=True)

out = "data/analysis/r6/c7g_positive_control.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["delta_kcal", "sigma_kcal", "detect_rate", "n_worlds"])
    w.writeheader()
    w.writerows(rows)
print("saved", out, len(rows), "cells")
