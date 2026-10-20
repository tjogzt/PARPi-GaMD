#!/usr/bin/env python3
"""s64_positive_control_extra.py — extra grid points at the actual panel
separations (2.6, 3.3, 5.3 kcal/mol) for the positive-control benchmark.

Independent RNG stream (seed 1049) and 3000 worlds/cell; the original 2-kcal
grid in c7g_positive_control.csv is left untouched (its values are asserted in
the manuscript and audit). Output: c7g_positive_control_extra.csv
"""
import csv
import itertools
import math

import numpy as np

rng = np.random.default_rng(1049)
Y = np.array([5.0, 4.0, 2.5, 2.5, 1.0])
DELTAS = [2.6, 3.3, 5.3]
SIGMAS = [2.0, 4.6, 8.0, 15.0, 18.4]
N_WORLDS = 3000
N_PERM = math.factorial(5)


def midranks(v):
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
    n, _ = X.shape
    R = np.empty_like(X)
    for i in range(n):
        R[i] = midranks(X[i])
    return R


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
    rho_perm = (rx_c @ RY_C.T) / (rx_norm[:, None] * RY_NORM[None, :])
    rho_obs = (rx_c @ ry_obs_c) / (rx_norm * ry_obs_norm)
    cnt = (np.abs(rho_perm) >= np.abs(rho_obs)[:, None] - 1e-12).sum(axis=1)
    return float((cnt / N_PERM < 0.05).mean())


rows = []
for sigma in SIGMAS:
    for delta in DELTAS:
        p = detect_rate(delta, sigma)
        mcse = math.sqrt(max(p * (1.0 - p), 0.0) / N_WORLDS)
        rows.append(dict(delta_kcal=delta, sigma_kcal=sigma, detect_rate=p,
                         n_worlds=N_WORLDS, mcse=mcse))
        print(f"sigma={sigma:>4g} delta={delta:>4g}  rate={p:.3f}  mcse={mcse:.3f}", flush=True)

out = "data/analysis/r6/c7g_positive_control_extra.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["delta_kcal", "sigma_kcal", "detect_rate",
                                      "n_worlds", "mcse"])
    w.writeheader()
    w.writerows(rows)
print("saved", out)
