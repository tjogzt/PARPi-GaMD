#!/usr/bin/env python3
#
# Purpose:  Niraparib-protonated panel: veliparib single-point perturbation grid + q90 LOO + exact permutation p.
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s29_counterfactual_grid.py   (from the repository root)

"""Niraparib-protonated panel: veliparib single-point perturbation grid + q90 LOO + exact permutation p."""
import itertools
import math
import numpy as np

spans = {  # S2 CV1 C3, corrected panel
    "tala": 53.1, "nira": 62.0, "ola": 58.5, "ruca": 59.3, "veli": 58.9,
}
ranks = {"tala": 5, "nira": 4, "ola": 2.5, "ruca": 2.5, "veli": 1}


def midrank(v):
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


def rho(xs, ys):
    rx = midrank(np.array(xs, float))
    ry = midrank(np.array(ys, float))
    if rx.std() == 0 or ry.std() == 0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def exact_p(xs, ys):
    """tie-aware: enumerate all permutations of spans against fixed ranks,
    count |rho_perm| >= |rho_obs|."""
    obs = rho(xs, ys)
    cnt, tot, mx = 0, 0, 0.0
    for perm in itertools.permutations(xs):
        tot += 1
        r = rho(list(perm), ys)
        mx = max(mx, abs(r))
        if abs(r) >= abs(obs) - 1e-12:
            cnt += 1
    return obs, cnt / tot


names = list(spans)
base_x = [spans[k] for k in names]
rank_y = [ranks[k] for k in names]
r0, p0 = exact_p(base_x, rank_y)
print(f"base: rho={r0:.4f} p={p0:.4f}")

# veliparib perturbation grid
for d in [-8, -6, -4, -2, 2, 4, 6, 8]:
    xs = [spans[k] + (d if k == "veli" else 0) for k in names]
    r, p = exact_p(xs, rank_y)
    print(f"veli {d:+d}: rho={r:+.4f} p={p:.4f}")

# LOO on q90 estimator (n=5, values from matrix)
q90 = {"tala": 52.9, "nira": 54.6, "ola": 53.2, "ruca": 52.9, "veli": 55.8}
# placeholders from new matrix -- recompute below from real q90 spans
print("---")
# leave-one-out for base panel (min-max)
for k in names:
    xs = [spans[m] for m in names if m != k]
    ys = [ranks[m] for m in names if m != k]
    # tie-aware permutation for n=4
    obs = rho(xs, ys)
    cnt, tot = 0, 0
    for perm in itertools.permutations(xs):
        tot += 1
        if abs(rho(list(perm), ys)) >= abs(obs) - 1e-12:
            cnt += 1
    print(f"LOO excl {k}: rho={obs:+.4f} p={cnt/tot:.4f}")
