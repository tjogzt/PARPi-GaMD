#!/usr/bin/env python3
#
# Purpose:  D4-1: Bayesian ordinal model + LOO/counterfactual fragility analysis.
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s35_bayesian_ordinal.py   (from the repository root)

"""D4-1: Bayesian ordinal model + LOO/counterfactual fragility analysis.

Bayesian latent model: true spans theta_i ~ Normal(x_i, sd_i), rho = cor(rank(theta), y).
MCMC (random-walk Metropolis on theta) -> posterior P(rho < 0) for the three
associations (S2 CV1 primary, S1, AAI), n = 5.

Counterfactuals:
  (a) add AZD5305 as a 6th system with assumed trapping rank (worst case):
      rho/p shift of the primary association.
  (b) flip grid: move veliparib's S2 CV1 span by d in [-8, +8] kcal/mol;
      find the flip threshold where rho becomes non-negative.
  (c) LOO on all 5 systems (report min/max rho).
"""
import numpy as np

rng = np.random.default_rng(42)
ks = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
trap_rank = {"talazoparib": 5.0, "niraparib": 4.0, "olaparib": 2.5,
             "rucaparib": 2.5, "veliparib": 1.0}

s2cv1 = {"talazoparib": 53.1, "niraparib": 62.0, "olaparib": 58.5,
         "rucaparib": 59.3, "veliparib": 58.9}
s1 = {"talazoparib": 51.4, "niraparib": 59.3, "olaparib": 51.6,
      "rucaparib": 48.7, "veliparib": 62.9}
aai = {"talazoparib": 0.830, "niraparib": 1.174, "olaparib": 1.075,
       "rucaparib": 0.926, "veliparib": 0.863}
# per-system uncertainty: half the C1-C3 SD for S2 (convergence diag window SD
# gives ~8-16 at 2ns; use per-system 24ns window SD estimates ~3-4). Conservative:
sd_s2 = 4.0
sd_s1 = {"talazoparib": 7.6, "niraparib": 13.4, "olaparib": 8.9,
         "rucaparib": 7.8, "veliparib": 4.8}
sd_aai = {"talazoparib": 0.15, "niraparib": 0.10, "olaparib": 0.19,
          "rucaparib": 0.15, "veliparib": 0.21}

# rucaparib-protonated panel (upper-bound perturbation of the ~9% minority state)
s2cv1_ruca = {**s2cv1, "rucaparib": 79.7}
s1_ruca = {**s1, "rucaparib": 61.0}
aai_ruca = {**aai, "rucaparib": 61.0 / 75.4}
sd_s2_ruca = 4.0
sd_s1_ruca = {**sd_s1, "rucaparib": 14.2}
sd_aai_ruca = {**sd_aai, "rucaparib": 0.20}


def spearman(x, y):
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return np.corrcoef(rx, ry)[0, 1]


def bayesian_p_rho_neg(x_obs, sd, y, n_iter=60000, burn=5000, prior="flat"):
    """Posterior P(rho < 0) for latent ordinal model.

    prior="flat":   flat prior on latent theta (inducing approx. uniform rho)
    prior="arcsine": Jeffrey's arcsine prior on rho, density ~ 1/sqrt(1-rho^2),
                    applied in the acceptance step via log-prior difference.
    """
    theta = np.array([x_obs[k] for k in ks], dtype=float)
    yv = np.array([y[k] for k in ks])
    sdv = np.array([sd[k] if isinstance(sd, dict) else sd for k in ks], dtype=float)
    logpost = lambda t: -0.5 * np.sum(((t - np.array([x_obs[k] for k in ks])) / sdv) ** 2)

    def log_prior(r):
        if prior == "arcsine":
            return -0.5 * np.log(max(1e-6, 1.0 - r * r))
        return 0.0

    rho_cur = spearman(theta, yv)
    acc = 0
    samples = []
    for i in range(n_iter):
        prop = theta + rng.normal(0, sdv * 0.5)
        dlog = logpost(prop) - logpost(theta)
        rho_prop = spearman(prop, yv)
        dlog += log_prior(rho_prop) - log_prior(rho_cur)
        if np.log(rng.random()) < dlog:
            theta, rho_cur = prop, rho_prop
            acc += 1
        if i >= burn and i % 5 == 0:
            samples.append(rho_cur)
    samples = np.array(samples)
    return samples, acc / n_iter


out = []
for name, x_obs, sd in [("S2_CV1", s2cv1, sd_s2), ("S1", s1, sd_s1), ("AAI", aai, sd_aai)]:
    samples, acc = bayesian_p_rho_neg(x_obs, sd, trap_rank)
    p_neg = (samples < 0).mean()
    q = np.percentile(samples, [2.5, 50, 97.5])
    out.append(f"{name}: posterior P(rho<0)={p_neg:.3f}, rho median={q[1]:+.3f} [95% CI {q[0]:+.3f}, {q[2]:+.3f}], acc={acc:.2f}")
    print(out[-1], flush=True)

print("\nRucaparib-protonated panel:")
for name, x_obs, sd in [("S2_CV1", s2cv1_ruca, sd_s2_ruca), ("S1", s1_ruca, sd_s1_ruca), ("AAI", aai_ruca, sd_aai_ruca)]:
    samples, acc = bayesian_p_rho_neg(x_obs, sd, trap_rank)
    p_neg = (samples < 0).mean()
    q = np.percentile(samples, [2.5, 50, 97.5])
    out.append(f"{name}_ruca_prot: posterior P(rho<0)={p_neg:.3f}, rho median={q[1]:+.3f} [95% CI {q[0]:+.3f}, {q[2]:+.3f}], acc={acc:.2f}")
    print(out[-1], flush=True)

# (a) add AZD5305 (S2 CV1 = 79.4; assumed worst trapping rank 0.5 of 6 -> rank 0.5)
print("\nCounterfactual (a): n=6 with AZD5305 (span 79.4, assumed rank 0.5 of 6)")
x6 = dict(s2cv1); x6["AZD5305"] = 79.4
y6 = dict(trap_rank); y6["AZD5305"] = 0.5
k6 = list(x6.keys())
vals = np.array([x6[k] for k in k6]); ranks = np.array([y6[k] for k in k6])
rho6 = spearman(vals, ranks)
perms = np.array([p for p in __import__("itertools").permutations(range(6))])
from collections import Counter
# exact permutation p via rank-permutation of the 6 systems
rho_perm = np.array([np.corrcoef(np.argsort(np.argsort(np.array([vals[i] for i in p]))), np.argsort(np.argsort(ranks)))[0,1]
                     for p in perms])
p6 = (np.abs(rho_perm) >= abs(rho6)).mean()
print(f"n=6 (incl. AZD): rho={rho6:+.3f}, exact p={p6:.4f}")

# (b) flip grid on veliparib S2 CV1
print("\nCounterfactual (b): veliparib S2 CV1 flip grid")
base = np.array([s2cv1[k] for k in ks]); yv = np.array([trap_rank[k] for k in ks])
for d in [-8, -6, -4, -2, 0, 2, 4, 6, 8]:
    v = base.copy(); v[ks.index("veliparib")] += d
    # exact permutation p (n=5, 120 perms)
    perms5 = np.array([p for p in __import__("itertools").permutations(range(5))])
    rho = spearman(v, yv)
    rho_p = np.array([np.corrcoef(np.argsort(np.argsort(v[p])), np.argsort(np.argsort(yv)))[0, 1]
                      for p in perms5])
    p = (np.abs(rho_p) >= abs(rho)).mean()
    print(f"  veli span {s2cv1['veliparib']+d:5.1f}: rho={rho:+.3f}, p={p:.3f}")

# (c) LOO
print("\nCounterfactual (c): LOO (n=5, primary)")
for k in ks:
    sub = {kk: vv for kk, vv in s2cv1.items() if kk != k}
    keys = list(sub.keys())
    v = np.array([sub[kk] for kk in keys]); yv2 = np.array([trap_rank[kk] for kk in keys])
    print(f"  excl {k:12s}: rho={spearman(v, yv2):+.3f}")

with open("results/analysis/d4_bayesian_ordinal.txt", "w") as f:
    f.write("\n".join(out) + "\n")
print("\nWROTE results/analysis/d4_bayesian_ordinal.txt")
