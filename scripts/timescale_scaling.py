#!/usr/bin/env python3
"""Convergence-timescale scaling: window-SD of the S1 span vs window length,
using the existing tala/veli 200-400 ns replicate trajectories.
Answers: marginal value of longer sampling, directly measured.

Purpose:  Convergence-timescale scaling (diagnostic 3): window-SD of the S1 span versus window length on the talazoparib/veliparib replicate trajectories; feeds the D4 timescale figure (code/33).
Created:  packaged 2026-10-06
Depends:  common.paths, numpy
Run:      python3 scripts/timescale_scaling.py   (from the repository root)
"""
import numpy as np
from pathlib import Path
import sys

from common.paths import data_root

DATA_ROOT = data_root()
REPL = DATA_ROOT / "review_runs"
OUT = Path("results/analysis")

def load_prod(tag, seed):
    d = REPL / f"sys1_{tag}_{seed}"
    lg = np.loadtxt(d / "gamd.log", comments="#")
    prod = lg[:, 1] >= 6.0e6
    dv = lg[:, 7][prod]
    cv_full = np.loadtxt(Path("results/analysis/replicates") / f"sys1_{tag}_{seed}" / "cv.dat")
    cv = cv_full[prod]
    n = min(len(cv), len(dv))
    return cv[:n], dv[:n]

def reweighted_pmf_span(cv, dv, nbin=40):
    lo, hi = cv.min(), cv.max()
    edges = np.linspace(lo, hi, nbin + 1)
    idx = np.clip(np.digitize(cv, edges) - 1, 0, nbin - 1)
    beta = 1.0 / 0.5961  # 300 K in kcal/mol
    w = np.exp(beta * dv)
    w /= w.sum()
    h = np.bincount(idx, weights=w, minlength=nbin)
    h = h / h.sum()
    nz = h > 0
    pmf = -np.log(h[nz] + 1e-30)
    pmf -= pmf.min()
    return pmf.max()

# window-length scaling: split production into blocks of length L ns, span SD across blocks
results = {}
for tag, seeds in [("talazoparib", ["s7", "s49", "s123"]),
                   ("veliparib", ["s7", "s49", "s123"])]:
    print(f"== {tag} ==", flush=True)
    for seed in seeds:
        cv, dv = load_prod(tag, seed)
        n = len(cv)
        frame_ps = 50.0
        total_ns = n * frame_ps / 1000
        print(f"  {seed}: n={n} ({total_ns:.1f} ns production)", flush=True)
        for L_ns in [2, 4, 6, 12, 24]:
            blk = int(L_ns * 1000 / frame_ps)
            if blk > n:
                continue
            spans = []
            for i in range(0, n - blk + 1, blk):
                spans.append(reweighted_pmf_span(cv[i:i+blk], dv[i:i+blk]))
            spans = np.array(spans)
            results[(tag, seed, L_ns)] = (spans.mean(), spans.std(ddof=1), len(spans))
            print(f"    {L_ns:>3}ns x{len(spans)} blocks: span {spans.mean():.1f} +/- {spans.std(ddof=1):.1f}",
                  flush=True)

with open(OUT / "timescale_scaling.csv", "w") as f:
    f.write("ligand,seed,window_ns,span_mean,span_sd,nblocks\n")
    for (tag, seed, L), (m, s, k) in sorted(results.items()):
        f.write(f"{tag},{seed},{L},{m:.3f},{s:.3f},{k}\n")
print("WROTE results/analysis/timescale_scaling.csv")
