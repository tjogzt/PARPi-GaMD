#!/usr/bin/env python3
"""S2 window-timescale diagnostic (R1 W3 closure): C3 span window-to-window
SD versus window length on the six rerun / four protonated S2 trajectories
(COM CV1, DBE reweighting), extending diagnostic (3) from S1 to S2.
Output: results/analysis/s2_window_timescale.csv

Purpose:  S2 window-timescale diagnostic (diagnostic 3 extended to S2): C3 span window-to-window SD versus window length on the rerun/protonated S2 trajectories; feeds the sampling-toolkit figure.
Created:  packaged 2026-10-06
Depends:  common.pmf, numpy, gzip
Run:      python3 scripts/s2_window_timescale.py   (from the repository root)
"""
import gzip
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis"
WINDOWS_NS = [2, 4, 6, 12, 24]
FPS = 20  # frames per ns (50 ps frames)

tags = ["sys2_niraparib", "sys2_olaparib", "sys2_rucaparib",
        "sys2_nira_prot", "sys2_ruca_prot"]

rows = []
print(f"{'tag':20s} {'win_ns':>6s} {'n_win':>6s} {'C3_sd':>7s} {'C3_mean':>8s}",
      flush=True)
for tag in tags:
    with gzip.open(f"results/archive/per_frame/{tag}_frames.csv.gz", "rt") as f:
        f.readline()
        arr = np.loadtxt(f)
    dv = arr[:, 2]
    cv = arr[:, 3]
    n = len(cv)
    for wns in WINDOWS_NS:
        wf = int(wns * FPS)
        if wf > n:
            continue
        spans = []
        for start in range(0, n - wf + 1, wf):
            seg_cv = cv[start:start + wf]
            seg_dv = dv[start:start + wf]
            w = np.column_stack([BETA * seg_dv, np.zeros(wf), seg_dv])
            r = run_pyrew(seg_cv, w)
            spans.append(r.c3)
        spans = np.array(spans)
        rows.append((tag, wns, len(spans), spans.std(), spans.mean()))
        print(f"{tag:20s} {wns:6d} {len(spans):6d} {spans.std():7.1f} "
              f"{spans.mean():8.1f}", flush=True)

with open(OUT / "s2_window_timescale.csv", "w") as f:
    f.write("tag,window_ns,n_windows,c3_span_sd,c3_span_mean\n")
    for r in rows:
        f.write(",".join(map(str, r)) + "\n")
print("WROTE", OUT / "s2_window_timescale.csv", flush=True)
