#!/usr/bin/env python3
#
# Purpose:  Average the three replicate PMF curves (s7/s49/s123, 24 ns truncated, DBE)
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s28_replicate_mean_pmf.py   (from the repository root)

"""Average the three replicate PMF curves (s7/s49/s123, 24 ns truncated, DBE)
for talazoparib and veliparib onto a common CV grid, and write the mean C3 PMF
xvg into results/analysis/ so figure scripts use the 3-seed average that
matches Table S1 (46.2 / 61.1 kcal/mol)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
DATA = data_root()
REPL = Path(__file__).resolve().parents[1] / "results" / "analysis" / "replicates"
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis"


def w_dbe(dv):
    return np.column_stack([BETA * dv, np.zeros(len(dv)), dv])


def pmf_curve(tmpdir):
    lines = [l for l in (tmpdir / "pmf-c3-cv.dat.xvg").read_text().splitlines()
             if not l.startswith(("#", "@")) and l.strip()]
    return np.array([[float(x) for x in l.split()] for l in lines])


def main():
    for lig in ["talazoparib", "veliparib"]:
        curves = []
        for s in ["s7", "s49", "s123"]:
            d = DATA / "review_runs" / f"sys1_{lig}_{s}"
            log = np.loadtxt(d / "gamd.log", comments="#")
            prod = log[:, 1] >= 6.0e6
            dv = log[:, 7][prod][:480]
            cv_full = np.loadtxt(REPL / f"sys1_{lig}_{s}" / "cv.dat")
            cv = cv_full[prod][:480]
            n = min(len(cv), len(dv))
            r = run_pyrew(cv[:n], w_dbe(dv[:n]))
            curves.append(pmf_curve(r.tmpdir))
            print(f"{lig}_{s}: C3={r.c3:.1f}")
        # common grid
        lo = min(c[:, 0].min() for c in curves)
        hi = max(c[:, 0].max() for c in curves)
        grid = np.arange(lo, hi + 0.05, 0.05)
        ys = np.array([np.interp(grid, c[:, 0], c[:, 1]) for c in curves])
        ymean = ys.mean(axis=0)
        np.savetxt(OUT / f"sys1_{lig}_pmf_c3.xvg",
                   np.column_stack([grid, ymean]), fmt="%.6f %.4f")
        span = ymean.max() - ymean.min()
        print(f"{lig}: mean-curve span={span:.1f} -> sys1_{lig}_pmf_c3.xvg")
    print("DONE")


if __name__ == "__main__":
    main()
