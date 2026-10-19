#!/usr/bin/env python3
"""s56_s10_crossseed_spanfix.py — cross-seed spans under the manuscript convention.

Purpose:  Recomputes the SI Table S17 quantities from the archived per-frame
          series: the CV1 C3 span of the protein-DNA PMF under the table
          convention (occupied-support span: max-min of the C3 cumulant PMF
          over bins with occupancy >= 2; empty and single-count bins excluded)
          and the within-trajectory 2-ns window standard deviation. Runs for
          the published seed-42 trajectories (niraparib S2, neutral and
          protonated; from data/per_frame/) and the two additional seeds per
          state (7 and 123; s10_* files). The script first validates the
          pipeline against the frozen seed-42 values (spanfix data layer:
          c3_valid 26.926038 / 34.112682 kcal/mol; 2-ns window SD 7.4725
          kcal/mol; occupancy-bin counts 44/5/3 and 49/3/3) and then asserts
          cell-for-cell agreement with the frozen summary table.
Conventions: PyReweighting-1D amdweight_CE, T = 300 K, bin width 0.1 A,
          Emax = 20 kcal/mol, cutoff 2; weights = (beta*dV_D, 0, dV_D).
Inputs:   data/per_frame/{sys2_niraparib,sys2_nira_prot}_frames.csv.gz
          (5 columns) and data/per_frame/s10_{neu,prot}_{s7,s123}_frames.csv.gz
          (6 columns; time in ns). Column lookup is by header name.
Outputs:  data/analysis/s10/s10_crossseed_recomputed.csv (the same six
          quantities as the frozen s10_crossseed_summary.csv, recomputed from
          first principles with compacted column names)
Depends:  numpy; vendored tools/PyReweighting/PyReweighting-1D.py (scipy)
Run:      python3 scripts/s56_s10_crossseed_spanfix.py   (from the repository root)
"""
import contextlib
import csv
import gzip
import importlib.util
import io
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "per_frame"
OUT = ROOT / "data" / "analysis" / "s10"
OUT.mkdir(parents=True, exist_ok=True)

_spec = importlib.util.spec_from_file_location(
    "pyrew1d", ROOT / "tools" / "PyReweighting" / "PyReweighting-1D.py")
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)  # main() is guarded by __main__

BETA = 1.0 / (0.001987 * 300.0)
DISC, EMAX, CUT, T = 0.1, 20.0, 2, 300.0
WIN = 40  # 2-ns windows at 50-ps cadence

FROZEN = {  # (state, seed): span, window_sd, n_windows
    ("neutral", "42"): (26.9, 7.47, 13),
    ("neutral", "7"): (25.9, 12.38, 13),
    ("neutral", "123"): (46.1, 14.73, 11),
    ("protonated", "42"): (34.1, 4.62, 13),
    ("protonated", "7"): (49.3, 9.57, 13),
    ("protonated", "123"): (82.2, 18.32, 13),
}
FILES = {
    ("neutral", "42"): "sys2_niraparib_frames.csv.gz",
    ("protonated", "42"): "sys2_nira_prot_frames.csv.gz",
    ("neutral", "7"): "s10_neu_s7_frames.csv.gz",
    ("neutral", "123"): "s10_neu_s123_frames.csv.gz",
    ("protonated", "7"): "s10_prot_s7_frames.csv.gz",
    ("protonated", "123"): "s10_prot_s123_frames.csv.gz",
}


def load_csv(p):
    with gzip.open(p, "rt") as f:
        cols = f.readline().strip().split(",")
        arr = np.loadtxt(f)
    return arr[:, cols.index("cv1_prot_dna_com")], arr[:, cols.index("dV_D")]


def _edges(data):
    max_data = DISC * (int(np.amax(data) / DISC) + 1)
    min_data = DISC * (int(np.amin(data) / DISC) - 1)
    return mod.assignbins([min_data, max_data], DISC)


def span_full(cv, dv):
    data = np.asarray(cv, float)
    dV = np.asarray(dv, float)
    binsX = _edges(data)
    with contextlib.redirect_stdout(io.StringIO()):
        hist, _e, c1, c2, c3 = mod.reweight_CE(data, CUT, binsX, DISC, dV, T, False)
    c1 = -c1 / BETA
    c2 = -c2 / BETA
    c3 = -c3 / BETA
    with contextlib.redirect_stdout(io.StringIO()):
        pmf = mod.prephist(hist, T, EMAX)
    pmf_c3 = mod.normalize(pmf + c1 + c2 + c3, EMAX)
    valid = hist >= CUT
    return float(pmf_c3[valid].max() - pmf_c3[valid].min()), hist


def main():
    rows = []
    for (state, seed), fn in FILES.items():
        cv, dv = load_csv(DATA / fn)
        span, hist = span_full(cv, dv)
        spans = [span_full(cv[s:s + WIN], dv[s:s + WIN])[0]
                 for s in range(0, len(cv) - WIN + 1, WIN)]
        sd = float(np.std(spans))
        f_span, f_sd, f_nwin = FROZEN[(state, seed)]
        ok_span = abs(span - f_span) < 0.05
        ok_sd = abs(sd - f_sd) < 0.02
        ok_nwin = len(spans) == f_nwin
        print(f"{state:10s} seed {seed:>3s}: span={span:7.1f} (frozen {f_span:5.1f}, "
              f"{'OK' if ok_span else 'MISMATCH'})  winSD={sd:6.2f} (frozen {f_sd:5.2f}, "
              f"{'OK' if ok_sd else 'MISMATCH'})  windows={len(spans)} "
              f"(frozen {f_nwin}, {'OK' if ok_nwin else 'MISMATCH'})", flush=True)
        assert ok_span and ok_sd and ok_nwin, (state, seed)
        rows.append(dict(state=state, seed=seed, c3_valid=round(span, 1),
                         win_sd=round(sd, 2), n_windows=len(spans)))

    with open(OUT / "s10_crossseed_recomputed.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["state", "seed", "c3_valid", "win_sd", "n_windows"])
        w.writeheader()
        w.writerows(rows)
    # exact seed-42 convention self-check against the frozen spanfix values
    cv, dv = load_csv(DATA / FILES[("neutral", "42")])
    span, hist = span_full(cv, dv)
    n_empty = int((hist == 0).sum())
    n_count1 = int((hist == 1).sum())
    assert abs(span - 26.926038) < 2e-3, span
    assert n_empty == 5 and n_count1 == 3 and len(hist) == 44, (n_empty, n_count1, len(hist))
    print(f"seed-42 spanfix self-check: span={span:.6f} (frozen 26.926038), "
          f"bins={len(hist)} empty={n_empty} count1={n_count1} — ALL PASS")
    print("DONE ->", OUT / "s10_crossseed_recomputed.csv")


if __name__ == "__main__":
    main()
