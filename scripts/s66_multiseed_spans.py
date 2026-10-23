#!/usr/bin/env python3
"""s66_multiseed_spans.py — spans for the new tala/veli seed-7 and seed-123 runs.

Reads per-frame CV CSVs extracted on the GPU host (50-ps frames, 12-ns
equilibration + 26-ns production) and computes Table-S5 cumulant spans
(min--max C3, q90, t95) on:
  * the harmonized 14-ns production window (frames 241-520)
  * the full 26-ns production window (frames 241-760)
Outputs data/analysis/r6/s66_multiseed_spans.csv
"""
import contextlib
import importlib.util
import io
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "pyrew1d", ROOT / "tools/PyReweighting/PyReweighting-1D.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

BETA = 1.0 / (0.001987 * 300.0)
DISC, EMAX, CUT, T = 0.1, 20.0, 2, 300.0


def spans(cv, dv):
    data, dV = np.asarray(cv, float), np.asarray(dv, float)
    mx = DISC * (int(np.amax(data) / DISC) + 1)
    mn = DISC * (int(np.amin(data) / DISC) - 1)
    binsX = mod.assignbins([mn, mx], DISC)
    with contextlib.redirect_stdout(io.StringIO()):
        hist, _e, c1, c2, c3 = mod.reweight_CE(data, CUT, binsX, DISC, dV, T, False)
    c1, c2, c3 = -c1 / BETA, -c2 / BETA, -c3 / BETA
    pmf = mod.prephist(hist, T, EMAX)
    p3 = mod.normalize(pmf + c1 + c2 + c3, EMAX)
    v = p3[hist >= CUT]
    return dict(mm=float(v.max() - v.min()),
                q90=float(np.quantile(v, 0.95) - np.quantile(v, 0.05)),
                t95=float(np.percentile(v, 97.5) - np.percentile(v, 2.5)))


rows = []
for system, seed in [("talazoparib", "seed7"), ("talazoparib", "seed123"),
                     ("veliparib", "seed7"), ("veliparib", "seed123")]:
    f = ROOT / f"data/per_frame/multiseed/{'tala' if system=='talazoparib' else 'veli'}_s2_{seed}_cv.csv"
    d = np.loadtxt(f, delimiter=",", skiprows=1)
    cv, dv = d[:, 3], d[:, 2]
    for label, sl in [("prod14", slice(240, 520)), ("prod26", slice(240, 760))]:
        s = spans(cv[sl], dv[sl])
        for est, val in s.items():
            rows.append((system, seed, label, est, val))
        print(system, seed, label, {k: round(v, 2) for k, v in s.items()})

out = ROOT / "data/analysis/r6/s66_multiseed_spans.csv"
with open(out, "w") as fh:
    fh.write("system,seed,window,estimator,span_kcal\n")
    for r in rows:
        fh.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]:.4f}\n")
print("wrote", out)
