#!/usr/bin/env python3
"""s72_tala_late_window.py — replace the unavailable tala late-window with
the new 26-ns seed-2026 run.

Reads the new tala per-frame CSV (production mask step>=6.0e6) plus the four
archived 26-ns neutral systems, truncates all to 520 frames (26 ns), and
computes C1/C2/C3 and 5-95% quantile spans for:
  first14  : frames   0-279
  late14_26: frames 280-519
  last10   : frames 320-519
Tie-aware exact Spearman associations are evaluated for C2, C2/q90, C3 and
C3/q90 in the neutral and fully protonated panels.

Outputs:
  data/analysis/r6/s72_tala_late_window_spans.csv
  data/analysis/r6/s72_tala_late_window_assoc.csv
"""
import contextlib
import csv
import gzip
import importlib.util
import io
import itertools
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "pyrew1d", ROOT / "tools/PyReweighting/PyReweighting-1D.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

BETA = 1.0 / (0.001987 * 300.0)
DISC, EMAX, CUT, T = 0.1, 20.0, 2, 300.0
OUT_SPANS = ROOT / "data/analysis/r6/s72_tala_late_window_spans.csv"
OUT_ASSOC = ROOT / "data/analysis/r6/s72_tala_late_window_assoc.csv"

# name -> (file stem, truncated-to-production flag, gz?)
FILES = {
    "talazoparib": ("tala_s2_26ns_frames.csv", True, False),
    "niraparib": ("sys2_niraparib_frames.csv.gz", False, True),
    "olaparib": ("sys2_olaparib_frames.csv.gz", False, True),
    "rucaparib": ("sys2_rucaparib_frames.csv.gz", False, True),
    "veliparib": ("orig_s2_veliparib_frames.csv.gz", False, True),
    "niraparib_prot": ("sys2_nira_prot_frames.csv.gz", False, True),
    "rucaparib_prot": ("sys2_ruca_prot_frames.csv.gz", False, True),
    "veliparib_prot": ("veli_prot_s2_frames.csv.gz", False, True),
}


def load(stem, prod_mask, gz):
    path = ROOT / "data/per_frame" / stem
    op = gzip.open if gz else open
    with op(path, "rt") as f:
        header = f.readline().strip().split(",")
        rows = [ln.replace(",", " ").split() for ln in f if ln.strip()]
    arr = np.asarray(rows, float)
    if prod_mask:
        arr = arr[arr[:, header.index("step")] >= 6_000_000]
    cv_col = "cv1_prot_dna_com" if "cv1_prot_dna_com" in header else "cv1"
    return arr[:520, header.index(cv_col)], arr[:520, header.index("dV_D")]


def spans(cv, dv):
    data, dV = np.asarray(cv, float), np.asarray(dv, float)
    mx = DISC * (int(np.amax(data) / DISC) + 1)
    mn = DISC * (int(np.amin(data) / DISC) - 1)
    binsX = mod.assignbins([mn, mx], DISC)
    with contextlib.redirect_stdout(io.StringIO()):
        hist, _e, c1, c2, c3 = mod.reweight_CE(data, CUT, binsX, DISC, dV, T, False)
    c1, c2, c3 = -c1 / BETA, -c2 / BETA, -c3 / BETA
    pmf = mod.prephist(hist, T, EMAX)
    out = {}
    for order, coef in (("C1", c1), ("C2", c1 + c2), ("C3", c1 + c2 + c3)):
        p = mod.normalize(pmf + coef, EMAX)
        vals = p[hist >= CUT]
        out[order] = float(vals.max() - vals.min())
        out[order + "_q90"] = float(np.quantile(vals, 0.95) - np.quantile(vals, 0.05))
        if order == "C3":
            out["t95"] = float(np.percentile(vals, 97.5) - np.percentile(vals, 2.5))
    return out


def midrank(v):
    v = list(v)
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for t in range(i, j + 1):
            r[order[t]] = (i + j + 2) / 2.0
        i = j + 1
    return r


def rho(xs, ys):
    n = len(xs)
    rx, ry = midrank(xs), midrank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    return sum((a - mx) * (b - my) for a, b in zip(rx, ry)) / math.sqrt(
        sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))


def exact(xs, ys):
    obs = rho(xs, ys)
    cnt = sum(1 for p in itertools.permutations(ys)
              if abs(rho(xs, list(p))) >= abs(obs) - 1e-12)
    return round(obs, 3), cnt, round(cnt / 120.0, 3)


def main():
    series = {}
    for name, (stem, pm, gz) in FILES.items():
        cv, dv = load(stem, pm, gz)
        series[name] = (cv, dv)
        print(f"{name:16s} n={len(cv)}", flush=True)

    windows = {
        "first14": lambda n: (0, 280),
        "late14_26": lambda n: (280, 520),
        "last10": lambda n: (320, 520),
    }
    rows = []
    print(f"\n{'system':16s} {'window':9s} {'n':>3s} {'C2':>7s} {'C2q90':>7s} {'C3':>7s} {'C3q90':>7s}")
    for s, (cv, dv) in series.items():
        for wn, fn in windows.items():
            a, b = fn(len(cv))
            sp = spans(cv[a:b], dv[a:b])
            rows.append((s, wn, b - a, sp["C2"], sp["C2_q90"], sp["C3"],
                         sp["C3_q90"], sp["C1"], sp["t95"]))
            print(f"{s:16s} {wn:9s} {b-a:3d} {sp['C2']:7.2f} {sp['C2_q90']:7.2f} "
                  f"{sp['C3']:7.2f} {sp['C3_q90']:7.2f}", flush=True)

    with open(OUT_SPANS, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["system", "window", "n_frames", "C2", "C2_q90", "C3",
                    "C3_q90", "C1", "t95"])
        w.writerows(rows)

    Y = [5.0, 4.0, 2.5, 2.5, 1.0]
    neutral = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
    prot = ["talazoparib", "niraparib_prot", "olaparib", "rucaparib_prot", "veliparib_prot"]
    assoc = []
    print("\nassociations")
    for wn in windows:
        for metric in ("C2", "C2_q90", "C3", "C3_q90"):
            vals = {(r[0], r[1]): r[3 + ("C2", "C2_q90", "C3", "C3_q90").index(metric)]
                    for r in rows}
            for panel, systems in (("neutral", neutral), ("fully_protonated", prot)):
                xs = [vals[(s, wn)] for s in systems]
                r, c, p = exact(xs, Y)
                assoc.append((wn, metric, panel, r, c, p))
                print(f"{wn:9s} {metric:7s} {panel:18s} rho={r:+.3f} p={p:.3f} ({c}/120)")
    with open(OUT_ASSOC, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["window", "metric", "panel", "rho", "cnt120", "p_exact"])
        w.writerows(assoc)
    print(f"\nwrote {OUT_SPANS} and {OUT_ASSOC}")


if __name__ == "__main__":
    main()
