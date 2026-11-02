#!/usr/bin/env python3
"""s69_late_window.py — estimand/window lock for the major revision.

SUPERSEDED for the late-window numbers by s72_tala_late_window.py, which uses
the same seed-2026 input and is the provenance of Table S22 / the late-window
text.  This script is retained for provenance and for its first/late/last
estimator inventory; do not cite its span values as current.

Computes C1/C2/C3, q90 and t95 spans for three windows:
  first14  : frames   0-279 (14 ns)
  late14_26: frames 280-519 (12 ns)
  last10   : final 200 frames (10 ns)
for the eight S2 systems in the neutral/fully-protonated association panels,
and evaluates tie-aware exact Spearman associations for C2, q90 and C3.

P0 fix (2026-11-01): the archived 26-ns talazoparib trajectory
(``orig_s2_talazoparib_full``) was rebuilt with a different CV1 extraction and
cannot enter a panel with the other systems (s61_archive_reconciliation.csv).
Its production segment is also only 14 ns, so it has no late window. This
script therefore uses the independent seed-2026 26-ns run
(``tala_s2_26ns_frames.csv``, production mask step >= 6.0e6, 520 frames),
the same input as s72_tala_late_window.py. Every archive is production-masked
before window slicing.

Outputs:
  data/analysis/r6/s69_late_window_spans.csv
  data/analysis/r6/s69_late_window_assoc.csv
No simulation is run.
"""
import contextlib
import gzip
import importlib.util
import io
import itertools
import math
from pathlib import Path

import numpy as np

from production_mask import mask_production, assert_veli_prot

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "pyrew1d", ROOT / "tools/PyReweighting/PyReweighting-1D.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

BETA = 1.0 / (0.001987 * 300.0)
DISC, EMAX, CUT, T = 0.1, 20.0, 2, 300.0
OUT_SPANS = ROOT / "data/analysis/r6/s69_late_window_spans.csv"
OUT_ASSOC = ROOT / "data/analysis/r6/s69_late_window_assoc.csv"

# name -> (file stem, production-mask flag, gzip?)
FILES = {
    "talazoparib": ("tala_s2_26ns_frames.csv", True, False),
    "niraparib": ("sys2_niraparib_frames.csv.gz", False, True),
    "olaparib": ("sys2_olaparib_frames.csv.gz", False, True),
    "rucaparib": ("sys2_rucaparib_frames.csv.gz", False, True),
    "veliparib": ("orig_s2_veliparib_frames.csv.gz", False, True),
    "niraparib_prot": ("sys2_nira_prot_frames.csv.gz", False, True),
    "rucaparib_prot": ("sys2_ruca_prot_frames.csv.gz", False, True),
    "veliparib_prot": ("veli_prot_s2_frames.csv.gz", True, True),
}


def load(stem, prod_mask, gz):
    path = ROOT / "data/per_frame" / stem
    op = gzip.open if gz else open
    with op(path, "rt") as f:
        cols = f.readline().strip().split(",")
        rows = [ln.replace(",", " ").split() for ln in f if ln.strip()]
    arr = np.asarray(rows, float)
    name = stem.replace("_frames.csv.gz", "").replace("_frames.csv", "")
    arr, _ = mask_production(name, cols, arr)
    if name == "veli_prot_s2":
        assert_veli_prot(len(arr))
    if prod_mask:
        step_col = "step" if "step" in cols else "step_ps"
        arr = arr[arr[:, cols.index(step_col)] >= 6_000_000]
    cv_col = "cv1_prot_dna_com" if "cv1_prot_dna_com" in cols else "cv1"
    return arr[:, cols.index(cv_col)], arr[:, cols.index("dV_D")]


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
        if order == "C3":
            out["q90"] = float(np.quantile(vals, 0.95) - np.quantile(vals, 0.05))
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


def rho_p(xs, ys):
    n = len(xs)
    rx, ry = midrank(xs), midrank(ys)
    my = sum(ry) / n
    syy = math.sqrt(sum((b - my) ** 2 for b in ry))

    def r_of(seq):
        m = sum(seq) / n
        s = math.sqrt(sum((a - m) ** 2 for a in seq))
        return sum((a - m) * (b - my) for a, b in zip(seq, ry)) / (s * syy)

    obs = r_of(rx)
    cnt = 0
    for perm in itertools.permutations(xs):
        if max(perm) == min(perm):
            continue
        if abs(r_of(midrank(perm))) >= abs(obs) - 1e-12:
            cnt += 1
    return round(obs, 3), cnt


def main():
    series = {k: load(*v) for k, v in FILES.items()}
    # truncate every series to 26 ns (520 frames at 50 ps)
    for k in series:
        series[k] = (series[k][0][:520], series[k][1][:520])
        print(k, len(series[k][0]), flush=True)

    windows = {
        "first14": lambda n: (0, min(280, n)),
        "late14_26": lambda n: (280, min(520, n)),
        "last10": lambda n: (max(0, n - 200), n),
    }
    rows = []
    print(f"{'system':16s} {'window':10s} {'n':>4s} {'C2':>7s} {'C3':>7s} {'q90':>7s} {'t95':>7s}")
    for s, (cv, dv) in series.items():
        for wname, fn in windows.items():
            a, b = fn(len(cv))
            if b - a < 100:
                continue
            sp = spans(cv[a:b], dv[a:b])
            rows.append((s, wname, b - a, sp["C2"], sp["C3"], sp["q90"], sp["t95"]))
            print(f"{s:16s} {wname:10s} {b-a:4d} {sp['C2']:7.2f} {sp['C3']:7.2f} {sp['q90']:7.2f} {sp['t95']:7.2f}")

    with open(OUT_SPANS, "w", newline="") as f:
        import csv
        w = csv.writer(f)
        w.writerow(["system", "window", "n_frames", "C2", "C3", "q90", "t95"])
        w.writerows(rows)

    Y = [5.0, 4.0, 2.5, 2.5, 1.0]
    neutral = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
    prot = ["talazoparib", "niraparib_prot", "olaparib", "rucaparib_prot", "veliparib_prot"]
    assoc = []
    print("\nassociations")
    for wname in windows:
        for metric in ("C2", "q90", "C3"):
            vals = {(r[0], r[1]): r[3 if metric == "C2" else 4 if metric == "C3" else 5] for r in rows}
            for panel, systems in (("neutral", neutral), ("fully_protonated", prot)):
                xs = [vals[(s, wname)] for s in systems]
                rho, cnt = rho_p(xs, Y)
                assoc.append((wname, metric, panel, rho, cnt, cnt / 120.0))
                print(f"{wname:10s} {metric:3s} {panel:18s} rho={rho:+.3f} p={cnt}/120={cnt/120:.3f}")

    with open(OUT_ASSOC, "w", newline="") as f:
        import csv
        w = csv.writer(f)
        w.writerow(["window", "metric", "panel", "rho", "cnt120", "p_exact"])
        w.writerows(assoc)
    print(f"\nwrote {OUT_SPANS} and {OUT_ASSOC}")


if __name__ == "__main__":
    main()
