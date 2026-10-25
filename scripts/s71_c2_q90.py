#!/usr/bin/env python3
"""s71_c2_q90.py — bounded q90 span under the C2 PMF.

Computes, for the eight S2 systems and the first-14-ns harmonized window,
the 5–95% quantile span of the C2-reweighted PMF, alongside C2 min–max and
C3 q90, then evaluates tie-aware exact associations for the neutral and
fully protonated panels.

Outputs:
  data/analysis/r6/s71_c2_q90_spans.csv
  data/analysis/r6/s71_c2_q90_assoc.csv
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
OUT_SPANS = ROOT / "data/analysis/r6/s71_c2_q90_spans.csv"
OUT_ASSOC = ROOT / "data/analysis/r6/s71_c2_q90_assoc.csv"

FILES = {
    "talazoparib": ("orig_s2_talazoparib", 50),
    "niraparib": ("sys2_niraparib", 1),
    "olaparib": ("sys2_olaparib", 1),
    "rucaparib": ("sys2_rucaparib", 1),
    "veliparib": ("orig_s2_veliparib", 1),
    "niraparib_prot": ("sys2_nira_prot", 1),
    "rucaparib_prot": ("sys2_ruca_prot", 1),
    "veliparib_prot": ("veli_prot_s2", 1),
}


def load(name, stride):
    with gzip.open(ROOT / "data/per_frame" / f"{name}_frames.csv.gz", "rt") as f:
        cols = f.readline().strip().split(",")
        rows = [ln.replace(",", " ").split() for ln in f if ln.strip()]
    arr = np.asarray(rows, float)
    return arr[::stride, cols.index("cv1_prot_dna_com")], arr[::stride, cols.index("dV_D")]


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
    for order, coef in (("C2", c1 + c2), ("C3", c1 + c2 + c3)):
        p = mod.normalize(pmf + coef, EMAX)
        vals = p[hist >= CUT]
        out[order + "_mm"] = float(vals.max() - vals.min())
        out[order + "_q90"] = float(np.quantile(vals, 0.95) - np.quantile(vals, 0.05))
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
    rows = []
    for s, (name, stride) in FILES.items():
        cv, dv = load(name, stride)
        cv, dv = cv[:280], dv[:280]  # first 14 ns at 50 ps
        sp = spans(cv, dv)
        rows.append((s, *[sp[k] for k in ("C2_mm", "C2_q90", "C3_mm", "C3_q90")]))
        print(f"{s:16s} C2mm={sp['C2_mm']:6.2f} C2q90={sp['C2_q90']:6.2f} "
              f"C3mm={sp['C3_mm']:6.2f} C3q90={sp['C3_q90']:6.2f}")

    with open(OUT_SPANS, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["system", "C2_mm", "C2_q90", "C3_mm", "C3_q90"])
        w.writerows(rows)

    V = {r[0]: dict(zip(("C2_mm", "C2_q90", "C3_mm", "C3_q90"), r[1:])) for r in rows}
    Y = [5.0, 4.0, 2.5, 2.5, 1.0]
    neutral = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
    prot = ["talazoparib", "niraparib_prot", "olaparib", "rucaparib_prot", "veliparib_prot"]
    out = []
    print("\nassociations")
    for metric in ("C2_mm", "C2_q90", "C3_mm", "C3_q90"):
        for panel, systems in (("neutral", neutral), ("fully_protonated", prot)):
            xs = [V[s][metric] for s in systems]
            r, c, p = exact(xs, Y)
            out.append((metric, panel, r, p, c))
            print(f"{metric:7s} {panel:18s} rho={r:+.3f} p={p:.3f} ({c}/120)")
    with open(OUT_ASSOC, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "panel", "rho", "p_exact", "cnt120"])
        w.writerows(out)
    print(f"\nwrote {OUT_SPANS} and {OUT_ASSOC}")


if __name__ == "__main__":
    main()
