#!/usr/bin/env python3
"""s70_estimand_lock.py — C2/q90 estimand-lock table for the major revision.

Reads the frozen harmonized spans (s62_harmonized_summary.csv), computes
tie-aware exact Spearman associations for C2, C3, q90 and t95 in the neutral
and fully protonated panels, and writes s70_estimand_lock.csv.

This does NOT choose the estimand; it produces the decision table required
before the revision is written.
"""
import csv
import itertools
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/analysis/r6/s62_harmonized_summary.csv"
OUT = ROOT / "data/analysis/r6/s70_estimand_lock.csv"

NEUTRAL = ["orig_s2_talazoparib_50ps", "sys2_niraparib", "sys2_olaparib",
           "sys2_rucaparib", "orig_s2_veliparib"]
PROT = ["orig_s2_talazoparib_50ps", "sys2_nira_prot", "sys2_olaparib",
        "sys2_ruca_prot", "veli_prot_s2"]
Y = [5.0, 4.0, 2.5, 2.5, 1.0]


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
    rows = {r["system"]: r for r in csv.DictReader(open(SRC))}
    out = []
    print(f"{'metric':5s} {'panel':18s} {'rho':>7s} {'p':>7s} {'cnt':>4s}")
    for metric in ("C2", "C3", "q90", "t95"):
        for panel, systems in (("neutral", NEUTRAL), ("fully_protonated", PROT)):
            xs = [float(rows[s][metric]) for s in systems]
            r, c, p = exact(xs, Y)
            out.append((metric, panel, r, p, c))
            print(f"{metric:5s} {panel:18s} {r:+7.3f} {p:7.3f} {c:4d}")
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "panel", "rho", "p_exact", "cnt120"])
        w.writerows(out)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
