#!/usr/bin/env python3
"""s67_c2_panel_sensitivity.py — C2 (second-order cumulant) panel correlations.

Reads the frozen harmonized summary (s62_harmonized_summary.csv), selects the
S2 CV1 C2 spans for the neutral and fully protonated five-system panels, and
recomputes the tie-aware Spearman correlation against the published trapping
rank vector with an exact two-tailed permutation test (n = 5, 120 relabelings,
one tied pair, midranks). Writes s67_c2_panel_association.csv.

This is a text-level sensitivity check on existing data; no simulation is run.
"""
import csv
import itertools
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/analysis/r6/s62_harmonized_summary.csv"
OUT = ROOT / "data/analysis/r6/s67_c2_panel_association.csv"

SYSTEMS = {
    "tala": "orig_s2_talazoparib_50ps",
    "nira": "sys2_niraparib",
    "ola": "sys2_olaparib",
    "ruca": "sys2_rucaparib",
    "veli": "orig_s2_veliparib",
    "nira_prot": "sys2_nira_prot",
    "ruca_prot": "sys2_ruca_prot",
    "veli_prot": "veli_prot_s2",
}
# weakest (veliparib) to strongest (talazoparib); olaparib/rucaparib tied
TRAPPING = [5, 4, 2.5, 2.5, 1]  # tala, nira, ola, ruca, veli


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
        mx = sum(seq) / n
        sxx = math.sqrt(sum((a - mx) ** 2 for a in seq))
        return sum((a - mx) * (b - my) for a, b in zip(seq, ry)) / (sxx * syy)

    obs = r_of(rx)
    cnt = 0
    for perm in itertools.permutations(ys):
        if max(perm) == min(perm):
            continue
        # permute the phenotype vector against the fixed span vector; ties in
        # the phenotype are carried through the midrank transform
        if abs(sum((a - sum(rx) / n) * (b - my) for a, b in zip(rx, midrank(perm))) /
               (math.sqrt(sum((a - sum(rx) / n) ** 2 for a in rx)) * syy)) >= abs(obs) - 1e-12:
            cnt += 1
    return round(obs, 3), round(cnt / math.factorial(n), 3), cnt


def main():
    rows = {r["system"]: r for r in csv.DictReader(open(SRC))}
    panels = {
        "neutral": ["tala", "nira", "ola", "ruca", "veli"],
        "fully_protonated": ["tala", "nira_prot", "ola", "ruca_prot", "veli_prot"],
    }
    out = []
    print(f"{'panel':18s} {'metric':5s} {'rho':>7s} {'p':>7s} {'cnt':>5s}")
    for panel, keys in panels.items():
        for metric in ("C2", "C3"):
            xs = [float(rows[SYSTEMS[k]][metric]) for k in keys]
            r, p, cnt = rho_p(xs, TRAPPING)
            out.append((panel, metric, r, p, cnt))
            print(f"{panel:18s} {metric:5s} {r:+7.3f} {p:7.3f} {cnt:5d}")
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["panel", "metric", "rho", "p_exact", "cnt120"])
        for row in out:
            w.writerow(row)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
