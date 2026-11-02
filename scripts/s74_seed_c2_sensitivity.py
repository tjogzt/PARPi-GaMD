#!/usr/bin/env python3
"""s74_seed_c2_sensitivity.py — C2/q90 span sensitivity to archived seeds.

The earlier seed summaries (s60, s66) reported C3/min--max spans.  Reviewer
item P1 asks for the revision-primary estimand (C2 with a 5--95% quantile
span) on every archived seed.  This script recomputes, under the harmonized
50-ps cadence and first-14-ns production window:

  * talazoparib S2: original archive (seed 42), seed 7, seed 123, and the
    independent 26-ns seed 2026;
  * niraparib S2 neutral and protonated: original archive (seed 42),
    seed 7, seed 123;
  * veliparib S2 neutral: original archive (seed 42), seed 7, seed 123.

It then evaluates the five-system exact association for every combination of
the available seeds (one value substituted per system) in the neutral and
state-corrected panels, and reports the range of rho and exact p.

Outputs
  data/analysis/r6/s74_c2_q90_seeds.csv
  data/analysis/r6/s74_c2_q90_seed_assoc.csv
"""
import contextlib
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
FRAMES_14NS = 280


def spans(cv, dv):
    data, dV = np.asarray(cv, float), np.asarray(dv, float)
    mx = DISC * (int(np.amax(data) / DISC) + 1)
    mn = DISC * (int(np.amin(data) / DISC) - 1)
    binsX = mod.assignbins([mn, mx], DISC)
    with contextlib.redirect_stdout(io.StringIO()):
        hist, _e, c1, c2, _c3 = mod.reweight_CE(data, CUT, binsX, DISC, dV, T, False)
        c1, c2 = -c1 / BETA, -c2 / BETA
        pmf = mod.prephist(hist, T, EMAX)
        p2 = mod.normalize(pmf + c1 + c2, EMAX)
    valid = hist >= CUT
    vals = p2[valid]
    # Population-weighted 5--95% interval: bin weights exp(-beta G), not
    # equal bin weights.  This is the alternative the reviewer requested in
    # place of a fixed-support DeltaG, which is undefined when supports differ.
    w = np.exp(-BETA * vals)
    w = w / w.sum()
    order = np.argsort(vals)
    cw = np.cumsum(w[order])
    q05 = vals[order][int(np.searchsorted(cw, 0.05))]
    q95 = vals[order][int(np.searchsorted(cw, 0.95))]
    return (float(vals.max() - vals.min()),
            float(np.quantile(vals, 0.95) - np.quantile(vals, 0.05)),
            float(q95 - q05),
            int(valid.sum()), int(hist.size))


def load_gz(name, stride=1, step_start=None):
    with gzip.open(ROOT / "data/per_frame" / f"{name}_frames.csv.gz", "rt") as f:
        cols = f.readline().strip().split(",")
        arr = np.array([[float(x) for x in ln.replace(",", " ").split()]
                        for ln in f if ln.strip()])
    if step_start is not None:
        step_col = "step" if "step" in cols else ("step_ps" if "step_ps" in cols else None)
        if step_col is None:
            raise KeyError(f"{name}: no MD-step column")
        arr = arr[arr[:, cols.index(step_col)] >= step_start]
    if stride > 1:
        arr = arr[::stride]
    cv = arr[:, cols.index("cv1_prot_dna_com")]
    dv = arr[:, cols.index("dV_D")]
    return cv[:FRAMES_14NS], dv[:FRAMES_14NS]


def load_s10(name):
    with gzip.open(ROOT / "data/per_frame" / f"{name}_frames.csv.gz", "rt") as f:
        cols = f.readline().strip().split(",")
        arr = np.array([[float(x) for x in ln.replace(",", " ").split()]
                        for ln in f if ln.strip()])
    keep = arr[:, cols.index("time_ns")] >= 12.0  # production begins at 12 ns
    arr = arr[keep]
    return (arr[:FRAMES_14NS, cols.index("cv1_prot_dna_com")],
            arr[:FRAMES_14NS, cols.index("dV_D")])


def load_multiseed(prefix, seed):
    d = np.loadtxt(ROOT / f"data/per_frame/multiseed/{prefix}_s2_{seed}_cv.csv",
                   delimiter=",", skiprows=1)
    # 760 rows, 50 ps, production frames 241--760 (12--38 ns)
    return d[240:520, 3], d[240:520, 2]


def load_tala26():
    path = ROOT / "data/per_frame/tala_s2_26ns_frames.csv"
    with path.open() as f:
        cols = f.readline().strip().split(",")
        arr = np.array([[float(x) for x in ln.split(",")] for ln in f if ln.strip()])
    arr = arr[arr[:, cols.index("step")] >= 6_000_000]
    return arr[:FRAMES_14NS, cols.index("cv1")], arr[:FRAMES_14NS, cols.index("dV_D")]


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
    return round(obs, 3), round(cnt / 120.0, 3)


def main():
    specs = {
        "talazoparib": [
            ("seed42", *load_gz("orig_s2_talazoparib", stride=50)),
            ("seed7", *load_multiseed("tala", "seed7")),
            ("seed123", *load_multiseed("tala", "seed123")),
            ("seed2026", *load_tala26()),
        ],
        "niraparib": [
            ("seed42", *load_gz("sys2_niraparib")),
            ("seed7", *load_s10("s10_neu_s7")),
            ("seed123", *load_s10("s10_neu_s123")),
        ],
        "niraparib_prot": [
            ("seed42", *load_gz("sys2_nira_prot")),
            ("seed7", *load_s10("s10_prot_s7")),
            ("seed123", *load_s10("s10_prot_s123")),
        ],
        "veliparib": [
            ("seed42", *load_gz("orig_s2_veliparib")),
            ("seed7", *load_multiseed("veli", "seed7")),
            ("seed123", *load_multiseed("veli", "seed123")),
        ],
        "veliparib_prot": [
            ("seed42", *load_gz("veli_prot_s2", step_start=6_000_000)),
        ],
    }

    rows = []
    values = {}
    pop_values = {}
    print(f"{'system':16s} {'seed':9s} {'C2mm':>7s} {'C2q90':>7s} {'q90pop':>7s} {'nvalid':>7s}/{'nbins':<6s}")
    for system, entries in specs.items():
        values[system] = {}
        pop_values[system] = {}
        for seed, cv, dv in entries:
            mm, q90, q90pop, nvalid, nbins = spans(cv, dv)
            values[system][seed] = q90
            pop_values[system][seed] = q90pop
            rows.append((system, seed, mm, q90, q90pop, nvalid, nbins, len(cv)))
            print(f"{system:16s} {seed:9s} {mm:7.2f} {q90:7.2f} {q90pop:7.2f} "
                  f"{nvalid:7d}/{nbins:<6d}")

    out = ROOT / "data/analysis/r6/s74_c2_q90_seeds.csv"
    with out.open("w") as fh:
        fh.write("system,seed,C2_minmax,C2_q90,C2_q90_pop,n_valid_bins,n_bins,n_frames\n")
        for r in rows:
            fh.write(f"{r[0]},{r[1]},{r[2]:.4f},{r[3]:.4f},{r[4]:.4f},{r[5]},{r[6]},{r[7]}\n")

    Y = [5.0, 4.0, 2.5, 2.5, 1.0]
    panels = {
        "neutral": ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"],
        "state_corrected": ["talazoparib", "niraparib_prot", "olaparib",
                            "rucaparib_prot", "veliparib_prot"],
    }
    fixed = {"olaparib": 11.23087818532426, "rucaparib": 8.781788783742979,
             "rucaparib_prot": 15.121140397267}
    combo_rows = []
    print("\nseed-combination association range")
    for panel, systems in panels.items():
        choices = []
        for s in systems:
            if s in values:
                choices.append(sorted(values[s].items()))
            else:
                choices.append([("seed42", fixed[s])])
        for metric, table in (("C2_q90", values), ("C2_q90_pop", pop_values)):
            rhos, ps = [], []
            for combo in itertools.product(*choices):
                xs = []
                for s, (seed, _v) in zip(systems, combo):
                    xs.append(table[s][seed] if s in table else fixed[s])
                r, p = exact(xs, Y)
                rhos.append(r)
                ps.append(p)
                combo_rows.append((panel, metric,
                                   ";".join(f"{s}:{seed}" for s, (seed, _v) in zip(systems, combo)),
                                   r, p))
            print(f"{panel:16s} {metric:11s} n_combos={len(rhos):3d} "
                  f"rho [{min(rhos):+.3f}, {max(rhos):+.3f}] "
                  f"p [{min(ps):.3f}, {max(ps):.3f}]")

    out2 = ROOT / "data/analysis/r6/s74_c2_q90_seed_assoc.csv"
    with out2.open("w") as fh:
        fh.write("panel,metric,seed_combination,rho,p_exact\n")
        for r in combo_rows:
            fh.write(f"{r[0]},{r[1]},\"{r[2]}\",{r[3]:+.3f},{r[4]:.3f}\n")
    print(f"\nwrote {out} and {out2}")


if __name__ == "__main__":
    main()
