#!/usr/bin/env python3
"""s57_manuscript_number_audit.py — manuscript <-> data-layer number audit (v2).

Validates every headline number of the JCIM manuscript + SI against the frozen
data layer, and re-derives the statistical claims (associations, LOO,
counterfactual placement, perturbation grid, permutation minima) from the
published table values with an independent tie-aware exact-permutation engine.

Table S5 truth assembly:
  - s2_dbe_final.csv        (source-tag aware: exact_log_dbe / nira_prot_dbe /
                             veli_prot_dbe / rerun_exact_dbe)
  - r6c_s2_spanfix.csv      (CV1 valid spans: niraparib neutral, rucaparib prot)
  - r6c_s2_cv2_spanfix.csv  (CV2 valid spans: same two systems)

Tiers: A = value check from the data layer; B = statistical re-derivation;
C = presence of the formatted number in the LaTeX source (drift tripwire).

Run (from the repository root; tex outside the repo):
  python3 scripts/s57_manuscript_number_audit.py --tex-dir ../../manuscripts/jcim
"""
import argparse
import csv
import itertools
import math
import re
import sys
from pathlib import Path
from statistics import stdev

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

ap = argparse.ArgumentParser()
ap.add_argument("--tex-dir", default=str(ROOT / "../../manuscripts/jcim"))
args = ap.parse_args()
TEX = Path(args.tex_dir).resolve()

results = []


def rec(tier, cid, ok, detail=""):
    results.append((tier, cid, "PASS" if ok else "FAIL", str(detail)))


# ---------------------------------------------------------------- statistics
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
    for perm in itertools.permutations(xs):
        if max(perm) == min(perm):
            continue
        if abs(r_of(midrank(perm))) >= abs(obs) - 1e-12:
            cnt += 1
    return round(obs, 3), round(cnt / math.factorial(n), 3), cnt


# ---------------------------------------------------------------- load layer
final = {}
for r in csv.DictReader(open(DATA / "analysis/s2_dbe_final.csv")):
    src = r["source"]
    lig = r["ligand"]
    key = {"nira_prot_dbe": "niraparib_prot", "veli_prot_dbe": "veliparib_prot"}.get(src, lig)
    final[(key, r["metric"])] = (float(r["C1"]), float(r["C2"]), float(r["C3"]), src)

estim = {}
for r in csv.DictReader(open(DATA / "analysis/r6/r6c_estimatrix.csv")):
    try:
        estim[(r["panel"], r["metric"], r["estimator"])] = (float(r["rho"]), float(r["p"]))
    except ValueError:
        estim[(r["panel"], r["metric"], r["estimator"])] = None

lat = list(csv.DictReader(open(DATA / "analysis/r6/r9_state_lattice.csv")))
win = {(r["tag"], int(r["window_ns"])): float(r["c3_span_sd"])
       for r in csv.DictReader(open(DATA / "analysis/r6/r6c_window_timescale.csv"))}
s10 = list(csv.DictReader(open(DATA / "analysis/s10/s10_crossseed_summary.csv")))
pc = {(float(r["delta_kcal"]), float(r["sigma_kcal"])): float(r["detect_rate"])
      for r in csv.DictReader(open(DATA / "analysis/r6/c7g_positive_control.csv"))}
crossdock = {r["ligand"].lower(): float(r["rmsd_A"])
             for r in csv.DictReader(open(DATA / "analysis/crossdock_rmsd_results.csv"))}
ordi = list(csv.DictReader(open(DATA / "analysis/r6/r6d_ordinal_sensitivity.csv")))
s1seeds = {r["ligand"]: (r["C1"], r["C2"], r["C3"])
           for r in csv.DictReader(open(DATA / "analysis/r6/r6c_s1_seeds_canonical.csv"))}
meta = list(csv.DictReader(open(DATA / "trajectory_metadata.csv")))

ms = (TEX / "jcim_manuscript.tex").read_text()
si = (TEX / "si.tex").read_text()

# ---------------------------------------------------------------- S5 parse
LIG = {"APO": "APO", "Talazoparib": "talazoparib", "Niraparib": "niraparib",
       "Niraparib (prot.)": "niraparib_prot", "Olaparib": "olaparib",
       "Rucaparib": "rucaparib", "Rucaparib (prot.)": "rucaparib_prot",
       "Veliparib": "veliparib", "Veliparib (prot.)": "veliparib_prot",
       "AZD5305": "AZD5305", "Fluzoparib": "fluzoparib",
       "Pamiparib": "pamiparib", "Senaparib": "senaparib"}

s5blk = si[si.find("Table S5: S2 spans"):si.find("\\end{table}", si.find("Table S5: S2 spans"))]
s5val = {}
pat = re.compile(r"^([A-Za-z0-9 .()]+?)\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&\s*"
                 r"([\d.]+)\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&")
for line in s5blk.splitlines():
    m = pat.match(line)
    if not m:
        continue
    lig = LIG.get(m.group(1).strip())
    if lig:
        s5val[lig] = tuple(float(m.group(i)) for i in range(2, 8))
rec("A", "S5 rows parsed", len(s5val) == 13, f"{sorted(s5val)}")

# truth assembly
truth = {}
for (lig, metric), v in final.items():
    if metric in ("CV1", "CV2"):
        truth[(lig, metric)] = v[:3]
cv1fix = {r["system"]: (float(r["c1_valid"]), float(r["c2_valid"]), float(r["c3_valid"]))
          for r in csv.DictReader(open(DATA / "analysis/r6/r6c_s2_spanfix.csv"))}
cv2fix = {r["system"]: (float(r["c1"]), float(r["c2"]), float(r["c3"]))
          for r in csv.DictReader(open(DATA / "analysis/r6/r6c_s2_cv2_spanfix.csv"))}
for sysid, lig in [("niraparib", "niraparib"), ("rucaparib_prot", "rucaparib_prot")]:
    truth.setdefault((lig, "CV1"), cv1fix[sysid])
    truth.setdefault((lig, "CV2"), cv2fix[sysid])

ncell = nok = 0
bad = []
for lig, vals in s5val.items():
    for j in range(6):
        metric = "CV1" if j < 3 else "CV2"
        t = truth.get((lig, metric))
        if t is None:
            bad.append(f"{lig}/{metric}: no truth")
            continue
        want = t[j % 3]
        got = vals[j]
        ncell += 1
        if abs(got - round(want, 1)) < 0.051:
            nok += 1
        else:
            bad.append(f"{lig}/{metric} col{j%3}: table {got} vs data {want:.2f}")
rec("A", "S5 cells vs data", nok == ncell and ncell == 78,
    f"{nok}/{ncell}; issues: {bad[:6]}")

# ------------------------------------------------- S5 caption numbers
c3_1 = [v[2] for v in s5val.values()]
rec("A", "S5 CV1 C3 range", abs(min(c3_1) - 19.5) < 0.06 and abs(max(c3_1) - 49.5) < 0.06,
    f"{min(c3_1):.1f}--{max(c3_1):.1f} (claimed 19.5--49.5)")
neu = ["APO", "talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib", "AZD5305"]
c2c3 = [s5val[k][5] for k in neu]
rec("A", "S5 CV2 C3 neutral range", abs(min(c2c3) - 13.5) < 0.06 and abs(max(c2c3) - 41.2) < 0.06,
    f"{min(c2c3):.1f}--{max(c2c3):.1f} (claimed 13.5--41.2)")
rec("A", "S5 ruca prot CV2 C3", abs(s5val["rucaparib_prot"][5] - 45.1) < 0.06,
    f"{s5val['rucaparib_prot'][5]} (claimed 45.1)")
c1_1 = [v[0] for v in s5val.values()]
c1_2 = [v[3] for v in s5val.values()]
rec("A", "S5 C1 clusters", abs(min(c1_1) - 3.9) < 0.06 and abs(max(c1_1) - 6.8) < 0.06
    and abs(min(c1_2) - 2.4) < 0.06 and abs(max(c1_2) - 9.2) < 0.06,
    f"CV1 {min(c1_1):.1f}-{max(c1_1):.1f}; CV2 {min(c1_2):.1f}-{max(c1_2):.1f} (claimed 3.9-6.8; 2.4-9.2)")

# --------------------------------------------- separations & panel tail
sub = sorted([s5val["talazoparib"][2], s5val["veliparib_prot"][2],
              s5val["olaparib"][2], s5val["niraparib_prot"][2]])
diffs = [round(b - a, 1) for a, b in zip(sub, sub[1:])]
rec("A", "separation band", abs(min(diffs) - 2.6) < 0.06 and abs(max(diffs) - 5.3) < 0.06,
    f"{sub} -> {diffs} (claimed 2.6--5.3)")
ap5 = sorted([s5val[k][2] for k in ("talazoparib", "niraparib_prot", "olaparib",
                                      "rucaparib_prot", "veliparib_prot")])
rec("A", "association-panel span", abs(ap5[0] - 22.9) < 0.06 and abs(ap5[-1] - 48.7) < 0.06
    and abs(ap5[-1] - ap5[-2] - 14.6) < 0.06,
    f"{ap5} (claimed 22.9--48.7, tail 14.6)")

# --------------------------------------------- C2-C3 rank correlation
def spearman_point(xs, ys):
    rx, ry = midrank(xs), midrank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    syy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return round(sxy / (sxx * syy), 3)

r23 = spearman_point([s5val[k][2] for k in s5val], [s5val[k][1] for k in s5val])
rec("B", "C2-C3 rank correlation", abs(r23 - 0.60) <= 0.03, f"{r23} (claimed +0.60)")

# --------------------------------------------- window band (five systems)
wsd = [v for (tag, w), v in win.items() if w == 2]
rec("A", "2-ns window band", abs(min(wsd) - 4.6) < 0.06 and abs(max(wsd) - 18.4) < 0.06,
    f"{min(wsd):.2f}--{max(wsd):.2f} (claimed 4.6--18.4)")

# --------------------------------------------- S17 cells & means
claims = {("neutral", "42"): (26.9, 7.47, 13), ("neutral", "7"): (25.9, 12.38, 13),
          ("neutral", "123"): (46.1, 14.73, 11), ("protonated", "42"): (34.1, 4.62, 13),
          ("protonated", "7"): (49.3, 9.57, 13), ("protonated", "123"): (82.2, 18.32, 13)}
detail = []
ok17 = True
for r in s10:
    key = (r["system"], r["seed"])
    if key in claims:
        want = claims[key]
        got = (float(r["c3_valid_kcal"]), float(r["win_sd_2ns"]), int(r["n_windows"]))
        good = abs(got[0] - want[0]) < 0.06 and abs(got[1] - want[1]) < 0.02 and got[2] == want[2]
        ok17 &= good
        if not good:
            detail.append(f"{key}: {got} vs {want}")
rec("A", "S17 six cells", ok17, "; ".join(detail) or "6/6")

means = {state: [float(r["c3_valid_kcal"]) for r in s10 if r["system"] == state]
         for state in ("neutral", "protonated")}
mu = {k: sum(v) / len(v) for k, v in means.items()}
sd = {k: stdev(v) for k, v in means.items()}
rng = {k: max(v) - min(v) for k, v in means.items()}
rec("B", "S10 means/SD", abs(mu["neutral"] - 33.0) < 0.06 and abs(sd["neutral"] - 11.4) < 0.06
    and abs(mu["protonated"] - 55.2) < 0.06 and abs(sd["protonated"] - 24.6) < 0.06,
    f"{mu['neutral']:.1f}±{sd['neutral']:.1f} / {mu['protonated']:.1f}±{sd['protonated']:.1f}")
rec("B", "S10 mean difference", abs((mu["protonated"] - mu["neutral"]) - 22.2) < 0.06,
    f"{mu['protonated']-mu['neutral']:.1f} (claimed 22.2)")
rec("B", "S10 ranges", abs(rng["neutral"] - 20.2) < 0.06 and abs(rng["protonated"] - 48.1) < 0.06,
    f"{rng['neutral']:.1f}/{rng['protonated']:.1f} (claimed 20.2/48.1)")
pct = {k: 100 * sd[k] / mu[k] for k in mu}
rec("B", "S10 percent semantics",
    abs(pct["neutral"] - 35) < 1.5 and abs(pct["protonated"] - 45) < 1.5,
    f"SD/mean={pct['neutral']:.1f}%/{pct['protonated']:.1f}% -> '35-45%' is the SD/mean")
rec("C", "S10 pct attribution fixed",
    "standard deviations of approximately 35--45\\% of the respective means" in si
    and "spreads of 20.2 and 48.1~kcal/mol, approximately 35--45" not in si,
    "SI now attributes 35-45% to the SD; old conflated phrasing absent")

# --------------------------------------------- associations from S5
Y = [5, 4, 2.5, 2.5, 1]
P = {"neutral": ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"],
     "audit": ["talazoparib", "niraparib_prot", "olaparib", "rucaparib_prot", "veliparib_prot"]}
for pname, keys in P.items():
    r, p, _ = rho_p([s5val[k][2] for k in keys], Y)
    exp = {"neutral": (-0.667, 0.267), "audit": (-0.205, 0.767)}[pname]
    rec("B", f"rho S2CV1 {pname}", abs(r - exp[0]) < 0.006 and abs(p - exp[1]) < 0.006,
        f"{r}/{p} (claimed {exp[0]}/{exp[1]})")

loo = {}
for i, k in enumerate(P["audit"]):
    xs = [s5val[kk][2] for j, kk in enumerate(P["audit"]) if j != i]
    loo[k] = rho_p(xs, [v for j, v in enumerate(Y) if j != i])
rec("B", "LOO range", abs(loo["talazoparib"][0] - 0.63) < 0.01 and abs(loo["veliparib_prot"][0] + 0.63) < 0.01,
    f"excl tala {loo['talazoparib'][:2]}, excl veli {loo['veliparib_prot'][:2]} (claimed +0.63/-0.63, p=0.5)")

azd = s5val["AZD5305"][2]
cf = {pname: rho_p([s5val[k][2] for k in keys] + [azd], Y + [0.5]) for pname, keys in P.items()}
rec("B", "AZD5305 counterfactual", abs(cf["neutral"][0] + 0.81) <= 0.02 and abs(cf["neutral"][1] - 0.061) <= 0.006
    and abs(cf["audit"][0] + 0.55) <= 0.02 and abs(cf["audit"][1] - 0.272) <= 0.006,
    f"neutral {cf['neutral'][:2]}; audit {cf['audit'][:2]} (claimed -0.81/0.061; -0.55/0.272)")
cft = {pname: rho_p([s5val[k][2] for k in keys] + [azd], Y + [6.0]) for pname, keys in P.items()}
rec("B", "AZD5305 top-placement", abs(cft["audit"][0] - 0.32) <= 0.02 and abs(cft["neutral"][0] - 0.06) <= 0.02,
    f"audit {cft['audit'][0]}, neutral {cft['neutral'][0]} (claimed +0.32 / +0.06)")


def pert(keys, idx, delta):
    xs = [s5val[k][2] for k in keys]
    xs[idx] += delta
    return rho_p(xs, Y)


pn1 = pert(P["neutral"], P["neutral"].index("veliparib"), -6.0)
pn2 = pert(P["neutral"], P["neutral"].index("veliparib"), 2.0)
pa1 = pert(P["audit"], P["audit"].index("veliparib_prot"), -6.0)
pa2 = pert(P["audit"], P["audit"].index("veliparib_prot"), 2.0)
rec("B", "veli perturb neutral", abs(pn1[0] - 0.05) < 0.03 and pn1[1] > 0.99
    and abs(pn2[0] + 0.98) < 0.03 and abs(pn2[1] - 0.033) < 0.006,
    f"-6 {pn1[:2]}; +2 {pn2[:2]} (claimed +0.05/1.000; -0.98/0.033)")
lo = min(pa1[0], pa2[0]); hi = max(pa1[0], pa2[0]); pmin = min(pa1[1], pa2[1])
rec("B", "veli perturb audit", lo >= -0.38 and hi <= 0.23 and pmin >= 0.627,
    f"-6 {pa1[:2]}; +2 {pa2[:2]} (claimed -0.36..+0.21, p>=0.633)")

# --------------------------------------------- lattice checks
r_neg = [float(r["rho"]) for r in lat if r["metric"] == "S2CV1"]
r_s1 = [float(r["rho"]) for r in lat if r["metric"] == "S1"]
r_rt = [float(r["rho"]) for r in lat if r["metric"] == "ratio"]
ps = [float(r["p"]) for r in lat]
rec("A", "lattice S2CV1 neg range", abs(min(r_neg) + 0.667) < 0.006 and abs(max(r_neg) + 0.051) < 0.006,
    f"{min(r_neg):.3f}..{max(r_neg):.3f} (claimed -0.67..-0.05)")
rec("A", "lattice S1 pos range", abs(min(r_s1) - 0.308) < 0.006 and abs(max(r_s1) - 0.616) < 0.006,
    f"{min(r_s1):.3f}..{max(r_s1):.3f} (claimed +0.31..+0.62)")
rec("A", "lattice ratio sign change", min(r_rt) < 0 < max(r_rt),
    f"ratio {min(r_rt):.3f}..{max(r_rt):.3f}")
rec("A", "lattice p floor", abs(min(ps) - 0.133) < 0.006, f"min p {min(ps):.3f} (claimed >=0.133)")
best = max(lat, key=lambda r: abs(float(r["rho"])))
rec("A", "lattice strongest", abs(float(best["rho"]) - 0.821) < 0.006,
    f"{best['panel']}/{best['metric']} {best['rho']}/{best['p']} (claimed +0.82/0.133)")

# --------------------------------------------- p_min enumeration
seq = [5, 4, 2.5, 2.5, 1]
perms_rank = {tuple(midrank(p)) for p in itertools.permutations(seq)}
rec("B", "distinct tied vectors", len(perms_rank) == 60, f"{len(perms_rank)} (claimed 60)")
def corr_point(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb)

ry = midrank(seq)
raw = [abs(corr_point(midrank(p), ry)) for p in itertools.permutations([1, 2, 3, 4, 5])]
max_rho = max(raw)
cnt975 = sum(1 for v in raw if v >= max_rho - 1e-9)
rec("B", "p_min", abs(max_rho - 0.9747) < 1e-3 and cnt975 == 4,
    f"max|rho|={max_rho:.4f}, {cnt975}/120 = {cnt975/120:.4f} (claimed 0.975, 4/120 = 0.033)")

# --------------------------------------------- ordinal claims (manuscript text)
OD = {(r["panel"], r["metric"], float(r["s"])): r for r in ordi}
obad = []


def ocheck(panel, metric, s, key, want, tol):
    r = OD.get((panel, metric, s))
    if r is None:
        obad.append(f"{panel}/{metric}/s={s}: missing")
        return
    if abs(float(r[key]) - want) > tol:
        obad.append(f"{panel}/{metric}/s={s} {key}={r[key]} vs {want}")


ocheck("combined", "S2CV1", 12.0, "P_rho_neg", 0.57, 0.006)
ocheck("combined", "S2CV1", 12.0, "median", -0.10, 0.006)
ocheck("combined", "S2CV1", 12.0, "ci_lo", -0.82, 0.006)
ocheck("combined", "S2CV1", 12.0, "ci_hi", 0.67, 0.006)
band = [float(OD[("combined", "S2CV1", s)]["P_rho_neg"]) for s in (6.0, 12.0, 18.0)]
if not (abs(min(band) - 0.55) < 0.006 and abs(max(band) - 0.63) < 0.006):
    obad.append(f"combined S2CV1 band {band} vs 0.55-0.63")
ocheck("combined", "S1", 12.0, "P_rho_neg", 0.10, 0.006)
ocheck("combined", "AAI", 0.3, "P_rho_neg", 0.44, 0.006)
ocheck("ruca_prot", "S2CV1", 12.0, "P_rho_neg", 0.67, 0.006)
ocheck("ruca_prot", "S1", 12.0, "P_rho_neg", 0.40, 0.006)
ocheck("ruca_prot", "AAI", 0.3, "P_rho_neg", 0.30, 0.006)
banded = [float(r["P_rho_neg"]) for r in ordi if float(r["s"]) != 0.0]
if abs(max(banded) - 0.78) > 0.006:
    obad.append(f"banded max {max(banded)} vs 0.78")
tilt = [float(r["P_rho_neg_arcsine_tilt"]) for r in ordi if float(r["s"]) != 0.0]
if not (abs(min(tilt) - 0.01) < 0.006 and abs(max(tilt) - 0.81) < 0.006):
    obad.append(f"tilt range {min(tilt)}-{max(tilt)} vs 0.01-0.81")
ocheck("combined", "S1", 6.0, "P_rho_neg", 0.01, 0.004)
ocheck("combined", "S1", 6.0, "ci_lo", 0.21, 0.006)
ocheck("combined", "S1", 6.0, "ci_hi", 0.87, 0.006)
ocheck("combined", "S1", 12.0, "P_rho_neg", 0.10, 0.006)
ocheck("combined", "S1", 18.0, "P_rho_neg", 0.19, 0.006)
oth = [float(OD[(p, "S1", s)]["P_rho_neg"]) for p in ("neutral", "nira_prot", "ruca_prot")
       for s in (6.0, 12.0, 18.0)]
if not (abs(min(oth) - 0.18) < 0.006 and abs(max(oth) - 0.43) < 0.006):
    obad.append(f"other-panel S1 range {min(oth):.2f}-{max(oth):.2f} vs 0.18-0.43")
rec("A", "ordinal claims", not obad, str(obad[:4]) or "all 14 sub-claims match r6d")

# --------------------------------------------- cross-dock + S1 seeds
xd_ok = all(abs(crossdock[k] - v) < 0.005 for k, v in
            [("rucaparib", 0.58), ("olaparib", 2.46), ("niraparib", 0.91),
             ("talazoparib", 1.45), ("veliparib", 1.32)])
rec("A", "crossdock RMSDs", xd_ok, f"{crossdock}")
rec("A", "S1 seed SDs", s1seeds.get("talazoparib", ("", "", ""))[2] == "42.2±26.6"
    and s1seeds.get("veliparib", ("", "", ""))[2] == "40.3±6.1", f"{s1seeds}")

# --------------------------------------------- S1 13-system ranges
s2blk = si[si.find("Table S2: Cumulant"):si.find("\\end{table}", si.find("Table S2: Cumulant"))]
s1vals = {}
ppat = re.compile(r"^([A-Za-z0-9 .()]+?)\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&")
for line in s2blk.splitlines():
    m = ppat.match(line)
    if m:
        nm = LIG.get(m.group(1).strip())
        if nm:
            s1vals[nm] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
ext = {r["ligand"].lower(): (float(r["C1"]), float(r["C2"]), float(r["C3"]))
       for r in csv.DictReader(open(DATA / "01_curated/extension_s1_wells.csv"))}
alls1 = dict(s1vals)
alls1.update(ext)
rec("A", "S1 table assembly", len(alls1) == 13, sorted(alls1))
c3s = [v[2] for v in alls1.values()]
sps = [stdev(v) for v in alls1.values()]
rec("A", "S1 C3 range", abs(min(c3s) - 12.5) < 0.06 and abs(max(c3s) - 59.2) < 0.06,
    f"{min(c3s):.1f}--{max(c3s):.1f} (claimed 12.5--59.2)")
rec("A", "S1 spread range", abs(min(sps) - 4.7) < 0.15 and abs(max(sps) - 28.6) < 0.06,
    f"{min(sps):.1f}--{max(sps):.1f} (claimed 4.7--28.6)")

# --------------------------------------------- matrix cells
matrix_claims = [
    ("audit S2CV1 mm", ("combined", "S2CV1", "mm"), (-0.21, 0.767)),
    ("audit S2CV1 q90", ("combined", "S2CV1", "q90"), (-0.36, 0.633)),
    ("audit S2CV1 t95", ("combined", "S2CV1", "t95"), (-0.21, 0.767)),
    ("audit S2CV2 mm", ("combined", "S2CV2", "mm"), (0.41, 0.500)),
    ("audit S2CV2 q90", ("combined", "S2CV2", "q90"), (0.67, 0.267)),
    ("audit S2CV2 t95", ("combined", "S2CV2", "t95"), (0.56, 0.367)),
    ("audit S1 mm", ("combined", "S1", "mm"), (0.62, 0.300)),
    ("audit ratio mm", ("combined", "AAI", "mm"), (-0.10, 0.900)),
    ("audit ratio q90", ("combined", "AAI", "q90"), (-0.21, 0.767)),
    ("audit ratio t95", ("combined", "AAI", "t95"), (0.46, 0.433)),
    ("neutral S2CV1 mm", ("neutral", "S2CV1", "mm"), (-0.667, 0.267)),
    ("neutral S2CV1 q90", ("neutral", "S2CV1", "q90"), (-0.359, 0.633)),
    ("neutral S1 mm", ("neutral", "S1", "mm"), (0.308, 0.667)),
    ("neutral ratio mm", ("neutral", "AAI", "mm"), (0.205, 0.767)),
]
mbad = []
for cid, key, exp in matrix_claims:
    got = estim.get(key)
    if got is None or abs(got[0] - exp[0]) > 0.006 or abs(got[1] - exp[1]) > 0.006:
        mbad.append(f"{cid} {got} vs {exp}")
rec("A", "estimatrix cells", not mbad, str(mbad[:4]) or "14/14")

# --------------------------------------------- tex presence (C tier)
tex_checks = [
    ("abstract rho", ms, "$\\rho = -0.67$ (exact $p = 0.267$"),
    ("abstract audit rho", ms, "$\\rho = -0.21$ ($p = 0.767$)"),
    ("abstract pmin", ms, "four of 120 data permutations"),
    ("abstract window band", ms, "4.6--18.4~kcal/mol"),
    ("abstract separation band", ms, "2.6--5.3~kcal/mol"),
    ("four-panel ratio", ms, "span ratio ranging from $-0.10$ to $+0.31$ ($p \\geq 0.667$)"),
    ("alt ordering", ms, "$\\rho = -0.45$ (exact two-tailed $p = 0.600$)"),
    ("four-panel S1", ms, "$+0.31$ to $+0.62$; $p \\geq 0.300$"),
    ("lattice sign", ms, "change sign with the panel"),
    ("lattice p floor", ms, "all $p \\geq 0.133$"),
    ("lattice strongest", ms, "for the S2 CV2 span under veliparib protonation alone"),
    ("estimatrix q90", ms, "$-0.36$ ($p = 0.633$) for the 5--95\\% quantile"),
    ("estimatrix t95", ms, "$-0.21$ ($p = 0.767$) for the 95\\% truncated range"),
    ("estimatrix ratio", ms, "$-0.21$ to $+0.46$ ($p \\geq 0.433$)"),
    ("S10 gap", si, "15.28 versus 14.73~kcal/mol"),
    ("S10 pct phrase", si + ms, "35--45\\% of the respective means"),
    ("power curve", ms, "serves as the corresponding power curve"),
    ("S17 caption", si, "Table S17: Cross-seed stability"),
    ("pc 32pct", ms, "recovery reaches 32\\% at the largest within-band separations"),
    ("pc sigma8", ms, "$\\sigma \\geq 8$~kcal/mol"),
    ("pc sigma4.6", ms, "$\\sigma = 4.6$~kcal/mol"),
    ("frames 25858", si, "25{,}858"),
    ("frames 521", si, "521 records"),
    ("S12 AZD span", si, "49.5~kcal/mol) inserted as a hypothetical sixth"),
]
for cid, tex, s in tex_checks:
    rec("C", cid, s in tex, s[:60])

# --------------------------------------------- positive control cells
pc_claims = [
    ((2.0, 2.0), 0.39, 0.03), ((4.0, 2.0), 0.87, 0.03), ((6.0, 2.0), 0.98, 0.02),
    ((2.0, 4.6), 0.10, 0.03), ((4.0, 4.6), 0.32, 0.03), ((10.0, 4.6), 0.90, 0.03),
    ((2.0, 8.0), 0.07, 0.03), ((4.0, 8.0), 0.14, 0.03),
    ((10.0, 8.0), 0.53, 0.03), ((20.0, 8.0), 0.95, 0.02),
    ((20.0, 15.0), 0.61, 0.04), ((20.0, 18.4), 0.45, 0.04),
]
pbad = []
for key, want, tol in pc_claims:
    got = pc.get(key)
    if got is None or abs(got - want) > tol:
        pbad.append(f"{key}: {got} vs {want}")
n0 = [v for (d, s), v in pc.items() if d == 0.0]
rec("A", "pc quoted cells", not pbad,
    str(pbad[:4]) or f"12/12; Delta=0 range {min(n0)*100:.1f}-{max(n0)*100:.1f}%")

# ---------------------------------------------------------------- report
fails = [r for r in results if r[2] != "PASS"]
print(f"{'TIER':4s} {'CHECK':34s} STATUS  DETAIL")
for t, cid, stt, det in results:
    print(f"{t:4s} {cid:34s} {stt:7s} {det[:100]}")
print(f"\n{len(results) - len(fails)}/{len(results)} checks pass")
sys.exit(0 if not fails else 1)
