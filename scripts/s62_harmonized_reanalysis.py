#!/usr/bin/env python3
"""s62_harmonized_reanalysis.py — Q-max step 1: dual harmonization
(frame cadence 50 ps AND matched 14-ns production window) for the S2 panel.

Systems and sources
  talazoparib   orig_s2_talazoparib (1 ps production archive -> every 50th frame)
  niraparib     sys2_niraparib (50 ps)
  olaparib      sys2_olaparib
  rucaparib     sys2_rucaparib
  veliparib     orig_s2_veliparib (50 ps, 517 rows)
  nira_prot     sys2_nira_prot
  ruca_prot     sys2_ruca_prot
  veli_prot     extracted here from DATA_ROOT/veli_prot/sys2_veliparib_prot (50 ps)

Harmonized window: first 280 production frames at 50 ps (= 14 ns, the minimum
production length in the panel; talazoparib). Every system is truncated to this
window; nothing is extrapolated.

Outputs
  data/analysis/r6/s62_harmonized_summary.csv     spans (C1/C2/C3, q90, t95),
                                                  2-ns window SD, block bootstrap
  data/analysis/r6/s62_harmonized_association.csv neutral / fully-protonated
                                                  panels under the harmonized
                                                  convention + seed-mean variants
  data/analysis/r6/s62_82_2_diagnosis.csv         edge-bin diagnosis of the
                                                  seed-123 protonated span
  data/per_frame/veli_prot_s2_frames.csv.gz       extracted series
"""
import contextlib
import gzip
import importlib.util
import io
import itertools
import math
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("DATA_ROOT", "/Volumes/tjogzt4T/PARPi_data"))
spec = importlib.util.spec_from_file_location(
    "pyrew1d", ROOT / "tools/PyReweighting/PyReweighting-1D.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

BETA = 1.0 / (0.001987 * 300.0)
DISC, EMAX, CUT, T = 0.1, 20.0, 2, 300.0
RNG = np.random.default_rng(49)
NWIN = 280           # first 280 frames at 50 ps = 14 ns
PER_WIN = 40         # 2 ns


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
        valid = hist >= CUT
        vals = p[valid]
        out[order] = float(vals.max() - vals.min())
        if order == "C3":
            out["q90"] = float(np.quantile(vals, 0.95) - np.quantile(vals, 0.05))
            out["t95"] = float(np.percentile(vals, 97.5) - np.percentile(vals, 2.5))
    out["n_valid_bins"] = int((hist >= CUT).sum())
    out["n_empty"] = int((hist == 0).sum())
    out["n_single"] = int((hist == 1).sum())
    return out


def load_arch(name):
    with gzip.open(ROOT / "data/per_frame" / f"{name}_frames.csv.gz", "rt") as f:
        cols = f.readline().strip().split(",")
        arr = np.loadtxt(f)
    return arr[:, cols.index("cv1_prot_dna_com")], arr[:, cols.index("dV_D")]


def extract_veli_prot():
    import MDAnalysis as mda
    d = DATA / "veli_prot/sys2_veliparib_prot"
    out_csv = ROOT / "data/per_frame/veli_prot_s2_frames.csv.gz"
    if out_csv.exists():
        cv, dv = load_arch("veli_prot_s2")
        return cv, dv
    u = mda.Universe(str(d / "veliparib_best_prot.prmtop"), str(d / "output.dcd"))
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    log = np.loadtxt(d / "gamd.log", comments="#")
    step_dv = dict(zip(log[:, 1].astype(np.int64), log[:, 7]))
    rows = []
    for k, ts in enumerate(u.trajectory, start=1):
        step = 25000 * k
        dv = step_dv.get(step)
        if dv is None:
            continue
        cv1 = float(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
        rows.append((k, step, dv, cv1))
    with gzip.open(out_csv, "wt") as f:
        f.write("frame,step_ps,dV_D,cv1_prot_dna_com\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]:.6f},{r[3]:.6f}\n")
    arr = np.array([[r[2], r[3]] for r in rows])
    return arr[:, 1], arr[:, 0]


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
    mx, my = sum(rx) / n, sum(ry) / n
    sxx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    syy = math.sqrt(sum((b - my) ** 2 for b in ry))

    def r_of(seq):
        m2 = sum(seq) / n
        s2 = math.sqrt(sum((a - m2) ** 2 for a in seq))
        return sum((a - m2) * (b - my) for a, b in zip(seq, ry)) / (s2 * syy)

    obs = r_of(rx)
    cnt = 0
    for perm in itertools.permutations(xs):
        if max(perm) == min(perm):
            continue
        if abs(r_of(midrank(perm))) >= abs(obs) - 1e-12:
            cnt += 1
    return round(obs, 3), cnt


# ---------- assemble harmonized series ----------
SYS = {}
for name in ["sys2_niraparib", "sys2_olaparib", "sys2_rucaparib",
             "sys2_nira_prot", "sys2_ruca_prot", "orig_s2_veliparib",
             "s10_neu_s7", "s10_neu_s123", "s10_prot_s7", "s10_prot_s123"]:
    cv, dv = load_arch(name)
    SYS[name] = (cv, dv)

# talazoparib: 1-ps production archive -> 50-ps
cv, dv = load_arch("orig_s2_talazoparib")
SYS["orig_s2_talazoparib_50ps"] = (cv[::50], dv[::50])

cv, dv = extract_veli_prot()
SYS["veli_prot_s2"] = (cv, dv)

print("── harmonized spans (first 280 frames at 50 ps = 14 ns) ──")
print(f"{'system':26s} {'n':>4s} {'C1':>6s} {'C2':>6s} {'C3':>6s} {'q90':>6s} {'t95':>6s} "
      f"{'winSD':>6s} {'nv':>4s} {'n0':>3s} {'n1':>3s} {'bootSD':>6s}")
summary = []
for name, (cv, dv) in SYS.items():
    cv, dv = cv[:NWIN], dv[:NWIN]
    s = spans(cv, dv)
    ws = [spans(cv[k:k + PER_WIN], dv[k:k + PER_WIN])["C3"]
          for k in range(0, len(cv) - PER_WIN + 1, PER_WIN)]
    wsd = float(np.std(ws)) if len(ws) > 1 else float("nan")
    nb = int(np.ceil(len(cv) / PER_WIN))
    boots = []
    for _ in range(400):
        st = RNG.integers(0, len(cv), nb)
        idx = np.concatenate([(t + np.arange(PER_WIN)) % len(cv) for t in st])[:len(cv)]
        boots.append(spans(cv[idx], dv[idx])["C3"])
    boots = np.array(boots)
    summary.append((name, len(cv), s["C1"], s["C2"], s["C3"], s["q90"], s["t95"],
                    wsd, s["n_valid_bins"], s["n_empty"], s["n_single"],
                    float(boots.std(ddof=1)), float(boots.mean()),
                    float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))))
    print(f"{name:26s} {len(cv):>4d} {s['C1']:6.2f} {s['C2']:6.2f} {s['C3']:6.2f} "
          f"{s['q90']:6.2f} {s['t95']:6.2f} {wsd:6.2f} {s['n_valid_bins']:>4d} "
          f"{s['n_empty']:>3d} {s['n_single']:>3d} {boots.std(ddof=1):6.2f}", flush=True)

# ---------- associations ----------
Y = [5.0, 4.0, 2.5, 2.5, 1.0]
V = {row[0]: row[4] for row in summary}
panels = {
    "harmonized neutral (tala, nira, ola, ruca, veli)":
        [V["orig_s2_talazoparib_50ps"], V["sys2_niraparib"], V["sys2_olaparib"],
         V["sys2_rucaparib"], V["orig_s2_veliparib"]],
    "harmonized fully protonated":
        [V["orig_s2_talazoparib_50ps"], V["sys2_nira_prot"], V["sys2_olaparib"],
         V["sys2_ruca_prot"], V["veli_prot_s2"]],
    "harmonized neutral, nira = 3-seed mean":
        [V["orig_s2_talazoparib_50ps"],
         float(np.mean([V["sys2_niraparib"], V["s10_neu_s7"], V["s10_neu_s123"]])),
         V["sys2_olaparib"], V["sys2_rucaparib"], V["orig_s2_veliparib"]],
    "harmonized protonated, nira = 3-seed mean":
        [V["orig_s2_talazoparib_50ps"],
         float(np.mean([V["sys2_nira_prot"], V["s10_prot_s7"], V["s10_prot_s123"]])),
         V["sys2_olaparib"], V["sys2_ruca_prot"], V["veli_prot_s2"]],
}
print("\n── exact permutation under the harmonized convention ──")
assoc = []
for pname, xs in panels.items():
    r, cnt = rho_p(xs, Y)
    assoc.append((pname, r, cnt, cnt / 120.0))
    print(f"{pname:52s} rho={r:+.3f} p={cnt}/120={cnt/120:.3f}")

# ---------- 82.2 edge-bin diagnosis (seed123 protonated, full 26 ns) ----------
print("\n── 82.2 diagnosis (s10_prot_s123) ──")
cv, dv = load_arch("s10_prot_s123")
data, dV = np.asarray(cv, float), np.asarray(dv, float)
mx = DISC * (int(np.amax(data) / DISC) + 1)
mn = DISC * (int(np.amin(data) / DISC) - 1)
binsX = mod.assignbins([mn, mx], DISC)
with contextlib.redirect_stdout(io.StringIO()):
    hist, _e, c1, c2, c3 = mod.reweight_CE(data, CUT, binsX, DISC, dV, T, False)
c1, c2, c3 = -c1 / BETA, -c2 / BETA, -c3 / BETA
pmf = mod.prephist(hist, T, EMAX)
p3 = mod.normalize(pmf + c1 + c2 + c3, EMAX)
valid = (hist >= CUT)
argmax_bin = int(np.argmax(np.where(valid, p3, -np.inf)))
argmin_bin = int(np.argmin(np.where(valid, p3, np.inf)))
diag = {
    "span_C3": float(p3[valid].max() - p3[valid].min()),
    "nbins": int(len(hist)), "n_valid": int(valid.sum()),
    "max_bin_index": argmax_bin, "max_bin_occ": int(hist[argmax_bin]),
    "max_bin_edge": bool(argmax_bin in (np.flatnonzero(valid)[0], np.flatnonzero(valid)[-1])),
    "min_bin_index": argmin_bin, "min_bin_occ": int(hist[argmin_bin]),
    "min_bin_edge": bool(argmin_bin in (np.flatnonzero(valid)[0], np.flatnonzero(valid)[-1])),
}
# leave-one-extreme-bin-out
for label, drop in (("drop_max_bin", argmax_bin), ("drop_min_bin", argmin_bin)):
    h2 = hist.copy()
    h2[drop] = 0
    # recompute weights excluding frames in that bin: approximate by occupancy check
    in_bin = (np.floor((data - binsX[0]) / DISC).astype(int) == drop)
    if (~in_bin).sum() > 50:
        keep = ~in_bin
        d2, v2 = data[keep], dV[keep]
        m2x = DISC * (int(np.amax(d2) / DISC) + 1)
        m2n = DISC * (int(np.amin(d2) / DISC) - 1)
        with contextlib.redirect_stdout(io.StringIO()):
            h, _e2, a1, a2, a3 = mod.reweight_CE(d2, CUT, mod.assignbins([m2n, m2x], DISC), DISC, v2, T, False)
        pm2 = mod.prephist(h, T, EMAX)
        pp = mod.normalize(pm2 + a1 + a2 + a3, EMAX)
        vv = h >= CUT
        diag[label + "_span"] = float(pp[vv].max() - pp[vv].min())
        diag[label + "_n_frames_removed"] = int(in_bin.sum())
print(diag)

# ---------- write ----------
out1 = ROOT / "data/analysis/r6/s62_harmonized_summary.csv"
with open(out1, "w") as f:
    f.write("system,n_frames,C1,C2,C3,q90,t95,win_sd_2ns,n_valid_bins,n_empty,n_single,"
            "boot_sd,boot_mean,boot_q025,boot_q975\n")
    for r in summary:
        f.write(",".join(str(x) for x in r) + "\n")
out2 = ROOT / "data/analysis/r6/s62_harmonized_association.csv"
with open(out2, "w") as f:
    f.write("panel,rho,p_exact,cnt120\n")
    for r in assoc:
        f.write(f"{r[0]},{r[1]},{r[3]:.4f},{r[2]}\n")
out3 = ROOT / "data/analysis/r6/s62_82_2_diagnosis.csv"
with open(out3, "w") as f:
    for k, v in diag.items():
        f.write(f"{k},{v}\n")
print(f"\nwrote {out1.name}, {out2.name}, {out3.name}")
