#!/usr/bin/env python3
"""s63_harmonized_matrix_lattice.py — R4: harmonized-convention estimator
matrix and chemical-state lattice (50-ps cadence, common 14-ns window).

Builds, for the eight S2 series (five neutral + three protonated), harmonized
CV1 and CV2 spans under three estimators (min--max C3, 5--95% quantile, 95%
truncated), combines them with the archival S1 spans and the harmonized span
ratio, and evaluates:

  (A) the 4-panel x 4-metric x 3-estimator matrix (exact tie-aware permutation)
  (B) the complete 8-combination chemical-state lattice for four metrics

Outputs
  data/analysis/r6/s63_harmonized_matrix.csv
  data/analysis/r6/s63_harmonized_lattice.csv
  data/per_frame/veli_prot_s2_frames.csv.gz   (re-extracted with CV2)
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
NWIN = 280


def spans_all(cv, dv):
    data, dV = np.asarray(cv, float), np.asarray(dv, float)
    mx = DISC * (int(np.amax(data) / DISC) + 1)
    mn = DISC * (int(np.amin(data) / DISC) - 1)
    binsX = mod.assignbins([mn, mx], DISC)
    with contextlib.redirect_stdout(io.StringIO()):
        hist, _e, c1, c2, c3 = mod.reweight_CE(data, CUT, binsX, DISC, dV, T, False)
    c1, c2, c3 = -c1 / BETA, -c2 / BETA, -c3 / BETA
    pmf = mod.prephist(hist, T, EMAX)
    p3 = mod.normalize(pmf + c1 + c2 + c3, EMAX)
    valid = hist >= CUT
    vals = p3[valid]
    return dict(mm=float(vals.max() - vals.min()),
                q90=float(np.quantile(vals, 0.95) - np.quantile(vals, 0.05)),
                t95=float(np.percentile(vals, 97.5) - np.percentile(vals, 2.5)))


def _read_table(f):
    cols = f.readline().strip().split(",")
    body = f.read()
    first = body.split("\n", 1)[0]
    delim = "," if "," in first else None
    arr = np.loadtxt(io.StringIO(body), delimiter=delim)
    return cols, arr


def load_arch(name, cvcol="cv1_prot_dna_com"):
    with gzip.open(ROOT / "data/per_frame" / f"{name}_frames.csv.gz", "rt") as f:
        cols, arr = _read_table(f)
    return arr[:, cols.index(cvcol)], arr[:, cols.index("dV_D")]


def extract_veli_prot():
    cache = ROOT / "data/per_frame/veli_prot_s2_frames.csv.gz"
    if cache.exists():
        with gzip.open(cache, "rt") as f:
            hdr = f.readline().strip()
            if "cv2_hd_art_com" in hdr:
                body = f.read()
                delim = "," if "," in body.split("\n", 1)[0] else None
                arr = np.loadtxt(io.StringIO(body), delimiter=delim)
                return arr[:, 3], arr[:, 4], arr[:, 2]
    import MDAnalysis as mda
    d = DATA / "veli_prot/sys2_veliparib_prot"
    u = mda.Universe(str(d / "veliparib_best_prot.prmtop"), str(d / "output.dcd"))
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    hd = u.select_atoms("resid 662:787 and name CA")
    art = u.select_atoms("resid 788:1011 and name CA")
    log = np.loadtxt(d / "gamd.log", comments="#")
    step_dv = dict(zip(log[:, 1].astype(np.int64), log[:, 7]))
    rows = []
    for k, ts in enumerate(u.trajectory, start=1):
        step = 25000 * k
        dv = step_dv.get(step)
        if dv is None:
            continue
        cv1 = float(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
        cv2 = float(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
        rows.append((k, step, dv, cv1, cv2))
    with gzip.open(ROOT / "data/per_frame/veli_prot_s2_frames.csv.gz", "wt") as f:
        f.write("frame,step_ps,dV_D,cv1_prot_dna_com,cv2_hd_art_com\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]:.6f},{r[3]:.6f},{r[4]:.6f}\n")
    arr = np.array([[r[2], r[3], r[4]] for r in rows])
    return arr[:, 1], arr[:, 2], arr[:, 0]


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


# ---------- build harmonized spans ----------
H = {}   # system -> dict(CV1=..., CV2=...)
SYSFILES = {
    "tala": ("orig_s2_talazoparib", 50), "nira": ("sys2_niraparib", 1),
    "ola": ("sys2_olaparib", 1), "ruca": ("sys2_rucaparib", 1),
    "veli": ("orig_s2_veliparib", 1),
    "nira_prot": ("sys2_nira_prot", 1), "ruca_prot": ("sys2_ruca_prot", 1),
}
for key, (name, stride) in SYSFILES.items():
    cv1, dv = load_arch(name, "cv1_prot_dna_com")
    cv2, _ = load_arch(name, "cv2_hd_art_com")
    if stride > 1:
        cv1, cv2, dv = cv1[::stride], cv2[::stride], dv[::stride]
    H[key] = dict(CV1=spans_all(cv1[:NWIN], dv[:NWIN]),
                  CV2=spans_all(cv2[:NWIN], dv[:NWIN]))

cv1, cv2, dv = extract_veli_prot()
H["veli_prot"] = dict(CV1=spans_all(cv1[:NWIN], dv[:NWIN]),
                      CV2=spans_all(cv2[:NWIN], dv[:NWIN]))

print("system      CV1(mm/q90/t95)            CV2(mm/q90/t95)")
for k, v in H.items():
    print(f"{k:10s} {v['CV1']['mm']:6.2f} {v['CV1']['q90']:6.2f} {v['CV1']['t95']:6.2f}   "
          f"{v['CV2']['mm']:6.2f} {v['CV2']['q90']:6.2f} {v['CV2']['t95']:6.2f}")

# archival S1 C3 spans (protocol-relative; unchanged)
S1 = {"tala": 42.2, "nira": 50.8, "nira_prot": 59.2, "ola": 55.4,
      "ruca": 27.6, "ruca_prot": 27.7, "veli": 40.3, "veli_prot": 23.7}

PANELS = {
    "neutral": dict(tala="tala", nira="nira", ola="ola", ruca="ruca", veli="veli",
                    s1={"tala": "tala", "nira": "nira", "ola": "ola",
                        "ruca": "ruca", "veli": "veli"}),
    "niraparib-protonated": dict(tala="tala", nira="nira_prot", ola="ola",
                                 ruca="ruca", veli="veli",
                                 s1={"tala": "tala", "nira": "nira_prot", "ola": "ola",
                                     "ruca": "ruca", "veli": "veli"}),
    "rucaparib-protonated": dict(tala="tala", nira="nira", ola="ola",
                                 ruca="ruca_prot", veli="veli",
                                 s1={"tala": "tala", "nira": "nira", "ola": "ola",
                                     "ruca": "ruca_prot", "veli": "veli"}),
    "fully protonated": dict(tala="tala", nira="nira_prot", ola="ola",
                             ruca="ruca_prot", veli="veli_prot",
                             s1={"tala": "tala", "nira": "nira_prot", "ola": "ola",
                                 "ruca": "ruca_prot", "veli": "veli_prot"}),
}
Y = [5.0, 4.0, 2.5, 2.5, 1.0]
ORDER = ["tala", "nira", "ola", "ruca", "veli"]

matrix_rows = []
print("\n── harmonized matrix (rho / exact p) ──")
for pname, p in PANELS.items():
    for metric in ("S2CV1", "S2CV2", "S1", "ratio"):
        for est in ("mm", "q90", "t95"):
            if metric in ("S2CV1", "S2CV2"):
                hkey = "CV1" if metric == "S2CV1" else "CV2"
                vals = [H[p[k]][hkey][est] for k in ORDER]
            elif metric == "S1":
                vals = [S1[p["s1"][k]] for k in ORDER]
            else:
                vals = [S1[p["s1"][k]] / H[p[k]]["CV2"][est] for k in ORDER]
            r, cnt = rho_p(vals, Y)
            matrix_rows.append((pname, metric, est, r, cnt, cnt / 120.0))
            print(f"{pname:22s} {metric:6s} {est:4s} rho={r:+.3f} p={cnt:3d}/120")
sig = [x for x in matrix_rows if x[5] < 0.05]
print(f"\nwidth: {len(matrix_rows)} cells; significant (<0.05): {len(sig)}")

lattice_rows = []
print("\n── harmonized lattice (8 state combinations x 4 metrics) ──")
for sn, sr, sv in itertools.product((0, 1), repeat=3):
    state = f"{sn}{sr}{sv}"
    for metric in ("S2CV1", "S2CV2", "S1", "ratio"):
        if metric == "S2CV1":
            vals = [H["tala"]["CV1"]["mm"], H["nira_prot" if sn else "nira"]["CV1"]["mm"],
                    H["ola"]["CV1"]["mm"], H["ruca_prot" if sr else "ruca"]["CV1"]["mm"],
                    H["veli_prot" if sv else "veli"]["CV1"]["mm"]]
        elif metric == "S2CV2":
            vals = [H["tala"]["CV2"]["mm"], H["nira_prot" if sn else "nira"]["CV2"]["mm"],
                    H["ola"]["CV2"]["mm"], H["ruca_prot" if sr else "ruca"]["CV2"]["mm"],
                    H["veli_prot" if sv else "veli"]["CV2"]["mm"]]
        elif metric == "S1":
            vals = [S1["tala"], S1["nira_prot" if sn else "nira"], S1["ola"],
                    S1["ruca_prot" if sr else "ruca"], S1["veli_prot" if sv else "veli"]]
        else:
            vals = [S1["tala"] / H["tala"]["CV2"]["mm"],
                    S1["nira_prot" if sn else "nira"] / H["nira_prot" if sn else "nira"]["CV2"]["mm"],
                    S1["ola"] / H["ola"]["CV2"]["mm"],
                    S1["ruca_prot" if sr else "ruca"] / H["ruca_prot" if sr else "ruca"]["CV2"]["mm"],
                    S1["veli_prot" if sv else "veli"] / H["veli_prot" if sv else "veli"]["CV2"]["mm"]]
        r, cnt = rho_p(vals, Y)
        lattice_rows.append((state, metric, r, cnt, cnt / 120.0))

cv1 = [x[2] for x in lattice_rows if x[1] == "S2CV1"]
s1 = [x[2] for x in lattice_rows if x[1] == "S1"]
ps = [x[4] for x in lattice_rows]
print(f"S2CV1 lattice range: [{min(cv1):+.3f}, {max(cv1):+.3f}]")
print(f"S1 lattice range:    [{min(s1):+.3f}, {max(s1):+.3f}]")
print(f"min exact p across 32 cells: {min(ps):.3f}")
print(f"max |rho| cell: {max(lattice_rows, key=lambda x: abs(x[2]) )}")

with open(ROOT / "data/analysis/r6/s63_harmonized_spans.csv", "w") as f:
    f.write("system,metric,estimator,span_kcal\n")
    for k, v in H.items():
        for m in ("CV1", "CV2"):
            for e in ("mm", "q90", "t95"):
                f.write(f"{k},{m},{e},{v[m][e]:.4f}\n")

with open(ROOT / "data/analysis/r6/s63_harmonized_matrix.csv", "w") as f:
    f.write("panel,metric,estimator,rho,p_exact,cnt120\n")
    for r in matrix_rows:
        f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[5]:.4f},{r[4]}\n")
with open(ROOT / "data/analysis/r6/s63_harmonized_lattice.csv", "w") as f:
    f.write("state,metric,rho,p_exact,cnt120\n")
    for r in lattice_rows:
        f.write(f"{r[0]},{r[1]},{r[2]},{r[4]:.4f},{r[3]}\n")
print("\nwrote s63_harmonized_matrix.csv / s63_harmonized_lattice.csv")
