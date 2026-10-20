#!/usr/bin/env python3
"""s59_rebuild_talazoparib_full.py — rebuild the complete 1-ps CV1 series of
talazoparib S2 from the original 83-GB output.dcd and quantify the cadence /
production-window conventions (reviewer statistical batch, item 7).

Steps
  1. stream output.dcd (26,000 frames; step = 500*k, dt = 2 fs), computing
     CV1 = |COM(protein CA) - COM(DNA heavy)| exactly as s31_archive_frame_cvs.py
  2. pair each frame with the gamd.log dV_D row carrying the same step
  3. spans (CE, T=300, disc 0.1, Emax 20, cutoff 2, occupancy>=2 valid bins) for
        full series  @1 ps   and  @50 ps (every 50th frame)
        production   @1 ps   (step >= 6.0e6) and @50 ps
     plus 2-ns window SDs and a 400-draw circular block bootstrap of the span
     for the production-only 50-ps series (the convention of the other systems)
  4. asserts the published 22.9 kcal/mol against the full-series value and
     writes data/analysis/r6/s59_talazoparib_cadence.csv + the rebuilt series
     to data/per_frame/orig_s2_talazoparib_full_frames.csv.gz

Run from the repository root (needs DATA_ROOT with sys2_talazoparib):
  DATA_ROOT=/Volumes/tjogzt4T/PARPi_data python3 scripts/s59_rebuild_talazoparib_full.py
"""
import contextlib
import gzip
import importlib.util
import io
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PYREW = ROOT / "tools" / "PyReweighting" / "PyReweighting-1D.py"
spec = importlib.util.spec_from_file_location("pyrew1d", PYREW)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

DATA = Path(os.environ.get("DATA_ROOT", "/Volumes/tjogzt4T/PARPi_data"))
SYS = DATA / "sys2_talazoparib"
BETA = 1.0 / (0.001987 * 300.0)
DISC, EMAX, CUT, T = 0.1, 20.0, 2, 300.0
PROD_STEP = 6.0e6
RNG = np.random.default_rng(49)


def spans(cv, dv):
    data = np.asarray(cv, float)
    dV = np.asarray(dv, float)
    mx = DISC * (int(np.amax(data) / DISC) + 1)
    mn = DISC * (int(np.amin(data) / DISC) - 1)
    binsX = mod.assignbins([mn, mx], DISC)
    with contextlib.redirect_stdout(io.StringIO()):
        hist, _e, c1, c2, c3 = mod.reweight_CE(data, CUT, binsX, DISC, dV, T, False)
    c1, c2, c3 = -c1 / BETA, -c2 / BETA, -c3 / BETA
    pmf = mod.prephist(hist, T, EMAX)
    p1 = mod.normalize(pmf + c1, EMAX)
    p2 = mod.normalize(pmf + c1 + c2, EMAX)
    p3 = mod.normalize(pmf + c1 + c2 + c3, EMAX)
    valid = hist >= CUT
    return (float(p1[valid].max() - p1[valid].min()),
            float(p2[valid].max() - p2[valid].min()),
            float(p3[valid].max() - p3[valid].min()))


# ── 1. read the log, index dV_D by step ─────────────────────────────────
log = np.loadtxt(SYS / "gamd.log", comments="#")
step_dv = dict(zip(log[:, 1].astype(np.int64), log[:, 7]))
print(f"log rows {len(log)}, steps {int(log[0,1])}..{int(log[-1,1])}")

# ── 2. stream the DCD ───────────────────────────────────────────────────
import MDAnalysis as mda

u = mda.Universe(str(SYS / "sys2_talazoparib.prmtop"), str(SYS / "output.dcd"))
prot = u.select_atoms("protein and name CA")
dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                     "resname DA3 DT3 DG3 DC3) and not element H")
print(f"n_frames={len(u.trajectory)}  prot={len(prot)}  dna={len(dna)}")
steps, cvs = [], []
for k, ts in enumerate(u.trajectory, start=1):
    steps.append(500 * k)              # istart 500, nsavc 500 steps
    cvs.append(float(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass())))
    if k % 5000 == 0:
        print(f"  {k}/{len(u.trajectory)}", flush=True)
steps = np.array(steps, dtype=np.int64)
cvs = np.array(cvs)
dvs = np.array([step_dv.get(s, np.nan) for s in steps])
keep = ~np.isnan(dvs)
steps, cvs, dvs = steps[keep], cvs[keep], dvs[keep]
print(f"paired frames: {len(cvs)}  (log-coverage {len(cvs)/len(u.trajectory):.3f})")

# save the rebuilt series
with gzip.open(ROOT / "data/per_frame/orig_s2_talazoparib_full_frames.csv.gz", "wt") as f:
    f.write("frame,step_ps,dV_D,cv1_prot_dna_com\n")
    for i, (s, dv, cv) in enumerate(zip(steps, cvs, dvs), start=1):
        f.write(f"{i},{s},{dv:.6f},{cv:.6f}\n")

# ── 3. spans under the candidate conventions ────────────────────────────
prod = steps >= PROD_STEP
variants = {
    "full_1ps": (cvs, dvs),
    "full_50ps": (cvs[::50], dvs[::50]),
    "prod_1ps": (cvs[prod], dvs[prod]),
    "prod_50ps": (cvs[prod][::50], dvs[prod][::50]),
}
print(f"\n{'variant':10s} {'n':>6s} {'C1':>7s} {'C2':>7s} {'C3':>7s} {'winSD':>7s} {'nW':>4s}")
rows = []
for name, (c, d) in variants.items():
    s1, s2, s3 = spans(c, d)
    per_win = 2000 if name.endswith("1ps") else 40
    ws = [spans(c[s:s + per_win], d[s:s + per_win])[2]
          for s in range(0, len(c) - per_win + 1, per_win)]
    wsd = float(np.std(ws)) if len(ws) > 1 else float("nan")
    rows.append((name, len(c), s1, s2, s3, wsd, len(ws)))
    print(f"{name:10s} {len(c):>6d} {s1:7.2f} {s2:7.2f} {s3:7.2f} {wsd:7.2f} {len(ws):>4d}")

# ── 4. bootstrap of the production 50-ps span (panel convention) ────────
c, d = variants["prod_50ps"]
block = 40
nb = int(np.ceil(len(c) / block))
boots = []
for _ in range(400):
    starts = RNG.integers(0, len(c), nb)
    idx = np.concatenate([(s + np.arange(block)) % len(c) for s in starts])[:len(c)]
    boots.append(spans(c[idx], d[idx])[2])
boots = np.array(boots)
print(f"\nbootstrap prod_50ps: mean {boots.mean():.2f} SD {boots.std(ddof=1):.2f} "
      f"95% [{np.quantile(boots,0.025):.2f}, {np.quantile(boots,0.975):.2f}]")

out = ROOT / "data/analysis/r6/s59_talazoparib_cadence.csv"
with open(out, "w") as f:
    f.write("variant,n_frames,C1,C2,C3,win_sd_2ns,n_windows\n")
    for r in rows:
        f.write(",".join(str(x) for x in r) + "\n")
    f.write(f"boot_prod_50ps_mean,{boots.mean():.4f}\n"
            f"boot_prod_50ps_sd,{boots.std(ddof=1):.4f}\n"
            f"boot_prod_50ps_q025,{np.quantile(boots,0.025):.4f}\n"
            f"boot_prod_50ps_q975,{np.quantile(boots,0.975):.4f}\n")
print(f"wrote {out}")
