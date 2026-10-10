#!/usr/bin/env python3
"""r6f convergence recompute — PART A: extension S1 cumulative convergence (audited convention).

For each extension inhibitor (fluzoparib, pamiparib, senaparib):
  - full 200-ns occupied-support C3 span (cross-check vs r6c master values)
  - 22-ns truncation span (matched-length sensitivity)
  - cumulative curve at t = 20..200 ns (20-ns steps), audited c3_v2 convention
Outputs: results/analysis/extension_s1_convergence.csv (new) + r6f_extension_summary.txt
"""
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path("/Users/taozhu/clacky_workspace/PARPi_design")
DATA = Path("/Volumes/tjogzt4T/PARPi_data")
AUD = Path(os.path.expanduser("~/.hermes/cache/scratch/audit"))
PYR = AUD / "PyRew-1D-dump.py"
OUTDIR = REPO / "results/analysis"
R6 = OUTDIR / "r6"


def run(cv, w, tag):
    tmp = Path(tempfile.mkdtemp(prefix=f"r6f_{tag}_"))
    np.savetxt(tmp / "cv.dat", cv, fmt="%.4f")
    np.savetxt(tmp / "weights.dat", w, fmt="%.4f")
    dump = AUD / f"r6f_{tag}.npz"
    env = dict(os.environ, PYREW_DUMP=str(dump))
    r = subprocess.run(["/opt/anaconda3/bin/python3", str(PYR), "-input", "cv.dat", "-T", "300",
                        "-disc", "0.1", "-Emax", "20", "-cutoff", "2", "-job", "amdweight_CE",
                        "-weight", "weights.dat"],
                       cwd=tmp, capture_output=True, text=True, env=env)
    assert (tmp / "pmf-c3-cv.dat.xvg").exists(), f"{tag}: {r.stdout[-200:]}|{r.stderr[-200:]}"
    z = np.load(dump)
    h = z["hist"]
    m = h >= 2
    return {c: (float(z[c][m].max() - z[c][m].min()) if m.sum() else np.nan)
            for c in ("pmf_c1", "pmf_c2", "pmf_c3")}


MASTER_FULL = {"fluzoparib": 22.8, "pamiparib": 16.3, "senaparib": 12.5}  # C3 occupied-support, full 200 ns
rows, summary = [], []
for lig in ["fluzoparib", "pamiparib", "senaparib"]:
    d = DATA / "s1_new_drugs" / f"s1_{lig}"
    cv = np.loadtxt(d / "analysis_cv.dat")
    w = np.loadtxt(d / "analysis_weights.dat")
    if w.ndim == 1:
        w = np.column_stack([w, np.zeros(len(w)), w])
    # cumulative curve at 20..200 ns (50 ps/frame)
    for t in range(20, 201, 20):
        n = t * 20
        r = run(cv[:n], w[:n], f"{lig}_t{t}")
        rows.append({"ligand": lig, "time_ns": t, "well_depth": round(r["pmf_c3"], 6)})
        print(f"  {lig} t={t:3d} ns  C3v={r['pmf_c3']:.2f}", flush=True)
    full = run(cv, w, f"{lig}_full")
    trunc = run(cv[:440], w[:440], f"{lig}_trunc22")  # 22 ns = 440 frames
    ok = abs(full["pmf_c3"] - MASTER_FULL[lig]) < 0.15
    summary.append(f"{lig}: full C3v={full['pmf_c3']:.2f} (master {MASTER_FULL[lig]}, {'OK' if ok else 'CHECK'})"
                   f" | trunc22 C3v={trunc['pmf_c3']:.2f}"
                   f" | C1v/C2v full = {full['pmf_c1']:.2f}/{full['pmf_c2']:.2f}")
    print(summary[-1], flush=True)

# backup + write
old = OUTDIR / "extension_s1_convergence.csv"
if old.exists():
    import shutil
    shutil.copy2(old, str(old) + ".pre_r6f")
with open(old, "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=["ligand", "time_ns", "well_depth"])
    wr.writeheader()
    wr.writerows(rows)
R6.mkdir(parents=True, exist_ok=True)
(R6 / "r6f_extension_summary.txt").write_text("\n".join(summary) + "\n")
print("WROTE extension_s1_convergence.csv + r6f_extension_summary.txt")
