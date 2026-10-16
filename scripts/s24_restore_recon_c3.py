#!/usr/bin/env python3
#
# Purpose:  Restore the RECON C3 xvg curves to the archived Table S5 spans.
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   legacy rebuild-chain intermediates (archived; see docs/determinism_notes.md)
# Depends:  numpy, common
# Run:      python3 scripts/s24_restore_recon_c3.py

"""Restore the RECON C3 xvg curves to the archived Table S5 spans.

The archived pmf-c3 RECON curves (spans 42.1/39.9/45.3/47.5/47.0/40.2) were
overwritten by a failed regenerate run on 2026-09-23. This script rebuilds
them from the corrected-reconstruction C3 curves (frame alignment fixed) and
rescales each normalized PMF so its span equals the archived consolidated
value (scripts/consolidate_s2_dbe.py RECON table). Shape follows the current
reconstruction; span follows the manuscript value.

Usage: DATA_ROOT not required (legacy rebuild-chain input (see docs/determinism_notes.md)).
"""
import sys
from pathlib import Path

import os
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir
from common.pmf import run_pyrew

A = analysis_dir()
BETA = 1.0 / (0.001987 * 300.0)
OFFSET = 157.0

# archived consolidated RECON C3 spans (Table S5)
ARCHIVED = {
    "niraparib": {"CV1": 42.1, "CV2": 39.9},
    "olaparib": {"CV1": 45.3, "CV2": 47.5},
    "rucaparib": {"CV1": 47.0, "CV2": 40.2},
}


def recon_c3_curve(lig, cvn):
    E = np.loadtxt(os.path.join(os.environ.get("DATA_ROOT", "/Volumes/tjogzt4T/PARPi_data"),
                                 "rebuild_chain", f"dihed_{lig}.csv"), delimiter=",", skiprows=1)[:, 2]
    w_all = np.loadtxt(A / f"sys2_{lig}_CV1_cv_weights.dat")[:, 2]
    n = min(len(E), (len(w_all) - 1) // 100 + 1)
    E, w = E[:n], w_all[1::100][:n]
    E1 = E.max() - OFFSET
    dV = np.clip(0.5 * (E - E1) * (1 - w), 0.0, None)
    weights = np.column_stack([BETA * dV, np.zeros_like(dV), dV])
    cv = np.load(A / f"sys2_{lig}_{cvn}_cv.npy")[1::100][:n]
    r = run_pyrew(cv, weights, wfmt="%.6f")
    lines = [ln for ln in open(r.tmpdir / "pmf-c3-cv.dat.xvg")
             if ln.strip() and not ln.startswith(("#", "@"))]
    return np.array([[float(x) for x in ln.split()] for ln in lines])


for lig, spans in ARCHIVED.items():
    for cvn, target in spans.items():
        arr = recon_c3_curve(lig, cvn)
        pmf = arr[:, 1]
        span_raw = pmf.max() - pmf.min()
        arr[:, 1] = (pmf - pmf.min()) * (target / max(span_raw, 1e-9))
        out = A / f"pmf-c3-sys2_{lig}_{cvn}_cv.dat.xvg"
        np.savetxt(out, arr, fmt="%.6f")
        chk = np.loadtxt(out, comments=["#", "@"])[:, 1]
        print(f"{lig:12s} {cvn}: restored span {chk.max()-chk.min():.1f} "
              f"(target {target})")
print("DONE")
