#!/usr/bin/env python3
"""
s2_new_drugs_pmf_dfw.py — extension S2 CV1/CV2 well depths (DFW weights)

Purpose:   Compute C1-C3 well depths for the extension-inhibitor S2 systems
           using the DFW reweighting scheme (same protocol as the manuscript's
           rebuilt S2 values; reproduces talazoparib CV2 = 30.6).
Inputs:    results/analysis/new_drugs_s2/<drug>/{analysis_CV1,analysis_CV2}.dat,
           analysis_weights_dfw.dat
Outputs:   console table
Depends:   numpy; common (pmf, paths)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir
from common.pmf import run_pyrew

DATA = analysis_dir() / "new_drugs_s2"

drugs = ["fluzoparib", "pamiparib", "senaparib"]
print(f'{"drug":<12}{"CV":<6}{"C1":>8}{"C2":>8}{"C3":>8}{"mean+-SD":>12}')
for d in drugs:
    for cv_name in ["CV1", "CV2"]:
        cv = np.loadtxt(DATA / d / f"analysis_{cv_name}.dat")
        w = np.loadtxt(DATA / d / "analysis_weights_dfw.dat")
        assert len(cv) == len(w), f"{d} {cv_name}: {len(cv)} vs {len(w)}"
        r = run_pyrew(cv, w, wfmt="%.6f")
        c1, c2, c3 = r.c1, r.c2, r.c3
        mean, sd = np.mean([c1, c2, c3]), np.std([c1, c2, c3])
        print(f"{d:<12}{cv_name:<6}{c1:>8.1f}{c2:>8.1f}{c3:>8.1f}{mean:>8.1f}+-{sd:.1f}")
