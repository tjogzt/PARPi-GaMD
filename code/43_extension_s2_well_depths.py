#!/usr/bin/env python3
"""
43_extension_s2_well_depths.py — extension S2 CV1/CV2 well depths (parameterised weights)

Purpose:   C1-C3 well depths for the three extension-inhibitor S2 systems
           (protein-DNA CV1 and HD-ART CV2) under both reweighting schemes:
             - DBE (textbook dihedral-boost-energy weights, analysis_weights.dat)
             - DFW (force-scaling scheme used for the rebuilt S2 values,
               analysis_weights_dfw.dat; reproduces talazoparib CV2 = 30.6).
           Merged from the former s2_new_drugs_pmf.py and s2_new_drugs_pmf_dfw.py
           (identical logic, weight file as a parameter).
Author:    Tao Zhu (tjogzt@gmail.com)
Created:   2026-09-16 (original scripts); merged 2026-10-05 (repo standardisation)
Inputs:    data/new_drugs_s2/<lig>/analysis_{CV1,CV2}.dat
           data/new_drugs_s2/<lig>/analysis_weights{,_dfw}.dat
Outputs:  console tables (no files written)
Depends:   numpy; common.pmf (run_pyrew)
Run:       python3 code/43_extension_s2_well_depths.py   (from the repo root)
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from common.pmf import run_pyrew  # noqa: E402

DATA = ROOT / "data" / "new_drugs_s2"
LIGANDS = ["fluzoparib", "pamiparib", "senaparib"]
SCHEMES = [("DBE", "analysis_weights.dat"), ("DFW", "analysis_weights_dfw.dat")]

for scheme, wfile in SCHEMES:
    print(f"=== Extension S2 well depths ({scheme} weights) ===")
    print(f'{"drug":<12}{"CV":<6}{"C1":>8}{"C2":>8}{"C3":>8}{"mean+-SD":>12}')
    for d in LIGANDS:
        for cv_name in ["CV1", "CV2"]:
            cv = np.loadtxt(DATA / d / f"analysis_{cv_name}.dat")
            w = np.loadtxt(DATA / d / wfile)
            assert len(cv) == len(w), f"{d} {cv_name}: {len(cv)} vs {len(w)}"
            r = run_pyrew(cv, w, wfmt="%.6f")
            c1, c2, c3 = r.c1, r.c2, r.c3
            mean, sd = np.mean([c1, c2, c3]), np.std([c1, c2, c3])
            print(f"{d:<12}{cv_name:<6}{c1:>8.1f}{c2:>8.1f}{c3:>8.1f}{mean:>8.1f}+-{sd:.1f}")
    print()
