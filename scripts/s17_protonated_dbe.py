#!/usr/bin/env python3
#
# Purpose:  Protonated-ligand re-run: exact DBE reweighting (S1 + S2). Merged from the former analyze_nira_prot.py / analyze_ruca_prot.py (2026-09-29 code audit; the two differed only in the ligand name).
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   <DATA_ROOT>/rerun_202609/runs_{lig}_prot/{sys1,sys2}/ (archived per-frame CV series)
# Outputs:  Output: results/analysis/{lig}_prot_dbe.csv
# Depends:  argparse, numpy, common
# Run:      python3 scripts/s17_protonated_dbe.py

"""Protonated-ligand re-run: exact DBE reweighting (S1 + S2).
Merged from the former analyze_nira_prot.py / analyze_ruca_prot.py
(2026-09-29 code audit; the two differed only in the ligand name).

Usage: python3 scripts/analyze_prot_dbe.py --lig {nira,ruca}

Data: <DATA_ROOT>/rerun_202609/runs_{lig}_prot/{sys1,sys2}/
Convention: identical to s18_rerun_exact_dbe.py -- dV_D = log col8 (index 7),
    production segment = production-start-step.txt (6.0e6), weight = exp(beta*dV_D).
Output: results/analysis/{lig}_prot_dbe.csv
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root
from common.pmf import run_pyrew
from common.cv_extract import load_dv, s1_cv, s2_cv

BETA = 1.0 / (0.001987 * 300.0)
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis"

LIGANDS = {
    "nira": {"name": "niraparib", "tag": "nira"},
    "ruca": {"name": "rucaparib", "tag": "ruca"},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lig", required=True, choices=["nira", "ruca"])
    args = ap.parse_args()
    lig = LIGANDS[args.lig]
    name = lig["name"]
    tag = lig["tag"]
    base = data_root() / f"rerun_202609/runs_{tag}_prot"

    rows = []
    # ---- S2 ----
    s2dir = base / "sys2"
    if (s2dir / "output.dcd").exists():
        dv, prod = load_dv(s2dir / "gamd.log")
        cv1, cv2 = s2_cv(s2dir / f"{name}_best_prot.prmtop", s2dir / "output.dcd")
        cv1, cv2 = cv1[prod], cv2[prod]
        n = min(len(cv1), len(dv))
        print(f"S2: dV_D mean {dv[:n].mean():.2f} "
              f"CV1 {cv1[:n].mean():.1f} CV2 {cv2[:n].mean():.1f} nprod={n}", flush=True)
        w = np.column_stack([BETA * dv[:n], np.zeros(n), dv[:n]])
        for cvn, cv in [("CV1", cv1), ("CV2", cv2)]:
            r = run_pyrew(cv[:n], w)
            rows.append([f"sys2_{name}_prot", cvn, round(r.c1, 1), round(r.c2, 1),
                         round(r.c3, 1), f"{tag}_prot_dbe"])
            print(f"    S2 {cvn}: C1={r.c1:.1f} C2={r.c2:.1f} C3={r.c3:.1f}", flush=True)
    else:
        print("S2 data missing — skipped")
    # ---- S1 ----
    s1dir = base / "sys1"
    if (s1dir / "output.dcd").exists():
        dv, prod = load_dv(s1dir / "gamd.log")
        cv = s1_cv(s1dir / f"sys1_{name}_prot.prmtop", s1dir / "output.dcd")
        cv = cv[prod]
        n = min(len(cv), len(dv))
        print(f"S1: dV_D mean {dv[:n].mean():.2f} CV mean {cv[:n].mean():.1f} nprod={n}",
              flush=True)
        w = np.column_stack([BETA * dv[:n], np.zeros(n), dv[:n]])
        r = run_pyrew(cv[:n], w)
        rows.append([f"sys1_{name}_prot", "CV", round(r.c1, 1), round(r.c2, 1),
                     round(r.c3, 1), f"{tag}_prot_dbe"])
        print(f"    S1 CV: C1={r.c1:.1f} C2={r.c2:.1f} C3={r.c3:.1f}", flush=True)
    else:
        print("S1 data missing — skipped")
    out = OUT / f"{tag}_prot_dbe.csv"
    with open(out, "w") as f:
        f.write("system,cv,c1,c2,c3,method\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    print(f"WROTE {out}")


if __name__ == "__main__":
    main()
