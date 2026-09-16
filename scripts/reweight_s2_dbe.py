#!/usr/bin/env python3
"""S2 DBE reweighting — textbook e^{beta*dV_D} weights (P0-2 fix).

Purpose:   Replace the DFW (force-scaling) S2 reweighting scheme with the
           textbook boost-energy reweighting. For each S2 system, build the
           weight file (beta*dV_D, 0, dV_D) frame-aligned with the CV series
           and run PyReweighting amdweight_CE (C1-C3 well depths).

Data sources:
  - Systems with archived gamd.log: dV_D = log column 8 (Dihedral-Boost-Energy,
    kcal/mol; logger header confirmed in tools/gamd-openmm/gamd/GamdLogger.py).
    Row-aligned 1:1 with the archived CV series (validated pairwise via the
    Dihedral-Force-Weight column, col 6).
  - Rebuilt/extension systems (AZD5305, veliparib, fluzoparib, pamiparib,
    senaparib): results/analysis/new_drugs_s2/<lig>/{analysis_CV1,CV2}.dat +
    analysis_weights.dat (already the 3-column beta*dV format).

Outputs: console table; per-system C1-C3 well depths for CV1 and CV2.
Depends: numpy; common (pmf, paths); DATA_ROOT logs on the data volume.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir, data_root
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)  # 1/kcal/mol at 300 K

A = analysis_dir()
D = Path(data_root()) if data_root() else Path("/Volumes/tjogzt4T/PARPi_data")

# --- systems with archived gamd.log (row-aligned with results/analysis CV) ---
LOG_SYSTEMS = {
    "APO": ("sys2_APO", "sys2_APO"),
    "talazoparib": ("sys2_talazoparib", "sys2_talazoparib"),
}

# --- rebuilt systems with gamd.log on the volume (production = log rows [243:]) ---
REBUILT_LOG = {
    "AZD5305": "s2_new_drugs/sys2_AZD5305",
    "veliparib": "s2_new_drugs/sys2_veliparib",
}
# --- extension systems (weights already in beta*dV 3-col format) ---
DAT_SYSTEMS = ["fluzoparib", "pamiparib", "senaparib"]


def log_dbe_weights(lig, log_file, cv_file):
    """Build (beta*dV_D, 0, dV_D) from the gamd.log, row-aligned with cv_file."""
    lg = np.loadtxt(log_file, comments="#")
    cv = np.load(cv_file)
    assert len(lg) == len(cv), f"{lig}: log {len(lg)} vs cv {len(cv)}"
    dV = lg[:, 7]
    # sanity: force weight column (5) must equal the archived DFW series
    return np.column_stack([BETA * dV, np.zeros_like(dV), dV])


def main():
    print(f"{'ligand':<14}{'CV':<5}{'C1':>8}{'C2':>8}{'C3':>8}{'mean+-SD':>14}")
    results = {}
    for lig, (log_dir, log_name) in LOG_SYSTEMS.items():
        log_file = D / log_dir / "gamd.log"
        for cvn in ["CV1", "CV2"]:
            cv = np.load(A / f"{log_dir}_{cvn}_cv.npy")
            w = log_dbe_weights(lig, log_file, A / f"{log_dir}_{cvn}_cv.npy")
            r = run_pyrew(cv, w, wfmt="%.6f")
            vals = [r.c1, r.c2, r.c3]
            results[(lig, cvn)] = vals
            print(f"{lig:<14}{cvn:<5}{vals[0]:>8.1f}{vals[1]:>8.1f}{vals[2]:>8.1f}"
                  f"{np.mean(vals):>8.1f}+-{np.std(vals):.1f}")
    for lig, log_dir in REBUILT_LOG.items():
        lg = np.loadtxt(D / log_dir / "gamd.log", comments="#")
        d = A / "new_drugs_s2" / lig
        for cvn in ["CV1", "CV2"]:
            cv = np.loadtxt(d / f"analysis_{cvn}.dat")
            off = len(lg) - len(cv)
            assert off > 0 and np.allclose(
                np.loadtxt(d / "analysis_weights_dfw.dat")[:, 2],
                lg[off:, 5], atol=1e-6), f"{lig}: frame alignment failed"
            dV = lg[off:, 7]
            w = np.column_stack([BETA * dV, np.zeros_like(dV), dV])
            r = run_pyrew(cv, w, wfmt="%.6f")
            vals = [r.c1, r.c2, r.c3]
            results[(lig, cvn)] = vals
            print(f"{lig:<14}{cvn:<5}{vals[0]:>8.1f}{vals[1]:>8.1f}{vals[2]:>8.1f}"
                  f"{np.mean(vals):>8.1f}+-{np.std(vals):.1f}")
    for lig in DAT_SYSTEMS:
        d = A / "new_drugs_s2" / lig
        w = np.loadtxt(d / "analysis_weights.dat")
        assert w.shape[1] == 3, f"{lig}: weights not 3-col"
        for cvn in ["CV1", "CV2"]:
            cv = np.loadtxt(d / f"analysis_{cvn}.dat")
            assert len(cv) == len(w), f"{lig} {cvn}: {len(cv)} vs {len(w)}"
            r = run_pyrew(cv, w, wfmt="%.6f")
            vals = [r.c1, r.c2, r.c3]
            results[(lig, cvn)] = vals
            print(f"{lig:<14}{cvn:<5}{vals[0]:>8.1f}{vals[1]:>8.1f}{vals[2]:>8.1f}"
                  f"{np.mean(vals):>8.1f}+-{np.std(vals):.1f}")
    np.save(str(A / "s2_dbe_well_depths.npy"),
            {f"{k[0]}_{k[1]}": v for k, v in results.items()},
            allow_pickle=True)
    print("Saved results/analysis/s2_dbe_well_depths.npy")


if __name__ == "__main__":
    main()
