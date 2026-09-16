#!/usr/bin/env python3
"""S2 DBE reconstruction for systems whose gamd.log was not archived (P0-2).

For niraparib/olaparib/rucaparib the per-frame dihedral boost energy dV_D was
not archived. We reconstruct it from first principles using the gamdRunner
lower-dual boost identity (validated frame-exact on APO/talazoparib, where both
the logged dV_D and the logged dihedral energy are available):

    dV_D = 0.5 * (E_dih - E1) * (1 - w_D)

where E_dih is the dihedral-group energy (PeriodicTorsionForce +
CMAPTorsionForce, the gamdRunner "dihedral" boost group per
tools/gamd-openmm/gamd/integrator_factory.py:40) computed with OpenMM from the
archived trajectory, w_D is the archived dihedral force-scaling weight, and E1
is the boost threshold (the prep-phase dihedral energy maximum; calibrated as
production Vmax - offset, with the offset measured on systems with complete
logs).

Pipeline per system:
  1. dihedral_group_energy.py -> E_dih per subsampled frame (stride 10).
  2. Build (beta*dV_D, 0, dV_D) weights aligned with the subsampled CV series.
  3. PyReweighting amdweight_CE -> C1-C3 well depths for CV1 and CV2.

Usage:
  python3 scripts/reconstruct_s2_dbe.py <ligand> <E_csv> <e1_offset>

Validated first on APO/talazoparib (reconstruction vs logged dV_D), then
applied to niraparib/olaparib/rucaparib.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir
from common.pmf import run_pyrew

BETA = 1.0 / (0.001987 * 300.0)
A = analysis_dir()


def main():
    lig = sys.argv[1]
    e_csv = sys.argv[2]
    offset = float(sys.argv[3])  # E1 = production Vmax - offset (kcal/mol)

    d = np.loadtxt(e_csv, delimiter=",", skiprows=1)
    E = d[:, 2]
    # E1 = Vmax - offset (offset calibrated on systems with complete logs)
    E1 = E.max() - offset

    # w_D series: archived DFW weights, subsampled at stride 10 from frame 0
    w_all = np.loadtxt(A / f"sys2_{lig}_CV1_cv_weights.dat")[:, 2]
    w = w_all[::10][: len(E)]
    assert len(w) == len(E), f"{lig}: w {len(w)} vs E {len(E)}"

    dV = 0.5 * (E - E1) * (1 - w)
    dV = np.clip(dV, 0.0, None)
    weights = np.column_stack([BETA * dV, np.zeros_like(dV), dV])
    print(f"{lig}: E1={E1:.1f}  dV_D range {dV.min():.2f}-{dV.max():.2f} "
          f"(mean {dV.mean():.2f}, active {np.mean(dV > 1e-6)*100:.1f}%)")

    out = {}
    for cvn in ["CV1", "CV2"]:
        cv = np.load(A / f"sys2_{lig}_{cvn}_cv.npy")[::10][: len(E)]
        assert len(cv) == len(E)
        r = run_pyrew(cv, weights, wfmt="%.6f")
        vals = [r.c1, r.c2, r.c3]
        out[cvn] = vals
        print(f"{lig:<14}{cvn:<5}{vals[0]:>8.1f}{vals[1]:>8.1f}{vals[2]:>8.1f}"
              f"{np.mean(vals):>8.1f}+-{np.std(vals):.1f}")
    np.save(str(A / f"s2_dbe_{lig}.npy"), out, allow_pickle=True)
    print(f"Saved results/analysis/s2_dbe_{lig}.npy")


if __name__ == "__main__":
    main()
