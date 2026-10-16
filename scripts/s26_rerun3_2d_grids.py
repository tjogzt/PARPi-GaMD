#!/usr/bin/env python3
#
# Purpose:  2D C3 PMF grids for the three re-simulated S2 systems (niraparib, olaparib, rucaparib) from the 2026-09 rerun trajectories, replacing the grids previously computed from reconstructed boost energies. Protocol identical to
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   (see docstring)
# Depends:  tempfile, numpy, common, MDAnalysis
# Run:      python3 scripts/s26_rerun3_2d_grids.py

"""2D C3 PMF grids for the three re-simulated S2 systems (niraparib, olaparib,
rucaparib) from the 2026-09 rerun trajectories, replacing the grids previously
computed from reconstructed boost energies. Protocol identical to
scripts/compute_2d_dbe.py (PyReweighting-2D, amdweight_CE, T=300, disc 0.5/0.5,
Emax 20, cutoff 4)."""
import subprocess
import sys
import tempfile
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir, data_root, pyrew

BETA = 1.0 / (0.001987 * 300.0)
A = Path(analysis_dir())
D = Path(data_root())
RERUN = D / 'rerun_202609'
PY2D = Path(pyrew()).parent / "PyReweighting-2D.py"
PROD_STEP = 6.0e6


def rerun_cvs(lig):
    import MDAnalysis as mda
    prmtop = D / f"sys2_{lig}" / f"sys2_{lig}.prmtop"
    u = mda.Universe(str(prmtop),
                     str(RERUN / "runs" / f"sys2_{lig}" / "gamd_out" / "output.dcd"))
    hd = u.select_atoms("resid 662:787 and name CA")
    art = u.select_atoms("resid 788:1011 and name CA")
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    cv1, cv2 = [], []
    for ts in u.trajectory:
        cv1.append(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
        cv2.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    return np.array(cv1), np.array(cv2)


def run_2d(cv1, cv2, dV):
    tmp = Path(tempfile.mkdtemp())
    np.savetxt(tmp / "cv2d.dat", np.column_stack([cv1, cv2]), fmt="%.4f")
    np.savetxt(tmp / "w.dat",
               np.column_stack([BETA * dV, np.zeros_like(dV), dV]),
               fmt="%.6f")
    r = subprocess.run(
        ["python3", str(PY2D), "-input", "cv2d.dat", "-T", "300",
         "-discX", "0.5", "-discY", "0.5", "-Emax", "20", "-cutoff", "4",
         "-job", "amdweight_CE", "-weight", "w.dat"],
        cwd=tmp, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stdout[-500:] + r.stderr[-500:])
    out = tmp / "pmf-c3-cv2d.dat.xvg"
    lines = [ln for ln in open(out) if ln.strip() and not ln.startswith(("#", "@"))]
    arr = np.array([[float(x) for x in ln.split()] for ln in lines])
    xs = np.unique(arr[:, 0])
    ys = np.unique(arr[:, 1])
    F = arr[:, 2].reshape(len(ys), len(xs))
    return xs, ys, F


def main():
    for lig in ["niraparib", "olaparib", "rucaparib"]:
        log = np.loadtxt(RERUN / "rerun_logs" / f"sys2_{lig}_gamd.log",
                         comments="#")
        prod = log[:, 1] >= PROD_STEP
        dV = log[:, 7][prod]
        cv1, cv2 = rerun_cvs(lig)
        n = min(len(cv1[prod]), len(dV))
        xs, ys, F = run_2d(cv1[prod][:n], cv2[prod][:n], dV[:n])
        np.savez(A / f"2d_c3_{lig}.npz", X=xs, Y=ys, F=F)
        print(f"{lig}: grid {F.shape}, F range {F.min():.1f}-{F.max():.1f} "
              f"kcal/mol, CV1 {cv1[prod][:n].min():.1f}-{cv1[prod][:n].max():.1f}")
    print("DONE")


if __name__ == "__main__":
    main()
