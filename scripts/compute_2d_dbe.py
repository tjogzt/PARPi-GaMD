#!/usr/bin/env python3
"""2D cumulant-expansion PMFs for the S2 landscape figure (P0-2 DBE).

Computes the 2D free-energy surfaces (PyReweighting-2D, amdweight_CE,
T = 300 K) for talazoparib, veliparib (rebuilt), and APO from the CV series
and the textbook DBE weights, and writes the C3 PMF grids as npy files for
code/16-2d_landscape.R.

Outputs: results/analysis/2d_c3_<lig>.npz (arrays: X, Y, F)
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir, data_root, pyrew

BETA = 1.0 / (0.001987 * 300.0)
A = analysis_dir()
D = Path(data_root())
PY2D = Path(pyrew()).parent / "PyReweighting-2D.py"


def dbe_weights(lig):
    """(CV1, CV2, dV) per system, matching regenerate_s2_pmf_xvgs.py."""
    if lig in ("AZD5305", "veliparib", "fluzoparib", "pamiparib", "senaparib"):
        log = D / "s2_new_drugs" / f"sys2_{lig}" / "gamd.log"
        lg = np.loadtxt(log, comments="#")
        cv1 = np.loadtxt(A / "new_drugs_s2" / lig / "analysis_CV1.dat")
        cv2 = np.loadtxt(A / "new_drugs_s2" / lig / "analysis_CV2.dat")
        off = len(lg) - len(cv1)
        dV = lg[off:, 7]
    elif lig in ("niraparib", "olaparib", "rucaparib"):
        d = np.loadtxt(f"/tmp/dihed_{lig}.csv", delimiter=",", skiprows=1)
        E = d[:, 2]
        w_all = np.loadtxt(A / f"sys2_{lig}_CV1_cv_weights.dat")[:, 2]
        n = min(len(E), len(w_all) // 10)
        w = w_all[::10][:n]
        E1 = E[:n].max() - 157.0
        dV = np.clip(0.5 * (E[:n] - E1) * (1 - w), 0.0, None)
        cv1 = np.load(A / f"sys2_{lig}_CV1_cv.npy")[::10][:n]
        cv2 = np.load(A / f"sys2_{lig}_CV2_cv.npy")[::10][:n]
    else:
        log = D / f"sys2_{lig}" / "gamd.log"
        lg = np.loadtxt(log, comments="#")
        cv1 = np.load(A / f"sys2_{lig}_CV1_cv.npy")
        cv2 = np.load(A / f"sys2_{lig}_CV2_cv.npy")
        dV = lg[: len(cv1), 7]
    assert len(cv1) == len(dV)
    return cv1, cv2, dV


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
    for lig in ["talazoparib", "veliparib", "APO", "AZD5305",
                "niraparib", "olaparib", "rucaparib"]:
        cv1, cv2, dV = dbe_weights(lig)
        xs, ys, F = run_2d(cv1, cv2, dV)
        np.savez(A / f"2d_c3_{lig}.npz", X=xs, Y=ys, F=F)
        print(f"{lig}: grid {F.shape}, F range {F.min():.1f}-{F.max():.1f} "
              f"kcal/mol, CV1 {cv1.min():.1f}-{cv1.max():.1f}, "
              f"CV2 {cv2.min():.1f}-{cv2.max():.1f}")


if __name__ == "__main__":
    main()
