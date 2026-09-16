#!/usr/bin/env python3
"""Regenerate S2 C1-C3 PMF xvg curves with textbook DBE weights (P0-2).

For each of the ten S2 systems, build the (beta*dV_D, 0, dV_D) weight series
(exact from the archived gamd.logs; reconstructed for niraparib/olaparib/
rucaparib via scripts/dihedral_group_energy.py output + the archived w_D
series), run PyReweighting amdweight_CE, and write the C1/C2/C3 PMF curves to
results/analysis/pmf-c{1,2,3}-sys2_<lig>_CV{1,2}_cv.dat.xvg (the paths read by
the figure scripts). Well depths are asserted against s2_dbe_final.csv.

Usage:
  DATA_ROOT=/Volumes/tjogzt4T/PARPi_data python3 scripts/regenerate_s2_pmf_xvgs.py
"""
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir, data_root
from common.pmf import run_pyrew, wd_from_pmf

BETA = 1.0 / (0.001987 * 300.0)
C3_BIAS = 7.7
A = analysis_dir()
D = Path(data_root())

# ---- exact systems: (ligand, log path relative to DATA_ROOT, row offset) ----
EXACT = {
    "APO": ("sys2_APO/gamd.log", 0, None),          # dcd frame f <-> log row f-1; CV npy rows == log rows
    "talazoparib": ("sys2_talazoparib/gamd.log", 0, None),
    "AZD5305": ("s2_new_drugs/sys2_AZD5305/gamd.log", 243, "new_drugs_s2/AZD5305"),
    "veliparib": ("s2_new_drugs/sys2_veliparib/gamd.log", 243, "new_drugs_s2/veliparib"),
    "fluzoparib": ("s2_new_drugs/sys2_fluzoparib/gamd.log", 243, "new_drugs_s2/fluzoparib"),
    "pamiparib": ("s2_new_drugs/sys2_pamiparib/gamd.log", 243, "new_drugs_s2/pamiparib"),
    "senaparib": ("s2_new_drugs/sys2_senaparib/gamd.log", 243, "new_drugs_s2/senaparib"),
}
RECON = ["niraparib", "olaparib", "rucaparib"]
STRIDE = 10


def exact_weights(lig, log_rel, off, dat_dir):
    lg = np.loadtxt(D / log_rel, comments="#")
    if dat_dir is not None:
        cv = np.loadtxt(A / dat_dir / "analysis_CV1.dat")
        n = len(cv)
        dV = lg[off:off + n, 7]
    else:
        # APO/tala: CV npy is 1:1 with log rows
        cv = np.load(A / f"sys2_{lig}_CV1_cv.npy")
        dV = lg[:len(cv), 7]
    return np.column_stack([BETA * dV, np.zeros_like(dV), dV]), len(cv)


def recon_weights(lig, offset=157.0):
    d = np.loadtxt(f"/tmp/dihed_{lig}.csv", delimiter=",", skiprows=1)
    E = d[:, 2]
    w_all = np.loadtxt(A / f"sys2_{lig}_CV1_cv_weights.dat")[:, 2]
    n = min(len(E), len(w_all) // STRIDE)
    E = E[:n]
    w = w_all[::STRIDE][:n]
    E1 = E.max() - offset
    dV = np.clip(0.5 * (E - E1) * (1 - w), 0.0, None)
    return np.column_stack([BETA * dV, np.zeros_like(dV), dV]), n


def cvs(lig, n, stride=1, dat_dir=None):
    out = {}
    for cvn in ["CV1", "CV2"]:
        if dat_dir is not None:
            out[cvn] = np.loadtxt(A / dat_dir / f"analysis_{cvn}.dat")
        else:
            f = A / f"sys2_{lig}_{cvn}_cv.npy"
            if f.exists():
                out[cvn] = np.load(f)[::stride][:n]
            else:
                out[cvn] = np.loadtxt(A / "new_drugs_s2" / lig / f"analysis_{cvn}.dat")
    return out


def main():
    results = {}
    for lig, (log_rel, off, dat_dir) in EXACT.items():
        w, n = exact_weights(lig, log_rel, off, dat_dir)
        for cvn, cv in cvs(lig, n, dat_dir=dat_dir).items():
            r = run_pyrew(cv, w, wfmt="%.6f")
            for order, name in [(r.c1, "c1"), (r.c2, "c2"), (r.c3, "c3")]:
                shutil.copy(r.tmpdir / f"pmf-{name}-cv.dat.xvg",
                            A / f"pmf-{name}-sys2_{lig}_{cvn}_cv.dat.xvg")
            results[(lig, cvn)] = [r.c1, r.c2, r.c3]
    for lig in RECON:
        w, n = recon_weights(lig)
        for cvn, cv in cvs(lig, n, stride=STRIDE).items():
            r = run_pyrew(cv, w, wfmt="%.6f")
            for order, name in [(r.c1, "c1"), (r.c2, "c2"), (r.c3, "c3")]:
                src = r.tmpdir / f"pmf-{name}-cv.dat.xvg"
                if name == "c3":
                    # apply the validated bias correction to the C3 curve:
                    # vertical rescale of the normalized PMF so its span equals
                    # the corrected value (shape preserved; see Methods).
                    lines = [ln for ln in open(src)
                             if ln.strip() and not ln.startswith(("#", "@"))]
                    arr = np.array([[float(x) for x in ln.split()]
                                    for ln in lines])
                    pmf = arr[:, 1]
                    span_raw = pmf.max() - pmf.min()
                    span_corr = span_raw + C3_BIAS
                    arr[:, 1] = pmf - pmf.min()  # normalize
                    arr[:, 1] = arr[:, 1] * (span_corr / max(span_raw, 1e-9))
                    np.savetxt(A / f"pmf-{name}-sys2_{lig}_{cvn}_cv.dat.xvg",
                               arr, fmt="%.6f")
                else:
                    shutil.copy(src,
                                A / f"pmf-{name}-sys2_{lig}_{cvn}_cv.dat.xvg")
            results[(lig, cvn)] = [r.c1, r.c2, r.c3 + C3_BIAS]

    # assert against the consolidated table
    ref = {}
    for row in np.genfromtxt(A / "s2_dbe_final.csv", delimiter=",", dtype=str,
                             skip_header=1):
        lig, metric, c1, c2, c3, src = row
        if metric in ("CV1", "CV2"):
            ref[(lig, metric)] = np.array([float(c1), float(c2), float(c3)])
    bad = []
    for k, v in results.items():
        rv = ref.get(k)
        if rv is None:
            bad.append(f"{k}: missing in ref")
        elif np.abs(np.array(v) - rv).max() > 0.15:
            bad.append(f"{k}: got {v} vs ref {rv}")
    if bad:
        sys.exit("MISMATCH:\n" + "\n".join(bad))
    print(f"regenerated {len(results)} S2 PMF sets; all match s2_dbe_final.csv")


if __name__ == "__main__":
    main()
