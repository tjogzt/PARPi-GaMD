#!/usr/bin/env python3
#
# Purpose:  Collect C1/C2/C3 well depths + C2-C3 PMF Pearson r for the 7 S1 systems (DBE weights, 24 ns) — feeds the revised Table S1.
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   (see docstring)
# Depends:  numpy, common, MDAnalysis
# Run:      python3 scripts/s22_collect_s1_table1.py

"""Collect C1/C2/C3 well depths + C2-C3 PMF Pearson r for the 7 S1 systems
(DBE weights, 24 ns) — feeds the revised Table S1."""
import sys
import subprocess
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root, pyrew as pyrew_path

BETA = 1.0 / (0.001987 * 300.0)
DATA = data_root()
REPL = Path(__file__).resolve().parents[1] / "results" / "analysis" / "replicates"
RERUN = data_root() / "rerun_202609"
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis"
import tempfile
TMP = Path(tempfile.mkdtemp(prefix="s1t1_collect_"))


def load_pmf3(tmp_dir):
    lines = [l for l in (tmp_dir / "pmf-c3-cv.dat.xvg").read_text().splitlines()
             if not l.startswith(("#", "@")) and l.strip()]
    return np.array([[float(x) for x in l.split()] for l in lines])


def run_and_r(cv, w, save_as=None):
    TMP.mkdir(parents=True, exist_ok=True)
    np.savetxt(TMP / "cv.dat", cv, fmt="%.4f")
    np.savetxt(TMP / "weights.dat", w, fmt="%.6f")
    subprocess.run(["python3", str(pyrew_path()), "-input", "cv.dat", "-T", "300",
                    "-disc", "0.1", "-Emax", "20", "-cutoff", "2",
                    "-job", "amdweight_CE", "-weight", "weights.dat"],
                   cwd=TMP, capture_output=True)
    def load_pmf(name):
        lines = [l for l in (TMP / name).read_text().splitlines()
                 if not l.startswith(("#", "@")) and l.strip()]
        return np.array([[float(x) for x in l.split()] for l in lines])

    pmf2 = load_pmf("pmf-c2-cv.dat.xvg")
    pmf3 = load_pmf("pmf-c3-cv.dat.xvg")
    if save_as is not None:
        np.savetxt(save_as, pmf3, fmt="%.6f %.4f")
    # interpolate c3 onto c2 grid for r
    n = min(len(pmf2), len(pmf3))
    x2, y2 = pmf2[:n, 0], pmf2[:n, 1]
    x3, y3 = pmf3[:n, 0], pmf3[:n, 1]
    y3i = np.interp(x2, x3, y3)
    r = np.corrcoef(y2, y3i)[0, 1]
    d2 = y2.max() - y2.min()
    d3 = y3.max() - y3.min()
    d1 = float(load_pmf("pmf-c1-cv.dat.xvg")[:, 1].max()
               - load_pmf("pmf-c1-cv.dat.xvg")[:, 1].min())
    return d1, d2, d3, r


def w_dbe(dv):
    return np.column_stack([BETA * dv, np.zeros(len(dv)), dv])


def main():
    TMP.mkdir(parents=True, exist_ok=True)
    import MDAnalysis as mda
    rows = {}
    for tag, prm in [("APO", "sys1_APO.prmtop"), ("AZD5305", "system.prmtop")]:
        d = DATA / f"sys1_{tag}"
        u = mda.Universe(str(d / prm), str(d / "output.dcd"))
        hd = u.select_atoms("resid 1:119 and name CA")
        art = u.select_atoms("resid 120:345 and name CA")
        cv = np.array([np.linalg.norm(hd.center_of_mass() - art.center_of_mass())
                       for ts in u.trajectory])
        log = np.loadtxt(d / "gamd.log", comments="#")
        dv = log[:, 7]
        n = min(len(cv), len(dv))
        d1, d2, d3, r = run_and_r(cv[:n], w_dbe(dv[:n]),
                              save_as=OUT / f"sys1_{tag}_pmf_c3.xvg")
        rows[tag] = (d1, d2, d3, r)
        print(f"{tag}: C1={d1:.1f} C2={d2:.1f} C3={d3:.1f} r={r:.2f}", flush=True)

    for tag in ["talazoparib", "veliparib"]:
        rs = []
        pmf3_seeds = []
        for s in ["s7", "s49", "s123"]:
            d = DATA / "review_runs" / f"sys1_{tag}_{s}"
            log = np.loadtxt(d / "gamd.log", comments="#")
            prod = log[:, 1] >= 6.0e6
            dv = log[:, 7][prod][:480]
            cv_full = np.loadtxt(REPL / f"sys1_{tag}_{s}" / "cv.dat")
            cv = cv_full[prod][:480]
            n = min(len(cv), len(dv))
            d1, d2, d3, r = run_and_r(cv[:n], w_dbe(dv[:n]))
            pmf3_seeds.append(load_pmf3(TMP))
            rs.append((d1, d2, d3, r))
            print(f"{tag}_{s}: C1={d1:.1f} C2={d2:.1f} C3={d3:.1f} r={r:.2f}", flush=True)
        # save the 3-seed mean C3 curve (mean-curve convention), not a single seed
        grid = np.linspace(min(p[:, 0].min() for p in pmf3_seeds),
                           max(p[:, 0].max() for p in pmf3_seeds), 200)
        interp = [np.interp(grid, p[:, 0], p[:, 1]) for p in pmf3_seeds]
        mean_curve = np.column_stack([grid, np.mean(interp, axis=0)])
        np.savetxt(OUT / f"sys1_{tag}_pmf_c3.xvg", mean_curve, fmt="%.6f %.4f")
        rows[tag] = tuple(np.mean(np.array(rs), axis=0))
        print(f"{tag} mean: {['%.1f' % x for x in rows[tag][:3]]} r={rows[tag][3]:.2f}")

    for lig in ["niraparib", "olaparib", "rucaparib"]:
        tag = f"sys1_{lig}"
        log = np.loadtxt(RERUN / "rerun_logs" / f"{tag}_gamd.log", comments="#")
        prod = log[:, 1] >= 6.0e6
        dv = log[:, 7][prod]
        u = mda.Universe(str(RERUN / "sys_files" / f"{tag}.prmtop"),
                         str(RERUN / "runs" / tag / "gamd_out" / "output.dcd"))
        hd = u.select_atoms("resid 1:119 and name CA")
        art = u.select_atoms("resid 120:345 and name CA")
        cv_all = np.array([np.linalg.norm(hd.center_of_mass() - art.center_of_mass())
                           for ts in u.trajectory])
        cv = cv_all[prod]
        n = min(len(cv), len(dv))
        d1, d2, d3, r = run_and_r(cv[:n], w_dbe(dv[:n]),
                              save_as=OUT / f"sys1_{lig}_pmf_c3.xvg")
        rows[lig] = (d1, d2, d3, r)
        print(f"{lig}: C1={d1:.1f} C2={d2:.1f} C3={d3:.1f} r={r:.2f}", flush=True)

    order = ["APO", "AZD5305", "niraparib", "olaparib", "rucaparib",
             "talazoparib", "veliparib"]
    with open(OUT / "s1_table1_rows.csv", "w") as f:
        f.write("ligand,C1,C2,C3,r\n")
        for k in order:
            d1, d2, d3, r = rows[k]
            f.write(f"{k},{d1:.2f},{d2:.2f},{d3:.2f},{r:.2f}\n")
    print("WROTE s1_table1_rows.csv")


if __name__ == "__main__":
    main()
