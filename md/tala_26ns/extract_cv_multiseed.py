#!/usr/bin/env python3
"""Extract per-frame CVs + GaMD boost energy from a completed run, on the GPU host.

Why: pulling 2.5 GB DCDs per seed just to compute spans is slow and unnecessary.
Run this on the instance, download the ~50 KB CSV, and do spans/associations locally.

Usage (PARPi S2 defaults):
  python extract_cv_multiseed.py PRMTOP DCD GAMDLOG OUT.csv
  # S1-style construct: override the selections
  python extract_cv_multiseed.py x.prmtop out.dcd gamd.log cv.csv \
      --sel1 'protein and resid 1:119 and name CA' \
      --sel2 'protein and resid 120:345 and name CA' --sel3 '' --sel4 ''

Frame->step mapping is step = interval_steps * k (k = 1-based DCD frame), which is
how gamdRunner/DCD frames align; do NOT pair by gamd.log line number (the log may
contain headers and buffered gaps). Missing log steps become NaN and are counted.
"""
import argparse

import numpy as np
import MDAnalysis as mda


def read_gamd_log(path, dv_col=7):
    steps, dvs = [], []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) <= dv_col:
                continue
            try:
                step = int(float(parts[1]))
                dv = float(parts[dv_col])
            except ValueError:
                continue
            steps.append(step)
            dvs.append(dv)
    return dict(zip(steps, dvs))


def com_distance(atoms_a, atoms_b):
    return float(np.linalg.norm(atoms_a.center_of_mass() - atoms_b.center_of_mass()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prmtop")
    ap.add_argument("dcd")
    ap.add_argument("gamdlog")
    ap.add_argument("out")
    ap.add_argument("--interval-steps", type=int, default=25000,
                    help="DCD reporting interval in MD steps (default 25000 = 50 ps at 2 fs)")
    ap.add_argument("--dv-col", type=int, default=7,
                    help="0-based gamd.log column of dV_D (default 7 = 8th column)")
    ap.add_argument("--sel1", default="protein and name CA")
    ap.add_argument("--sel2", default="(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 "
                                      "or resname DA3 DT3 DG3 DC3) and not element H")
    ap.add_argument("--sel3", default="resid 662:787 and name CA")
    ap.add_argument("--sel4", default="resid 788:1011 and name CA")
    args = ap.parse_args()

    dv_by_step = read_gamd_log(args.gamdlog, args.dv_col)
    u = mda.Universe(args.prmtop, args.dcd)
    a1, a2 = u.select_atoms(args.sel1), u.select_atoms(args.sel2)
    a3 = u.select_atoms(args.sel3) if args.sel3 else None
    a4 = u.select_atoms(args.sel4) if args.sel4 else None
    if len(a1) == 0 or len(a2) == 0:
        raise SystemExit("empty selection: %r / %r" % (args.sel1, args.sel2))
    dt = float(u.trajectory.dt) or args.interval_steps * 0.002

    missing = 0
    n = 0
    with open(args.out, "w") as fh:
        fh.write("frame,time_ps,step,dV_D,cv1,cv2\n")
        for k, _ts in enumerate(u.trajectory, start=1):
            step = args.interval_steps * k
            dv = dv_by_step.get(step)
            if dv is None:
                missing += 1
                dv = np.nan
            cv1 = com_distance(a1, a2)
            cv2 = com_distance(a3, a4) if a3 is not None and a4 is not None else np.nan
            fh.write(f"{k},{k * dt:.3f},{step},{dv:.6f},{cv1:.6f},{cv2:.6f}\n")
            n += 1

    print(f"wrote {args.out}: {n} frames, dt={dt} ps, missing dV_D={missing}")


if __name__ == "__main__":
    main()
