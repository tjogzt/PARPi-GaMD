#!/usr/bin/env python3
"""s68_trajectory_qc.py — per-frame CV/RMSD/Zn coordination for the S2 panel.

Reads the archived S2 trajectories and writes one tidy CSV:
  system, time_ns, rmsd_ca, zn1_min, zn2_min, zn3_min,
  zn1_ncoord, zn2_ncoord, zn3_ncoord

- CA RMSD is computed against frame 0 of the same trajectory (Kabsch fit on CA).
- Zn–Sγ distances are the minimum over all cysteine SG atoms for each of the
  three structural Zn ions; ncoord counts SG atoms within 2.8 Å.
- Neutral systems use the 26-ns production DCDs (1-ps frames, stride 100);
  protonated systems use the 50-ps DCDs (stride 2), truncated at 26 ns.
No simulation is run.
"""
import argparse
import csv
import warnings
from pathlib import Path

import numpy as np
import MDAnalysis as mda

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/Volumes/tjogzt4T/PARPi_data")
OUT = ROOT / "data/analysis/s68_trajectory_qc.csv"

SYSTEMS = {
    # name: (topology, dcd, frame_dt_ns, stride, max_ns)
    "talazoparib": (DATA / "sys2_talazoparib/sys2_talazoparib.prmtop",
                    DATA / "sys2_talazoparib/output.dcd", 0.001, 100, 26.0),
    "niraparib": (DATA / "sys2_niraparib/sys2_niraparib.prmtop",
                  DATA / "sys2_niraparib/output.dcd", 0.001, 100, 26.0),
    "olaparib": (DATA / "sys2_olaparib/sys2_olaparib.prmtop",
                 DATA / "sys2_olaparib/output.dcd", 0.001, 100, 26.0),
    "rucaparib": (DATA / "sys2_rucaparib/sys2_rucaparib.prmtop",
                  DATA / "sys2_rucaparib/output.dcd", 0.001, 100, 26.0),
    "veliparib": (DATA / "sys2_veliparib/sys2_veliparib.prmtop",
                  DATA / "sys2_veliparib/output.dcd", 26.0 / 1930.0, 7, 26.0),
    "niraparib_prot": (DATA / "rerun_202609/runs_nira_prot/sys2/niraparib_best_prot.prmtop",
                       DATA / "rerun_202609/runs_nira_prot/sys2/output.dcd", 0.05, 2, 26.0),
    "rucaparib_prot": (DATA / "rerun_202609/runs_ruca_prot/sys2/rucaparib_best_prot.prmtop",
                       DATA / "rerun_202609/runs_ruca_prot/sys2/output.dcd", 0.05, 2, 26.0),
    "veliparib_prot": (DATA / "veli_prot/sys2_veliparib_prot/veliparib_best_prot.prmtop",
                       DATA / "veli_prot/sys2_veliparib_prot/output.dcd", 0.05, 2, 26.0),
}


def kabsch_rmsd(P, Q):
    """RMSD of P onto Q after optimal rigid-body superposition (no reflection)."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _s, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    Prot = Pc @ R.T
    return float(np.sqrt(np.mean(np.sum((Prot - Qc) ** 2, axis=1))))


def run_system(name, top, dcd, dt_ns, stride, max_ns):
    if not top.exists() or not dcd.exists():
        raise FileNotFoundError(f"{name}: missing {top} or {dcd}")
    u = mda.Universe(str(top), str(dcd))
    ca = u.select_atoms("protein and name CA")
    zn = u.select_atoms("name ZN")
    sg = u.select_atoms("name SG")
    if len(ca) == 0 or len(zn) != 3 or len(sg) == 0:
        raise ValueError(f"{name}: CA={len(ca)} Zn={len(zn)} SG={len(sg)}")
    n_frames = len(u.trajectory)
    max_frame = min(n_frames - 1, int(round(max_ns / dt_ns)))
    frames = list(range(0, max_frame + 1, stride))
    if frames[-1] != max_frame:
        frames.append(max_frame)
    u.trajectory[0]
    ref = ca.positions.copy()
    rows = []
    for k, f in enumerate(frames):
        u.trajectory[f]
        rmsd = kabsch_rmsd(ca.positions, ref)
        zpos = zn.positions.copy()
        spos = sg.positions.copy()
        d = np.linalg.norm(zpos[:, None, :] - spos[None, :, :], axis=2)
        mins = d.min(axis=1)
        ncoord = (d <= 2.8).sum(axis=1)
        rows.append((name, round(f * dt_ns, 3), round(rmsd, 3),
                     *[round(float(x), 3) for x in mins],
                     *[int(x) for x in ncoord]))
        if (k + 1) % 25 == 0:
            print(f"  {name}: {k+1}/{len(frames)} frames", flush=True)
    print(f"{name}: {len(rows)} frames, {len(ca)} CA, Zn min {mins.min():.2f}-{mins.max():.2f} A, "
          f"min ncoord {ncoord.min()}", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", default=None, help="run one system only")
    args = ap.parse_args()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for name, spec in SYSTEMS.items():
        if args.system and name != args.system:
            continue
        print(f"\n=== {name} ===", flush=True)
        all_rows.extend(run_system(name, *spec))
    if args.system:
        print("single-system test complete; CSV not overwritten")
        return
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["system", "time_ns", "rmsd_ca", "zn1_min", "zn2_min", "zn3_min",
                    "zn1_ncoord", "zn2_ncoord", "zn3_ncoord"])
        w.writerows(all_rows)
    print(f"\nwrote {OUT} ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()
