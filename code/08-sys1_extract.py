#!/usr/bin/env python3
"""
08-sys1_extract.py — Robust CV + dV extraction for sys1 PMF.
Fixes: proper gamd.log column detection with leading tab handling.
Output: CV file (1 col) + weight file (3 col: dV*beta, 0, dV) for PyReweighting.
"""
from pathlib import Path
import numpy as np
import MDAnalysis as mda
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root

kT = 0.596  # kcal/mol at 300K, kT = 0.001987*300
STRIDE = 100

def log(msg):
    print(f"[cv] {msg}", flush=True)

def parse_gamd_log(path):
    """Parse gamd.log and extract Total-Boost-Energy-Potential for all frames.
    Returns: boost_arr (1D np.array), boost_idx (column index in data)"""
    with open(path) as f:
        content = f.read()

    # Find column header
    header = None
    for line in content.split('\n'):
        if line.strip().startswith('#') and 'Total-Boost-Energy-Potential' in line:
            header = line.strip().lstrip('#').strip()
            break
    if not header:
        raise ValueError(f"No header with Total-Boost-Energy-Potential in {path}")

    hdr_cols = [c.strip() for c in header.split(',')]
    total_boost_hdr_idx = hdr_cols.index('Total-Boost-Energy-Potential')

    # Parse data lines
    data_lines = []
    for line in content.split('\n'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 8:
            data_lines.append(parts)

    if not data_lines:
        raise ValueError(f"No data lines in {path}")

    # Detect offset: data may have leading empty element from leading tab
    n_hdr = len(hdr_cols)
    n_data = len(data_lines[0])
    offset = n_data - n_hdr
    boost_idx = total_boost_hdr_idx + offset

    log(f"  gamd.log: {n_hdr} header cols, {n_data} data cols, offset={offset}")
    log(f"  Total-Boost-Energy-Potential: header[{total_boost_hdr_idx}] → data[{boost_idx}]")

    boosts = np.array([float(line[boost_idx]) for line in data_lines])

    return boosts

def main():
    DATA = data_root()
    OUT = Path("results/analysis")
    OUT.mkdir(parents=True, exist_ok=True)

    systems = {
        "sys1_APO": (DATA / "sys1_APO" / "sys1_APO.prmtop",
                      DATA / "sys1_APO" / "output.dcd",
                      DATA / "sys1_APO" / "gamd.log"),
        "sys1_AZD5305": (DATA / "sys1_AZD5305" / "system.prmtop",
                          DATA / "sys1_AZD5305" / "output.dcd",
                          DATA / "sys1_AZD5305" / "gamd.log"),
    }

    for name, (top, traj, gamdlog) in systems.items():
        if not top.exists():
            log(f"SKIP {name}: no topology {top}")
            continue
        if not traj.exists():
            log(f"SKIP {name}: no trajectory {traj}")
            continue
        if not gamdlog.exists():
            log(f"SKIP {name}: no gamd.log {gamdlog}")
            continue

        log(f"\n{'='*40}")
        log(f"Processing {name}")

        # Extract CV from DCD
        log(f"Loading trajectory...")
        u = mda.Universe(str(top), str(traj))
        n_frames = len(u.trajectory)
        log(f"  {n_frames} frames, {len(u.atoms)} atoms")

        hd = u.select_atoms("resid 1:119 and name CA")
        art = u.select_atoms("resid 120:345 and name CA")
        log(f"  HD CA: {len(hd)}, ART CA: {len(art)}")

        n_analyze = n_frames // STRIDE
        cv_vals = np.zeros(n_analyze)

        for i, ts in enumerate(u.trajectory[::STRIDE]):
            cv_vals[i] = np.linalg.norm(hd.center_of_mass() - art.center_of_mass())
            if i % 50 == 0:
                log(f"  frame {i*STRIDE}: CV={cv_vals[i]:.2f}")

        log(f"  CV: mean={cv_vals.mean():.2f}±{cv_vals.std():.2f}, range=[{cv_vals.min():.2f}, {cv_vals.max():.2f}]")

        # Extract dV from gamd.log
        log(f"Parsing gamd.log...")
        boost_arr = parse_gamd_log(str(gamdlog))
        log(f"  {len(boost_arr)} boost entries")

        # Match CV frames with boost (both use stride)
        dv_vals = np.zeros(n_analyze)
        for i in range(n_analyze):
            frame_idx = i * STRIDE
            if frame_idx < len(boost_arr):
                dv_vals[i] = boost_arr[frame_idx]
            else:
                dv_vals[i] = boost_arr[-1]

        n_zero = np.sum(dv_vals < 0.001)
        log(f"  dV: mean={dv_vals.mean():.2f}±{dv_vals.std():.2f}, max={dv_vals.max():.2f}")
        log(f"  Zero boost frames: {n_zero}/{n_analyze} ({100*n_zero/n_analyze:.1f}%)")

        # Save
        cv_file = OUT / f"{name}_cv.dat"
        weight_file = OUT / f"{name}_weights.dat"

        np.savetxt(str(cv_file), cv_vals, fmt="%.4f")

        # 3-column weight format for PyReweighting
        beta = 1.0 / kT
        weights = np.column_stack([
            dv_vals * beta,            # col 0: dV/kT → exp(dV/kT)
            np.zeros(n_analyze),        # col 1: padding
            dv_vals                     # col 2: raw dV
        ])
        np.savetxt(str(weight_file), weights, fmt="%.6f")

        log(f"  Saved: {cv_file.name}, {weight_file.name}")

    log("\n[DONE]")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
