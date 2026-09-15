#!/usr/bin/env python3
"""
05-sys1_cv_extract.py — Extract HD-ART COM distance CV for sys1 GaMD trajectories.
sys1 = CAT domain only (residues 1-345, UniProt 662-1011).
CV: center-of-mass distance between HD (resid 1-119 CA) and ART (resid 120-345 CA).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np

def log(msg: str) -> None:
    print(f"[cv] {msg}", flush=True)

def main():
    import MDAnalysis as mda

    DATA = Path("/Volumes/tjogzt4T/PARPi_data")
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
        if not top.exists() or not traj.exists():
            log(f"SKIP {name}: missing files")
            continue

        log(f"Loading {name}...")
        u = mda.Universe(str(top), str(traj))
        n_frames = len(u.trajectory)
        log(f"  frames={n_frames} atoms={len(u.atoms)}")

        # Selections
        hd = u.select_atoms("resid 1:119 and name CA")
        art = u.select_atoms("resid 120:345 and name CA")
        log(f"  HD CA: {len(hd)}, ART CA: {len(art)}")

        # Extract CV, dV, and boost energies
        stride = 100
        n_analyze = n_frames // stride
        cv_vals = np.zeros(n_analyze)
        dv_vals = np.zeros(n_analyze)

        # Read gamd.log for boost potential
        boost_data = []
        with open(gamdlog) as f:
            # Skip comment lines to find the actual header
            cols = None
            for line in f:
                line = line.strip()
                if line.startswith("#") and "ntwx" in line:
                    cols = line.replace("# ", "").split(",")
                    break
            if cols is None:
                raise ValueError("Cannot find valid header in gamd.log")

            # Find boost energy column
            dih_idx = None
            for i, c in enumerate(cols):
                if "Boost-Energy" in c:
                    dih_idx = i
                    break
            if dih_idx is None:
                for i, c in enumerate(cols):
                    if "Boost-Energy-Potential" in c:
                        dih_idx = i
                        break
            if dih_idx is None:
                raise ValueError(f"Cannot find boost energy column in: {cols}")

            log(f"  gamd.log boost column: {cols[dih_idx]} (idx={dih_idx})")
            for line in f:
                if line.startswith("#"):
                    continue
                # Data lines use whitespace (tabs), header uses commas
                parts = line.strip().split()
                if len(parts) > dih_idx:
                    try:
                        boost_data.append(float(parts[dih_idx]))
                    except ValueError:
                        continue

        boost_arr = np.array(boost_data)
        log(f"  gamd.log entries: {len(boost_arr)}")

        for i, ts in enumerate(u.trajectory[::stride]):
            hd_com = hd.center_of_mass()
            art_com = art.center_of_mass()
            cv_vals[i] = np.linalg.norm(hd_com - art_com)
            # dV from boost potential at this frame
            frame_idx = i * stride
            if frame_idx < len(boost_arr):
                dv_vals[i] = boost_arr[frame_idx]
            if i % 50 == 0:
                log(f"  frame {frame_idx}: CV={cv_vals[i]:.2f} dV={dv_vals[i]:.2f}")

        log(f"  CV range: {cv_vals.min():.2f} - {cv_vals.max():.2f} Å")
        log(f"  CV mean: {cv_vals.mean():.2f} ± {cv_vals.std():.2f} Å")
        log(f"  dV mean: {dv_vals.mean():.2f} ± {dv_vals.std():.2f} kcal/mol")

        # Save as: CV-only file + 3-column weight file for PyReweighting
        cv_file = OUT / f"{name}_cv.dat"
        weight_file = OUT / f"{name}_weights.dat"

        np.savetxt(str(cv_file), cv_vals, fmt="%.4f")

        # PyReweighting format: dV*beta, 0.0, dV  (beta = 1/kT, kT≈0.596 at 300K)
        beta = 1.0 / 0.596
        weights = np.column_stack([dv_vals * beta,
                                    np.zeros(n_analyze),
                                    dv_vals])
        np.savetxt(str(weight_file), weights, fmt="%.4f")

        log(f"  Saved: {cv_file.name} ({n_analyze} frames)")
        log(f"  Saved: {weight_file.name}")

    log("\n[DONE]")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
