#!/usr/bin/env python3
"""Extract production start frames (frame 0) from the three archived S2 DCDs
as .rst7 restart files for the rerun (same prmtop, restart from frame 0).

Frame 0 of output.dcd is the equilibrated structure at production start
(deploy protocol: NVT 0.2 ns + NPT 1.0 ns then GaMD production). The rerun
uses the same prmtop and re-equilibrates from this structure (gamdRunner
conventional-MD stage regenerates velocities), so only coordinates + box are
needed.

Usage: python3 extract_rerun_starts.py [DATA_ROOT]
Output: results/rerun_starts/sys2_<lig>_start.rst7
"""
import sys
from pathlib import Path

import numpy as np

DATA_ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/Volumes/tjogzt4T/PARPi_data")

import MDAnalysis as mda

OUT = Path(__file__).resolve().parents[1] / "results" / "rerun_starts"
OUT.mkdir(parents=True, exist_ok=True)

LIGANDS = ["niraparib", "olaparib", "rucaparib"]

for lig in LIGANDS:
    d = DATA_ROOT / f"sys2_{lig}"
    prmtop = d / f"sys2_{lig}.prmtop"
    dcd = d / "output.dcd"
    u = mda.Universe(str(prmtop), str(dcd))
    u.trajectory[0]
    xyz = u.atoms.positions.astype(np.float64)
    box = u.dimensions  # [A, B, C, alpha, beta, gamma]
    natoms = len(u.atoms)

    # Amber rst7: title line, natoms, then 2 atoms per line (3f12.7 each,
    # 72 chars/line), then a 6-value box line. OpenMM's AmberAsciiRestart
    # parses fixed 12-char fields, so pairs per line are mandatory.
    out = OUT / f"sys2_{lig}_start.rst7"
    with open(out, "w") as f:
        f.write(f"rerun start frame 0 of sys2_{lig} (natoms {natoms})\n")
        f.write(f"{natoms:8d}\n")
        for i in range(0, natoms - 1, 2):
            f.write(f"{xyz[i, 0]:12.7f}{xyz[i, 1]:12.7f}{xyz[i, 2]:12.7f}"
                    f"{xyz[i+1, 0]:12.7f}{xyz[i+1, 1]:12.7f}{xyz[i+1, 2]:12.7f}\n")
        if natoms % 2 == 1:
            i = natoms - 1
            f.write(f"{xyz[i, 0]:12.7f}{xyz[i, 1]:12.7f}{xyz[i, 2]:12.7f}\n")
        f.write(f"{box[0]:12.7f}{box[1]:12.7f}{box[2]:12.7f}"
                f"{box[3]:12.7f}{box[4]:12.7f}{box[5]:12.7f}\n")
    print(f"{lig}: {natoms} atoms, box {box[:3].round(2)} -> {out}")
print("DONE")
