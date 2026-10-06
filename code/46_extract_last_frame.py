#!/usr/bin/env python3
"""Extract last-frame coordinates from a legacy S2 DCD to .npy (local; for cloud parmed PDB assembly)

Purpose:   One-off helper: pull the final frame of the AZD5305/veliparib S2
           trajectories into .npy files used for ligand-pose PDB assembly on
           the cloud machine (parmed). Not part of the analysis pipeline.
Inputs:    $DATA_ROOT/sys2_<drug>/output.dcd (CHARMM/X-PLOR DCD)
Outputs:   results/<drug>_last_xyz.npy (natoms x 3 float32)
Depends:   MDAnalysis; common.paths (DATA_ROOT)

Created:  2026-09-15 (header standardised 2026-10-05)
Run:      python3 code/46_extract_last_frame.py   (from the repository root)
"""
import os
import sys
os.makedirs("results", exist_ok=True)
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root

import MDAnalysis as mda


def parse_dcd_last(path):
    """Return (n_frames, natoms, last_frame_xyz) using MDAnalysis."""
    u = mda.Universe(path)
    n_frames = len(u.trajectory)
    natoms = len(u.atoms)
    u.trajectory[-1]
    return n_frames, natoms, u.atoms.positions.astype(np.float32)


for drug, dcd, out in [
    ('AZD5305', str(data_root() / 'sys2_AZD5305' / 'output.dcd'), 'results/azd5305_last_xyz.npy'),
    ('veliparib', str(data_root() / 'sys2_veliparib' / 'output.dcd'), 'results/veliparib_last_xyz.npy'),
]:
    n, natoms, xyz = parse_dcd_last(dcd)
    np.save(out, xyz)
    print(f'{drug}: {n} frames, last frame {natoms} atoms -> {out} ({os.path.getsize(out)} bytes)')
