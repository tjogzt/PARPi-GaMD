#!/usr/bin/env python3
"""
control_tala_rmsf.py — legacy talazoparib DCD control

Purpose:   Reproduce the legacy-table RMSF values (HD 3.12 / ART 3.65) with the
           old protocol applied to the original talazoparib S2 trajectory.
Inputs:    $DATA_ROOT/sys2_talazoparib/{sys2_talazoparib.prmtop, output.dcd}
Outputs:   console RMSF summary
Depends:   MDAnalysis, numpy; DATA_ROOT environment variable
"""
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.kabsch import kabsch_align
from common.paths import data_root

DR = data_root()
u = mda.Universe(str(DR / 'sys2_talazoparib' / 'sys2_talazoparib.prmtop'),
                 str(DR / 'sys2_talazoparib' / 'output.dcd'))
ca = u.select_atoms('protein and name CA')
resids = ca.residues.resids
n = len(u.trajectory)
print(f'frames {n}, CA atoms {len(ca)}, resid {resids[0]}-{resids[-1]}')

u.trajectory[0]
ref = ca.positions.copy()
sp = np.zeros((len(ca), 3))
sp2 = np.zeros((len(ca), 3))
for ts in u.trajectory:
    al = kabsch_align(ca.positions.copy(), ref)
    sp += al
    sp2 += al**2
mp = sp / n
rmsf = np.sqrt((sp2 / n - mp**2).sum(axis=1))
hd = (resids >= 662) & (resids <= 787)
art = (resids >= 788) & (resids <= 1014)
print(f'HD RMSF = {rmsf[hd].mean():.2f} (legacy table 3.12)')
print(f'ART RMSF = {rmsf[art].mean():.2f} (legacy table 3.65)')
