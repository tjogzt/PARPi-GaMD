#!/usr/bin/env python3
"""
control_tala_rmsf.py — legacy-protocol RMSF control for the talazoparib S2 DCD

Purpose:   Diagnostic: whole-protein-aligned, all-frames RMSF for the
           talazoparib S2 trajectory (protocol check only; produces no
           manuscript numbers).
           NOTE: the "legacy table" values quoted by earlier versions of this
           script (HD 3.12 / ART 3.65) were produced by the July-2025 run of
           19-rmsf_all_systems.py from an earlier version of the on-disk DCD.
           The current output.dcd (26,000 frames) yields ~6.1/~8.6 A under the
           same whole-protein protocol (production frames give the same), so
           the legacy values are not reproducible from the present trajectory
           file. The manuscript RMSF values come from code/recompute_all_s2.py
           (production frames, domain self-alignment; talazoparib HD 3.75 /
           ART 3.88), which is reproducible from the current data.
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
print(f'HD RMSF = {rmsf[hd].mean():.2f} A  (manuscript unified protocol: 3.75)')
print(f'ART RMSF = {rmsf[art].mean():.2f} A  (manuscript unified protocol: 3.88)')
