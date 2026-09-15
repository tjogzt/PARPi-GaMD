#!/usr/bin/env python3
"""
20-extract_ca_dcd.py — Extract CA-only DCDs using MDAnalysis.
Much more reliable than cpptraj for our stripped trajectories.
"""
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import warnings
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root
warnings.filterwarnings("ignore")

DATA = data_root()
OUT = data_root() / "md_analysis"
OUT.mkdir(parents=True, exist_ok=True)
APO_TOP = DATA / "sys2_APO" / "sys2_APO.prmtop"

systems = {}
for name in ["sys2_APO", "sys2_talazoparib"]:
    top = DATA / name / f"{name}.prmtop"
    dcd = DATA / name / "output.dcd"
    if top.exists() and dcd.exists():
        systems[name] = (str(top), str(dcd))

for name in ["sys2_AZD5305", "sys2_veliparib", "sys2_niraparib",
             "sys2_olaparib", "sys2_rucaparib"]:
    dcd = DATA / name / "output_stripped.dcd"
    if dcd.exists() and dcd.stat().st_size > 1000:
        systems[name] = (str(APO_TOP), str(dcd))

for name, (top_path, dcd_path) in systems.items():
    ca_dcd = OUT / f"{name}_ca.dcd"
    if ca_dcd.exists() and ca_dcd.stat().st_size > 102400:
        print(f"[SKIP] {name}: exists ({ca_dcd.stat().st_size/(1024**2):.1f} MB)")
        continue

    print(f"\n{name}: extracting CA atoms...")
    u = mda.Universe(top_path, dcd_path)
    ca = u.select_atoms("protein and name CA")
    n_ca = len(ca)
    n_frames = len(u.trajectory)
    expected_mb = n_frames * n_ca * 3 * 4 / (1024**2)
    print(f"  {n_frames} frames × {n_ca} CA = ~{expected_mb:.1f} MB")

    with mda.Writer(str(ca_dcd), n_atoms=n_ca) as W:
        for i, ts in enumerate(u.trajectory):
            W.write(ca)
            if (i + 1) % 5000 == 0:
                print(f"    {i+1}/{n_frames}")

    actual = ca_dcd.stat().st_size / (1024**2)
    print(f"  DONE: {actual:.1f} MB")

    # Free memory
    del u

print("\nAll CA DCDs ready!")
