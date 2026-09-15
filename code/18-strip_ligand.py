#!/usr/bin/env python3
"""
18-strip_ligand.py — Strip ligand atoms from all S2 DCDs.
Strips atoms from index APO_N (285775) onwards, keeping only
protein+DNA+water+ions matching APO prmtop.
"""
import MDAnalysis as mda
from pathlib import Path
import warnings
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root
warnings.filterwarnings("ignore")

DATA = data_root()
APO_N = 285775

systems = {
    "sys2_AZD5305":    {"stride": 10,  "atoms": 285792},
    "sys2_veliparib":  {"stride": 10,  "atoms": 285788},
    "sys2_niraparib":  {"stride": 100, "atoms": 285793},
    "sys2_olaparib":   {"stride": 100, "atoms": 285799},
    "sys2_rucaparib":  {"stride": 100, "atoms": 285793},
}

for name, cfg in systems.items():
    stride = cfg["stride"]
    n_atoms = cfg["atoms"]
    lig_atoms = n_atoms - APO_N

    dcd_in = DATA / name / "output.dcd"
    dcd_out = DATA / name / "output_stripped.dcd"

    if not dcd_in.exists():
        print(f"[SKIP] {name}: no DCD")
        continue

    # Remove old stripped file if stride changed
    if dcd_out.exists():
        # Check if already done with correct stride
        old_sz = dcd_out.stat().st_size
        if old_sz > 0:
            print(f"[SKIP] {name}: exists ({old_sz/(1024**3):.2f} GB)")
            continue

    sz_in = dcd_in.stat().st_size / (1024**3)
    print(f"\n{'='*50}")
    print(f"{name}: {sz_in:.1f} GB, {n_atoms} atoms, "
          f"stripping {lig_atoms} ligand atoms, stride={stride}")

    # Create universe matching DCD atom count
    u = mda.Universe.empty(n_atoms=n_atoms, trajectory=True)
    u.load_new(str(dcd_in))
    n_frames = len(u.trajectory)
    n_out = (n_frames + stride - 1) // stride
    print(f"  {n_frames} frames → ~{n_out} output frames")

    ag = u.atoms[:APO_N]  # keep only protein+DNA+water+ions

    with mda.Writer(str(dcd_out), n_atoms=APO_N) as W:
        cnt = 0
        for i, ts in enumerate(u.trajectory):
            if i % stride != 0:
                continue
            W.write(ag)
            cnt += 1
            if cnt % 200 == 0:
                print(f"    {cnt}/{n_out}")

    sz_out = dcd_out.stat().st_size / (1024**3)
    print(f"  DONE: {dcd_out.name} ({sz_out:.2f} GB, {cnt} frames)")

print("\n" + "="*50)
print("All stripped DCDs ready!")
