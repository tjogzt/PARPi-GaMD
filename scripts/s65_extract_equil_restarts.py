#!/usr/bin/env python3
"""s65_extract_equil_restarts.py — prepare ③ (multi-seed extension) locally.

Extracts the equilibration-endpoint frame from the two deposited S2
trajectories that lack a saved restart, and writes Amber rst7 files that can
seed new 26-ns GaMD productions on a GPU box (seeds 7 and 123), with fresh
velocities drawn at 300 K (shared initial structure, independent random
history).

  talazoparib : sys2_talazoparib/output.dcd, 1-ps frames, step = 500*k,
                frame 12000 == step 6,000,000 == 12 ns (equilibration end)
  veliparib   : sys2_veliparib/output.dcd, 10-ps frames, frame 1200 == 12 ns

Outputs
  data/restarts/{tala,veli}_s2_equil12ns.rst7
  data/restarts/restart_provenance.json
Sanity checks: atom count equals the prmtop; CV1/CV2 recomputed on the frame.
"""
import gzip
import hashlib
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("DATA_ROOT", "/Volumes/tjogzt4T/PARPi_data"))
OUT = ROOT / "data/restarts"
OUT.mkdir(parents=True, exist_ok=True)

CASES = [
    dict(key="tala", name="talazoparib",
         d=DATA / "sys2_talazoparib",
         prmtop="sys2_talazoparib.prmtop",
         frame=12000, step=6_000_000, cadence_ps=1.0),
    dict(key="veli", name="veliparib",
         d=DATA / "sys2_veliparib",
         prmtop="sys2_veliparib.prmtop",
         frame=1200, step=6_000_000, cadence_ps=10.0),
]


def write_rst7(path, xyz, box):
    """Minimal Amber rst7 (coordinates + box, no velocities)."""
    with open(path, "w") as f:
        f.write("equilibration endpoint restart (s65)\n")
        f.write(f"{len(xyz):5d}\n")
        flat = xyz.reshape(-1)
        for i in range(0, len(flat), 6):
            f.write("".join(f"{v:12.7f}" for v in flat[i:i + 6]) + "\n")
        if box is not None:
            f.write("".join(f"{v:12.7f}" for v in box) + "\n")


prov = {}
for c in CASES:
    import MDAnalysis as mda
    prm = c["d"] / c["prmtop"]
    u = mda.Universe(str(prm), str(c["d"] / "output.dcd"))
    n = len(u.trajectory)
    if c["frame"] >= n:
        raise SystemExit(f"{c['key']}: frame {c['frame']} >= n_frames {n}")
    u.trajectory[c["frame"]]
    xyz = u.atoms.positions.copy()
    box = u.dimensions.copy() if u.dimensions is not None else None
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    hd = u.select_atoms("resid 662:787 and name CA")
    art = u.select_atoms("resid 788:1011 and name CA")
    cv1 = float(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
    cv2 = float(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    out = OUT / f"{c['key']}_s2_equil12ns.rst7"
    write_rst7(out, xyz, box)
    prov[c["key"]] = dict(
        system=c["name"], source=str(c["d"] / "output.dcd"),
        prmtop=str(prm), prmtop_sha256=hashlib.sha256(prm.read_bytes()).hexdigest()[:16],
        n_frames=n, frame_index=c["frame"], step=c["step"],
        time_ns=c["frame"] * c["cadence_ps"] / 1000.0,
        n_atoms=int(len(u.atoms)), box=[float(v) for v in box] if box is not None else None,
        cv1_A=round(cv1, 2), cv2_A=round(cv2, 2), rst7=str(out))
    print(f"{c['key']}: frame {c['frame']} (t={prov[c['key']]['time_ns']:.2f} ns) "
          f"-> {out.name}  n_atoms={len(u.atoms)}  CV1={cv1:.1f}  CV2={cv2:.1f}")

(OUT / "restart_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
print("wrote", OUT / "restart_provenance.json")
