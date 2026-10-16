#!/usr/bin/env python3
#
# Purpose:  veli_prot sys2 CV extraction — mirrors extract_s2_cv.py, veli_prot paths.
# Created:  2026-10 (header standardised 2026-10-11)
# Run:      python3 code/pipeline/p09_extract_veli_prot_cv.py   (from the repository root)

"""veli_prot sys2 CV extraction — mirrors extract_s2_cv.py, veli_prot paths.

CV1 = protein CA COM − DNA heavy-atom COM (Å); CV2 = HD(662-787) CA − ART(788-1011) CA (Å).
Production mask read from gamd_out/production-start-step.txt (6.0e6).
Outputs: analysis_CV1.dat, analysis_CV2.dat, analysis_weights.dat (= (beta*dV_D, 0, dV_D)).
Usage: python veli_prot_s2_extract.py <base_dir e.g. .../runs_veli_prot/sys2_veliparib_prot>
"""
import sys
from pathlib import Path
import mdtraj as md
import numpy as np

base = Path(sys.argv[1])
beta = 1.0 / 0.596  # 300 K, matches frozen protocol
mask = int((base / "gamd_out/production-start-step.txt").read_text().strip())

top = md.load_prmtop(str(base / "build/veliparib_best_prot.prmtop"))
hd = top.select("resid 662 to 787 and name CA")
art = top.select("resid 788 to 1011 and name CA")
prot = top.select("protein and name CA")
dna = top.select(
    "(resname DA or resname DT or resname DG or resname DC or "
    "resname DG5 or resname DC3) and (not element H)"
)
assert len(hd) == 126, f"HD CA={len(hd)} != 126"
assert len(art) == 224, f"ART CA={len(art)} != 224"
assert len(dna) > 900, f"DNA heavy={len(dna)} unexpected"
print(f"HD={len(hd)} ART={len(art)} ProtCA={len(prot)} DNAheavy={len(dna)} mask={mask}", flush=True)

cv1c, cv2c = [], []
for chunk in md.iterload(str(base / "gamd_out/output.dcd"), top=top, chunk=100):
    pc = chunk.xyz[:, prot, :].mean(axis=1)
    dc = chunk.xyz[:, dna, :].mean(axis=1)
    cv1c.append(np.linalg.norm(pc - dc, axis=1) * 10.0)
    cv2c.append(
        np.linalg.norm(
            chunk.xyz[:, hd, :].mean(axis=1) - chunk.xyz[:, art, :].mean(axis=1), axis=1
        ) * 10.0
    )
cv1 = np.concatenate(cv1c)
cv2 = np.concatenate(cv2c)
print(f"frames={len(cv1)}", flush=True)

log = np.loadtxt(base / "gamd_out/gamd.log")
steps, mode, dbe = log[:, 1], log[:, 0], log[:, 7]
n = min(len(cv1), len(log))
prod = (steps[:n] >= mask) & (mode[:n] == 1)
print(f"log_rows={len(log)} production_frames={int(prod.sum())}", flush=True)
assert prod.sum() >= 500, f"too few production frames: {prod.sum()}"

cv1p, cv2p, dbep = cv1[:n][prod], cv2[:n][prod], dbe[:n][prod]
np.savetxt(base / "analysis_CV1.dat", cv1p, fmt="%.4f")
np.savetxt(base / "analysis_CV2.dat", cv2p, fmt="%.4f")
np.savetxt(
    base / "analysis_weights.dat",
    np.column_stack([dbep * beta, np.zeros(len(dbep)), dbep]),
    fmt="%.4f",
)
print(
    f"SAVED CV1={cv1p.mean():.2f}±{cv1p.std():.2f} "
    f"CV2={cv2p.mean():.2f}±{cv2p.std():.2f} DBE={dbep.mean():.2f}±{dbep.std():.2f}",
    flush=True,
)
print("ALL DONE", flush=True)
