#!/usr/bin/env python3
#
# Purpose:  docking receptor preparation + box definition
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   Usage: python prepare_receptor.py <PDB_ID> <chain> <ligand_resname> [padding]
# Outputs:  - receptor -> obabel -xr -> receptor.pdbqt (Gasteiger, rigid)
# Depends:  numpy, Bio
# Run:      python3 scripts/s11_prepare_receptor.py

"""
prepare_receptor.py -- docking receptor preparation + box definition

Extracts from the PDB: the receptor (protein, selected chain, water/buffer removed) + the reference ligand (crystallographically bound inhibitor).
- receptor -> obabel -xr -> receptor.pdbqt (Gasteiger, rigid)
- reference ligand -> ref_ligand.pdb/.pdbqt (for redocking RMSD validation)
- docking box defined from the reference-ligand centroid (ligand bounding box + padding)

Usage: python prepare_receptor.py <PDB_ID> <chain> <ligand_resname> [padding]
Example: python prepare_receptor.py 6VKK A RPB 8.0
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = ROOT / "data" / "00_raw" / "pdb"
REC_DIR = ROOT / "docking" / "receptors"
BOX_DIR = ROOT / "docking" / "box"

WATER = {"HOH", "WAT"}
BUFFERS = {"SO4", "GOL", "EDO", "PO4", "ACT", "DMS", "PEG", "1PE",
           "FMT", "MPD", "GLC", "FRU", "CO3", "NA", "CL", "K"}
KEEP_METAL = {"ZN", "MG", "MN", "CA"}  # structural metals retained


class ProteinSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, res):
        het = res.id[0]
        rn = res.get_resname().strip()
        if het == " ":
            return 1  # standard residue
        if rn in KEEP_METAL:
            return 1
        return 0  # water/ligand/buffer


class LigandSelect(Select):
    def __init__(self, chain_id, resname):
        self.chain_id, self.resname = chain_id, resname

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, res):
        return res.get_resname().strip() == self.resname


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def main() -> int:
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    pid, chain, lig = sys.argv[1], sys.argv[2], sys.argv[3]
    pad = float(sys.argv[4]) if len(sys.argv) > 4 else 8.0
    src = PDB_DIR / f"{pid}.pdb"
    if not src.exists():
        sys.exit(f"missing {src}")
    REC_DIR.mkdir(parents=True, exist_ok=True)
    BOX_DIR.mkdir(parents=True, exist_ok=True)

    s = PDBParser(QUIET=True).get_structure(pid, str(src))
    # keep only the first model
    model = next(iter(s))
    io = PDBIO()

    # receptor
    rec_pdb = REC_DIR / f"{pid}_{chain}_receptor.pdb"
    io.set_structure(model)
    io.save(str(rec_pdb), ProteinSelect(chain))

    # reference ligand
    ref_pdb = REC_DIR / f"{pid}_{chain}_{lig}_ref.pdb"
    io.set_structure(model)
    io.save(str(ref_pdb), LigandSelect(chain, lig))

    # ligand coordinates -> box
    coords = []
    for line in ref_pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not coords:
        sys.exit(f"no atoms for ligand {lig} in chain {chain}")
    arr = np.array(coords)
    center = arr.mean(axis=0)
    size = (arr.max(axis=0) - arr.min(axis=0)) + 2 * pad

    # receptor -> pdbqt (rigid)
    rec_pdbqt = REC_DIR / f"{pid}_{chain}_receptor.pdbqt"
    run(["obabel", str(rec_pdb), "-O", str(rec_pdbqt),
         "-xr", "--partialcharge", "gasteiger"])
    # reference ligand -> pdbqt (for redocking)
    ref_pdbqt = REC_DIR / f"{pid}_{chain}_{lig}_ref.pdbqt"
    run(["obabel", str(ref_pdb), "-O", str(ref_pdbqt),
         "--partialcharge", "gasteiger"])

    box = BOX_DIR / f"{pid}_{chain}_{lig}.box"
    box.write_text(
        f"center_x = {center[0]:.3f}\ncenter_y = {center[1]:.3f}\ncenter_z = {center[2]:.3f}\n"
        f"size_x = {size[0]:.1f}\nsize_y = {size[1]:.1f}\nsize_z = {size[2]:.1f}\n")

    print(f"[OK] receptor  -> {rec_pdbqt.name} ({rec_pdbqt.stat().st_size} B)")
    print(f"[OK] ref ligand-> {ref_pdbqt.name} ({len(coords)} atoms)")
    print(f"[OK] box center=({center[0]:.1f},{center[1]:.1f},{center[2]:.1f}) "
          f"size=({size[0]:.1f},{size[1]:.1f},{size[2]:.1f})  -> {box.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
