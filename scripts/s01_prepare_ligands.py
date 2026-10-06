#!/usr/bin/env python3
"""
prepare_ligands.py -- ligand 3D preparation (docking-ready)

Reads isomeric SMILES from data/01_curated/trapping_seed_dataset.csv,
RDKit add hydrogens -> ETKDGv3 multi-conformer embedding -> MMFF94 optimization -> lowest-energy conformer selected,
writes SDF + PDB to data/01_modeling/docking_ligands/, then converts to PDBQT with obabel (Gasteiger charges).

Stereochemistry is preserved from the isomeric SMILES (talazoparib dual stereocenters are critical).

Purpose:  Ligand 3D preparation (docking-ready): isomeric SMILES -> RDKit conformers -> SDF/PDB/PDBQT.
Created:  2026-09-16 (header standardised 2026-10-05)
Depends:  __future__, pandas, rdkit, rdkit.Chem
Run:      python3 scripts/s01_prepare_ligands.py   (from the repository root)
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "data" / "01_curated" / "trapping_seed_dataset.csv"
OUT = ROOT / "data" / "01_modeling" / "docking_ligands"
N_CONF = 20
SEED = 42


def embed_best(mol: Chem.Mol) -> tuple[Chem.Mol, float] | tuple[None, None]:
    """Multi-conformer embedding + MMFF optimization, returns the lowest-energy conformer."""
    molH = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.useSmallRingTorsions = True
    cids = AllChem.EmbedMultipleConfs(molH, numConfs=N_CONF, params=params)
    if not cids:
        return None, None
    energies = []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(
            molH, AllChem.MMFFGetMoleculeProperties(molH), confId=cid)
        if ff is None:
            continue
        ff.Minimize(maxIts=2000)
        energies.append((ff.CalcEnergy(), cid))
    if not energies:
        return None, None
    energies.sort()
    best_e, best_cid = energies[0]
    # keep only the best conformer
    keep = Chem.Mol(molH)
    keep.RemoveAllConformers()
    keep.AddConformer(molH.GetConformer(best_cid), assignId=True)
    return keep, best_e


def to_pdbqt(sdf: Path, pdbqt: Path) -> bool:
    try:
        subprocess.run(
            ["obabel", str(sdf), "-O", str(pdbqt), "--partialcharge", "gasteiger"],
            check=True, capture_output=True, text=True)
        return pdbqt.exists() and pdbqt.stat().st_size > 0
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  [obabel FAIL] {e}", file=sys.stderr)
        return False


def main() -> int:
    if not CSV.exists():
        sys.exit(f"missing {CSV}")
    df = pd.read_csv(CSV)
    OUT.mkdir(parents=True, exist_ok=True)
    ok = fail = 0
    summary = []
    for _, row in df.iterrows():
        name, smi = row["name"], row["isomeric_smiles"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"[FAIL] {name}: bad SMILES"); fail += 1; continue
        confmol, energy = embed_best(mol)
        if confmol is None:
            print(f"[FAIL] {name}: embed failed"); fail += 1; continue
        confmol.SetProp("_Name", name)
        sdf = OUT / f"{name}.sdf"
        pdb = OUT / f"{name}.pdb"
        pdbqt = OUT / f"{name}.pdbqt"
        Chem.MolToMolFile(confmol, str(sdf))
        Chem.MolToPDBFile(confmol, str(pdb))
        pq = to_pdbqt(sdf, pdbqt)
        # stereocenter verification
        sc = Chem.FindMolChiralCenters(confmol, useLegacyImplementation=False,
                                       includeUnassigned=True)
        flag = "OK" if pq else "PDBQT_FAIL"
        print(f"[{flag}] {name:12s} E={energy:8.1f} kcal/mol  chiral={len(sc)}  -> {sdf.name}")
        summary.append((name, round(energy, 1), len(sc), pq))
        ok += 1 if pq else 0
        fail += 0 if pq else 1

    print(f"\n[DONE] ligands prepared: {ok} ok / {fail} fail  -> {OUT}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
