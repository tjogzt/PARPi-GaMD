#!/usr/bin/env python3
"""
prepare_ligands.py — 配体 3D 准备 (对接就绪)

读 data/01_curated/trapping_seed_dataset.csv 的 isomeric SMILES,
RDKit 加氢 -> ETKDGv3 多构象嵌入 -> MMFF94 优化 -> 选最低能构象,
输出 SDF + PDB (docking/ligands/), 再用 obabel 转 PDBQT (Gasteiger 电荷)。

立体化学从 isomeric SMILES 保留 (talazoparib 双立体中心关键)。
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
OUT = ROOT / "docking" / "ligands"
N_CONF = 20
SEED = 42


def embed_best(mol: Chem.Mol) -> tuple[Chem.Mol, float] | tuple[None, None]:
    """多构象嵌入 + MMFF 优化, 返回最低能构象。"""
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
    # 仅保留最优构象
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
        # 立体中心核验
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
