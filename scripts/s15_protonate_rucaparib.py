#!/usr/bin/env python3
#
# Purpose:  Rucaparib protonation: build the correct structure from the protonated SMILES, map the docked pose coordinates by MCS. Protonation site = benzylamine N (CH2-NH-CH3 -> CH2-NH2(+)-CH3, net charge +1; pKa ~8.4).
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   (see docstring)
# Outputs:  protonated pose files (S1 .sdf / S2 .mol2; see docstring)
# Depends:  rdkit
# Run:      python3 scripts/s15_protonate_rucaparib.py

"""Rucaparib protonation: build the correct structure from the protonated SMILES, map the docked pose coordinates by MCS.
Protonation site = benzylamine N (CH2-NH-CH3 -> CH2-NH2(+)-CH3, net charge +1; pKa ~8.4).
Outputs:
  docking/results_system2/mol2/rucaparib_best_prot.mol2  (S2 pose coordinates, correct structure, +1)
  docking/results/rucaparib_docked_prot.sdf              (S1 pose coordinates, correct structure, +1)
"""
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

ROOT = Path(__file__).resolve().parents[1]

# rucaparib protonated: benzylamine [NH2+]
SMILES_PROT = "C[NH2+]Cc1ccc(-c2[nH]c3cc(F)cc4c3c2CCNC4=O)cc1"

def build_correct_mol():
    m = Chem.MolFromSmiles(SMILES_PROT)
    assert m is not None
    m = Chem.AddHs(m)
    for a in m.GetAtoms():
        if a.GetSymbol() == "H":
            a.SetNoImplicit(True)
    Chem.SanitizeMol(m)
    ti = None
    for a in m.GetAtoms():
        if a.GetSymbol() == "N" and a.GetFormalCharge() == 1:
            ti = a.GetIdx()
            break
    assert ti is not None, "protonated amine N not found"
    return m, ti

def map_pose_coords(ref_mol, pose_path):
    pose = Chem.MolFromMolFile(str(pose_path), removeHs=False, sanitize=False)
    for a in pose.GetAtoms():
        a.SetNoImplicit(False)
    Chem.SanitizeMol(pose)
    pose_heavy = Chem.RemoveHs(pose)
    ref_heavy = Chem.RemoveHs(Chem.Mol(ref_mol))
    res = rdFMCS.FindMCS([pose_heavy, ref_heavy], atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareOrder, timeout=10)
    patt = Chem.MolFromSmarts(res.smartsString)
    ref_match = ref_heavy.GetSubstructMatch(patt)
    pose_match = pose_heavy.GetSubstructMatch(patt)
    conf_pose = pose.GetConformer()
    conf_ref = Chem.Conformer(ref_heavy.GetNumAtoms())
    pos_map = dict(zip(ref_match, pose_match))
    ref_un = [i for i in range(ref_heavy.GetNumAtoms()) if i not in pos_map]
    pose_un = [i for i in range(pose_heavy.GetNumAtoms()) if i not in pos_map.values()]
    assert len(ref_un) == len(pose_un), f"unmatched count mismatch {len(ref_un)} vs {len(pose_un)}"
    assert all(ref_heavy.GetAtomWithIdx(i).GetSymbol() == pose_heavy.GetAtomWithIdx(j).GetSymbol()
               for i, j in zip(ref_un, pose_un)), "element mismatch in tail mapping"
    for i, j in zip(ref_un, pose_un):
        pos_map[i] = j
    for i in range(ref_heavy.GetNumAtoms()):
        p = conf_pose.GetAtomPosition(pos_map[i]) if i in pos_map else Chem.rdGeometry.Point3D(0, 0, 0)
        conf_ref.SetAtomPosition(i, p)
    ref_heavy.RemoveAllConformers()
    ref_heavy.AddConformer(conf_ref)
    ref_full = Chem.AddHs(ref_heavy, addCoords=True)
    for a in ref_full.GetAtoms():
        if a.GetSymbol() == "H":
            a.SetNoImplicit(True)
    Chem.SanitizeMol(ref_full)
    props = AllChem.MMFFGetMoleculeProperties(ref_full)
    ff = AllChem.MMFFGetMoleculeForceField(ref_full, props)
    ff.Minimize(maxIts=1000)
    conf_ff = ref_full.GetConformer()

    def centroid(mol, ids):
        x = y = z = 0.0
        for i in ids:
            p = mol.GetConformer().GetAtomPosition(i)
            x += p.x; y += p.y; z += p.z
        n = len(ids)
        return Chem.rdGeometry.Point3D(x/n, y/n, z/n)
    c_pose = centroid(ref_heavy, list(pos_map.keys()))
    c_ff = centroid(ref_full, list(pos_map.keys()))
    d = Chem.rdGeometry.Point3D(c_pose.x - c_ff.x, c_pose.y - c_ff.y, c_pose.z - c_ff.z)
    conf_out = Chem.Conformer(ref_full.GetNumAtoms())
    for i in range(ref_full.GetNumAtoms()):
        p = conf_ff.GetAtomPosition(i)
        conf_out.SetAtomPosition(i, Chem.rdGeometry.Point3D(p.x + d.x, p.y + d.y, p.z + d.z))
    ref_full.RemoveAllConformers()
    ref_full.AddConformer(conf_out)
    Chem.SanitizeMol(ref_full)
    print(f"  MCS matched {len(pos_map)} heavy atoms; full MMFF + centroid-aligned ({ref_full.GetNumAtoms()} atoms)")
    return ref_full

def write_mol2(mol, ti, out):
    conf = mol.GetConformer()
    lines = ["@<TRIPOS>MOLECULE", "RUCA_PROT",
             f"  {mol.GetNumAtoms()} {mol.GetNumBonds()} 0 0 0", "SMALL", "USER_CHARGES", "",
             "@<TRIPOS>ATOM"]
    for a in mol.GetAtoms():
        i = a.GetIdx()
        pos = conf.GetAtomPosition(i)
        sym = a.GetSymbol()
        if a.GetIsAromatic():
            t = {'C': 'C.ar', 'N': 'N.ar', 'O': 'O.ar'}.get(sym, sym)
        elif sym == 'N' and i == ti:
            t = 'N.4'
        else:
            t = {'C': 'C.3', 'N': 'N.3', 'O': 'O.2', 'H': 'H'}.get(sym, sym)
        lines.append(f"{i+1:>6d} {sym:<4s} {pos.x:10.4f} {pos.y:10.4f} {pos.z:10.4f} "
                     f"{t:<6s} 1  UNL1        {a.GetFormalCharge():.4f}")
    lines.append("@<TRIPOS>BOND")
    for b in mol.GetBonds():
        bt = {Chem.BondType.AROMATIC: 'ar', Chem.BondType.SINGLE: '1',
              Chem.BondType.DOUBLE: '2'}.get(b.GetBondType(), '1')
        lines.append(f"{b.GetIdx()+1:>6d} {b.GetBeginAtomIdx()+1:5d} {b.GetEndAtomIdx()+1:5d} {bt}")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")

def main():
    m2, ti2 = build_correct_mol()
    m2 = map_pose_coords(m2, ROOT / "docking/results_system2/rucaparib_best.sdf")
    write_mol2(m2, ti2, ROOT / "docking/results_system2/mol2/rucaparib_best_prot.mol2")
    w = Chem.SDWriter(str(ROOT / "docking/results_system2/rucaparib_best_prot.sdf"))
    w.write(m2); w.close()
    m1, ti1 = build_correct_mol()
    m1 = map_pose_coords(m1, ROOT / "docking/results/rucaparib_docked.fixed.sdf")
    w = Chem.SDWriter(str(ROOT / "docking/results/rucaparib_docked_prot.sdf"))
    w.write(m1); w.close()
    for tag, m in (("S2", m2), ("S1", m1)):
        chg = sum(a.GetFormalCharge() for a in m.GetAtoms())
        nH = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "H")
        print(f"{tag}: {m.GetNumAtoms()} atoms, {nH} H, charge {chg}")

if __name__ == "__main__":
    main()
