#!/usr/bin/env python3
"""p06_rdkit_to_mol2.py — convert docked RDKit poses to TRIPOS mol2

Purpose:   RDKit has no mol2 writer; this script hand-writes the TRIPOS
           format. Bond orders: intra-ring aromatic bonds are forced to 'ar';
           atom types are mapped simply (antechamber reassigns them with
           '-at gaff2' later). Aromatic-aliphatic DOUBLE bonds left over from
           Kekule input are demoted to SINGLE.
Created:   2026-09-15
Inputs:    data/01_modeling/docking_system2/<lig>_best.sdf (or _bestH.sdf)
Outputs:   data/01_modeling/docking_system2/mol2/<lig>_best.mol2
Depends:   rdkit
Run:       python3 code/pipeline/p06_rdkit_to_mol2.py   (from the repo root)
"""
import os
from rdkit import Chem

ATOMTYPE = {
    'C': {Chem.HybridizationType.SP3: 'C.3', Chem.HybridizationType.SP2: 'C.2', Chem.HybridizationType.SP: 'C.1'},
    'N': {Chem.HybridizationType.SP3: 'N.3', Chem.HybridizationType.SP2: 'N.2', Chem.HybridizationType.SP: 'N.1'},
    'O': {Chem.HybridizationType.SP3: 'O.3', Chem.HybridizationType.SP2: 'O.2'},
    'S': {Chem.HybridizationType.SP3: 'S.3', Chem.HybridizationType.SP2: 'S.2'},
    'F': 'F', 'Cl': 'Cl', 'Br': 'Br', 'I': 'I', 'P': 'P.3', 'H': 'H', 'B': 'B.2',
}
BONDTYPE = {Chem.BondType.AROMATIC: 'ar', Chem.BondType.SINGLE: '1',
            Chem.BondType.DOUBLE: '2', Chem.BondType.TRIPLE: '3'}


def rdkit_to_mol2(mol, name, out):
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()
    lines = [f'@<TRIPOS>MOLECULE', name, f'  {n_atoms} {n_bonds} 0 0 0', 'SMALL', 'USER_CHARGES', '',
             '@<TRIPOS>ATOM']
    conf = mol.GetConformer()
    for a in mol.GetAtoms():
        i = a.GetIdx()
        sym = a.GetSymbol()
        pos = conf.GetAtomPosition(i)
        if a.GetIsAromatic():
            t = {'C': 'C.ar', 'N': 'N.ar', 'S': 'S.ar', 'O': 'O.ar'}.get(sym, sym)
        else:
            t = ATOMTYPE_MAP(sym, a)
        lines.append(f'{i+1:>6d} {sym:<4s} {pos.x:10.4f} {pos.y:10.4f} {pos.z:10.4f} '
                     f'{t:<6s} 1  UNL1        0.0000')
    lines.append('@<TRIPOS>BOND')
    for b in mol.GetBonds():
        bt = BONDTYPE.get(b.GetBondType(), '1')
        lines.append(f'{b.GetIdx()+1:>6d} {b.GetBeginAtomIdx()+1:5d} {b.GetEndAtomIdx()+1:5d} {bt}')
    with open(out, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def ATOMTYPE_MAP(sym, atom):
    d = ATOMTYPE.get(sym, {})
    if isinstance(d, dict):
        return d.get(atom.GetHybridization(), sym)
    return d


if __name__ == '__main__':
    # The AZD5305/veliparib docked SDFs lack hydrogens (older docking); run
    # code/47_add_hydrogens_to_pose.py first to produce _bestH.sdf
    for lig in ['fluzoparib', 'pamiparib', 'senaparib', 'AZD5305', 'veliparib']:
        sdf = f'data/01_modeling/docking_system2/{lig}_best.sdf'
        sdf_h = f'data/01_modeling/docking_system2/{lig}_bestH.sdf'
        if os.path.exists(sdf_h):
            sdf = sdf_h
        mol = Chem.MolFromMolFile(sdf, removeHs=False)
        Chem.SanitizeMol(mol)
        # Bond-order fixes:
        # 1) intra-ring aromatic-aromatic bonds -> 'ar'
        # 2) DOUBLE bonds between an aromatic and a non-aromatic atom -> SINGLE
        #    (a Kekule alternation artefact on exocyclic bonds, e.g. aryl-carbonyl)
        for b in mol.GetBonds():
            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
            if a1.GetIsAromatic() and a2.GetIsAromatic():
                b.SetBondType(Chem.BondType.AROMATIC)
            elif (a1.GetIsAromatic() or a2.GetIsAromatic()) and b.GetBondType() == Chem.BondType.DOUBLE:
                b.SetBondType(Chem.BondType.SINGLE)
        out = f'data/01_modeling/docking_system2/mol2/{lig}_best.mol2'
        rdkit_to_mol2(mol, lig, out)
        # self-check
        mixed = [b for b in mol.GetBonds() if b.GetBeginAtom().GetIsAromatic()
                 and b.GetEndAtom().GetIsAromatic() and b.GetBondType() != Chem.BondType.AROMATIC]
        print(f'{lig}: {mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds -> {out} '
          f'(leftover mixed aromatic bonds: {len(mixed)})')
