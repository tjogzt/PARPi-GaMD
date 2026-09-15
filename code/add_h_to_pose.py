#!/usr/bin/env python3
"""Add hydrogens to a docking pose (RDKit AddHs addCoords=True, heavy atoms untouched)"""
from rdkit import Chem

for lig in ['AZD5305', 'veliparib']:
    src = f'docking/results_system2/{lig}_best.sdf'
    mol = Chem.MolFromMolFile(src, removeHs=False)
    # record heavy-atom coordinates
    conf = mol.GetConformer()
    heavy_xyz = {a.GetIdx(): list(conf.GetAtomPosition(a.GetIdx())) for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    mol_h = Chem.AddHs(mol, addCoords=True)
    # verify heavy atoms did not move
    conf2 = mol_h.GetConformer()
    moved = 0
    for idx, xyz in heavy_xyz.items():
        p = conf2.GetAtomPosition(idx)
        d = ((p.x-xyz[0])**2 + (p.y-xyz[1])**2 + (p.z-xyz[2])**2) ** 0.5
        if d > 1e-4:
            moved += 1
    out = f'docking/results_system2/{lig}_bestH.sdf'
    Chem.MolToMolFile(mol_h, out)
    print(f'{lig}: {mol.GetNumAtoms()} -> {mol_h.GetNumAtoms()} atoms (heavy-atom shift: {moved})')
