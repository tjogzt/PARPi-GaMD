#!/usr/bin/env python3
"""给对接 pose 补全氢 (RDKit AddHs addCoords=True, 不动重原子)"""
from rdkit import Chem

for lig in ['AZD5305', 'veliparib']:
    src = f'docking/results_system2/{lig}_best.sdf'
    mol = Chem.MolFromMolFile(src, removeHs=False)
    # 记录重原子坐标
    conf = mol.GetConformer()
    heavy_xyz = {a.GetIdx(): list(conf.GetAtomPosition(a.GetIdx())) for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    mol_h = Chem.AddHs(mol, addCoords=True)
    # 验证重原子未动
    conf2 = mol_h.GetConformer()
    moved = 0
    for idx, xyz in heavy_xyz.items():
        p = conf2.GetAtomPosition(idx)
        d = ((p.x-xyz[0])**2 + (p.y-xyz[1])**2 + (p.z-xyz[2])**2) ** 0.5
        if d > 1e-4:
            moved += 1
    out = f'docking/results_system2/{lig}_bestH.sdf'
    Chem.MolToMolFile(mol_h, out)
    print(f'{lig}: {mol.GetNumAtoms()} → {mol_h.GetNumAtoms()} 原子 (重原子移动: {moved})')
