#!/usr/bin/env python3
"""RDKit → mol2 转换器 (RDKit 无 mol2 输出器, 手写 TRIPOS 格式)
键级: 芳香键统一 ar; 原子类型简单映射(antechamber 会用 -at gaff2 重赋)
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
    # AZD5305/veliparib 的对接 SDF 缺氢 (老对接), 先用 add_h_to_pose.py 补全 → _bestH.sdf
    for lig in ['fluzoparib', 'pamiparib', 'senaparib', 'AZD5305', 'veliparib']:
        sdf = f'docking/results_system2/{lig}_best.sdf'
        sdf_h = f'docking/results_system2/{lig}_bestH.sdf'
        if os.path.exists(sdf_h):
            sdf = sdf_h
        mol = Chem.MolFromMolFile(sdf, removeHs=False)
        Chem.SanitizeMol(mol)
        # 键级修复:
        # 1) 环内芳香-芳香键 → ar
        # 2) 芳香原子与非芳香原子之间的 DOUBLE → SINGLE (Kekulé 交替双键落在环外键的假象, 如芳环-羰基)
        for b in mol.GetBonds():
            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
            if a1.GetIsAromatic() and a2.GetIsAromatic():
                b.SetBondType(Chem.BondType.AROMATIC)
            elif (a1.GetIsAromatic() or a2.GetIsAromatic()) and b.GetBondType() == Chem.BondType.DOUBLE:
                b.SetBondType(Chem.BondType.SINGLE)
        out = f'docking/results_system2/mol2/{lig}_best.mol2'
        rdkit_to_mol2(mol, lig, out)
        # 自检
        mixed = [b for b in mol.GetBonds() if b.GetBeginAtom().GetIsAromatic()
                 and b.GetEndAtom().GetIsAromatic() and b.GetBondType() != Chem.BondType.AROMATIC]
        print(f'{lig}: {mol.GetNumAtoms()} 原子 {mol.GetNumBonds()} 键 → {out} (芳香混合键残留: {len(mixed)})')
