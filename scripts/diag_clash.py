#!/usr/bin/env python3
"""
diag_clash.py — 定位平衡前体系的灾难性高能来源 (PE ~1e8 kJ/mol)。

1) 按 force group 拆分单点能量 -> 区分 bond/angle/torsion/nonbonded 哪个爆。
2) cKDTree 扫最近原子对 (<1.0 Å) -> 锁定重叠原子, 标注 residue/atom/是否成键。
3) 报告配体残基的近邻接触, 判断 docked pose 是否撞蛋白。

用法: python diag_clash.py sys.prmtop sys.inpcrd
"""
import sys
import numpy as np
from openmm import app, unit, Platform


def main():
    prmtop, inpcrd = sys.argv[1], sys.argv[2]
    print(f"[diag] 读 {prmtop} / {inpcrd}", flush=True)
    prm = app.AmberPrmtopFile(prmtop)
    crd = app.AmberInpcrdFile(inpcrd)
    # 无约束, 全键参数保留, 用于真实单点能
    system = prm.createSystem(nonbondedMethod=app.PME,
                              nonbondedCutoff=1.0 * unit.nanometer,
                              constraints=None, rigidWater=False)

    # ---- 1) 按 force group 拆能量 ----
    forces = system.getForces()
    for i, f in enumerate(forces):
        f.setForceGroup(i)
    from openmm import LangevinIntegrator
    integ = LangevinIntegrator(300, 1.0, 0.001)
    plat = Platform.getPlatformByName("CPU")
    ctx = app.Simulation(prm.topology, system, integ, plat).context
    pos = crd.positions
    ctx.setPositions(pos)
    if crd.boxVectors is not None:
        ctx.setPeriodicBoxVectors(*crd.boxVectors)
    tot = ctx.getState(getEnergy=True).getPotentialEnergy()
    print(f"[diag] 总 PE = {tot.value_in_unit(unit.kilojoule_per_mole):.3e} kJ/mol")
    print("[diag] 按力分组能量:")
    for i, f in enumerate(forces):
        e = ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy()
        print(f"    {f.__class__.__name__:28s} = "
              f"{e.value_in_unit(unit.kilojoule_per_mole):.4e} kJ/mol")

    # ---- 2) 最近原子对扫描 ----
    xyz = np.array(crd.positions.value_in_unit(unit.nanometer))  # (N,3) nm
    atoms = list(prm.topology.atoms())
    # 成键对集合 (用于标注)
    bonded = set()
    for b in prm.topology.bonds():
        i, j = b[0].index, b[1].index
        bonded.add((min(i, j), max(i, j)))

    print(f"\n[diag] 原子数 {len(atoms)}; 扫描 <1.0 Å 接触对 ...", flush=True)
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(xyz)
        pairs = tree.query_pairs(r=0.10)  # 0.10 nm = 1.0 Å
    except Exception as e:
        print(f"  [warn] scipy 不可用 ({e}), 跳过对扫描")
        pairs = set()

    rows = []
    for i, j in pairs:
        d = float(np.linalg.norm(xyz[i] - xyz[j])) * 10.0  # Å
        is_bond = (min(i, j), max(i, j)) in bonded
        rows.append((d, i, j, is_bond))
    rows.sort()

    def desc(idx):
        a = atoms[idx]
        return f"{a.residue.name}{a.residue.id}/{a.name}#{idx}"

    nonbond_clash = [r for r in rows if not r[3]]
    print(f"[diag] 非成键的 <1.0 Å 重叠对: {len(nonbond_clash)} (这些就是病灶)")
    print("[diag] 最近 30 对 (d=距离Å, B=是否成键):")
    for d, i, j, isb in rows[:30]:
        flag = "BOND" if isb else "CLASH"
        print(f"    {d:6.3f} Å  [{flag:5s}]  {desc(i):28s} -- {desc(j)}")

    # ---- 3) 配体残基近邻 ----
    res_names = {}
    for a in atoms:
        res_names.setdefault(a.residue.name, 0)
        res_names[a.residue.name] += 1
    # 非标准残基 (配体/离子/水之外)
    std = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HID","HIE",
           "HIP","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR",
           "VAL","HOH","WAT","NA","CL","Na+","Cl-"}
    lig = [n for n in res_names if n not in std]
    print(f"\n[diag] 非标准残基(配体候选): {lig}")
    for ln in lig:
        idxs = [a.index for a in atoms if a.residue.name == ln]
        if not idxs:
            continue
        lo, hi = min(idxs), max(idxs)
        print(f"  {ln}: {len(idxs)} 原子, index {lo}-{hi}")


if __name__ == "__main__":
    main()
