#!/usr/bin/env python3
"""
make_docked_sdf.py — 由 Vina docked pdbqt 生成化学正确的 docked-pose SDF

pdbqt 丢失键级且并掉非极性氢, 直接喂 openff 会参数化错误。
策略(坐标迁移, 最稳健): 化学完全取自"原始配体 SDF"(键级/电荷/质子化正确),
仅把 docked pose 的重原子坐标迁移过去, 再几何加氢 -> 分子式天然正确。
  1) obabel 取 docked 第1个 pose -> 仅用其重原子坐标
  2) 用"骨架(全单键去芳香)"做元素+连接性匹配, 得 ref<->docked 原子映射
  3) 把 docked 坐标写到 ref_heavy 的新 conformer
  4) AddHs(addCoords) 按几何补氢 (保留 docked 手性)
  5) 分子式自检: 必须 == 模板加氢后

用法:
  python make_docked_sdf.py --ref docking/ligands/talazoparib.sdf \
      --pdbqt docking/results/talazoparib_docked.pdbqt \
      --out  docking/results/talazoparib_docked.fixed.sdf
"""
from __future__ import annotations
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def skeleton(m):
    """返回全单键/去芳香/去电荷的骨架 copy, 仅用于元素+连接性子结构匹配。"""
    from rdkit import Chem
    em = Chem.RWMol(m)
    for b in em.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in em.GetAtoms():
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
        a.SetFormalCharge(0)
        a.SetNumExplicitHs(0)
    mm = em.GetMol()
    Chem.SanitizeMol(mm, catchErrors=True)
    return mm


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="原始配体 SDF (键级模板, 化学来源)")
    ap.add_argument("--pdbqt", required=True, help="Vina docked pdbqt (坐标来源)")
    ap.add_argument("--out", required=True, help="输出 docked SDF")
    ap.add_argument("--pose", type=int, default=1, help="取第几个 pose (默认1=best)")
    args = ap.parse_args()

    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from rdkit.Geometry import Point3D
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.warning")

    # 1) obabel: docked pdbqt -> 单 pose sdf (只取坐标)
    tmp = Path(tempfile.mkdtemp()) / "pose.sdf"
    r = subprocess.run(
        ["obabel", args.pdbqt, "-O", str(tmp), "-f", str(args.pose), "-l", str(args.pose)],
        capture_output=True, text=True)
    if not tmp.exists() or tmp.stat().st_size == 0:
        sys.exit(f"obabel 转换失败: {r.stderr.strip()[:300]}")

    # 2) 参考 (化学正确, 重原子) + docked (坐标, 键级不可靠)
    ref_heavy = Chem.MolFromMolFile(args.ref, removeHs=True)
    if ref_heavy is None:
        sys.exit(f"无法读取参考 SDF: {args.ref}")
    docked = Chem.MolFromMolFile(str(tmp), removeHs=False, sanitize=False)
    if docked is None:
        sys.exit("无法读取 obabel 转换的 docked sdf")
    docked.UpdatePropertyCache(strict=False)
    # sanitize=False 时 removeHs 不可靠 -> 手动剔除 H (原子序数=1), RWMol 会同步 conformer 坐标
    rw = Chem.RWMol(docked)
    for idx in sorted([a.GetIdx() for a in docked.GetAtoms() if a.GetAtomicNum() == 1], reverse=True):
        rw.RemoveAtom(idx)
    docked = rw.GetMol()
    docked.UpdatePropertyCache(strict=False)

    if ref_heavy.GetNumAtoms() != docked.GetNumAtoms():
        sys.exit(f"重原子数不一致: ref={ref_heavy.GetNumAtoms()} docked={docked.GetNumAtoms()}")

    # 3) 骨架匹配 -> ref 原子 i 对应 docked 原子 match[i]
    match = skeleton(docked).GetSubstructMatch(skeleton(ref_heavy))
    if len(match) != ref_heavy.GetNumAtoms():
        sys.exit("骨架子结构匹配失败 (ref 与 docked 重原子连接性不一致?)")

    # 4) 把 docked 坐标迁移到 ref_heavy 的新 conformer
    dconf = docked.GetConformer()
    newconf = Chem.Conformer(ref_heavy.GetNumAtoms())
    for i in range(ref_heavy.GetNumAtoms()):
        p = dconf.GetAtomPosition(match[i])
        newconf.SetAtomPosition(i, Point3D(p.x, p.y, p.z))
    ref_heavy.RemoveAllConformers()
    ref_heavy.AddConformer(newconf, assignId=True)

    # 5) 几何加氢 (保留 docked 重原子坐标与手性)
    fixed = Chem.AddHs(ref_heavy, addCoords=True)
    Chem.SanitizeMol(fixed)

    # 6) 分子式自检: 必须与"模板加氢后"完全一致
    f_fixed = rdMolDescriptors.CalcMolFormula(fixed)
    f_templ = rdMolDescriptors.CalcMolFormula(Chem.AddHs(Chem.MolFromMolFile(args.ref, removeHs=True)))
    if f_fixed != f_templ:
        sys.exit(f"[FAIL] 分子式不符: docked={f_fixed} vs 模板={f_templ}")

    # 7) 写 SDF
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(args.out)
    w.write(fixed)
    w.close()

    print(f"[OK] {args.out}  重原子={ref_heavy.GetNumAtoms()} 总原子={fixed.GetNumAtoms()} 分子式={f_fixed}")
    print(f"     SMILES(docked) = {Chem.MolToSmiles(Chem.RemoveHs(fixed))}")
    print(f"     SMILES(模板)   = {Chem.MolToSmiles(Chem.MolFromMolFile(args.ref, removeHs=True))}")
    print(f"     分子式自检通过: {f_fixed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
