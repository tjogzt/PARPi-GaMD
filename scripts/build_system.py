#!/usr/bin/env python3
"""
build_system.py — GaMD 体系搭建 (蛋白-配体-水)

流程: pdbfixer 补缺失残基/原子+加氢(pH7.4)+去水 -> 读对接配体(SDF,保留 pose/键级)
-> openmmforcefields.SystemGenerator (蛋白 ff19SB + 配体 GAFF2/OpenFF + 水 TIP3P)
-> Modeller 合并+溶剂化(padding)+0.15M NaCl 中和 -> 序列化 system.xml + 拓扑/坐标 pdb。

适配体系1 (6VKK 等 CAT 域+抑制剂, 无金属/DNA)。
体系2 (4DQY: 含 Zn4 + DNA): 见文末 NOTE, 需 Zn 键合/非键模型 + amber DNA(OL15) 力场, 单独处理。

用法:
  python build_system.py --receptor rec.pdb --ligand pose.sdf --out md/systems/sys1_talazoparib
  python build_system.py --receptor rec.pdb --apo --out md/systems/apo   # 无配体对照
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path


def log(msg: str) -> None:
    print(f"[build] {msg}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, help="受体 PDB")
    ap.add_argument("--ligand", help="对接配体 SDF (保留 pose); 省略=apo")
    ap.add_argument("--apo", action="store_true", help="无配体体系")
    ap.add_argument("--out", required=True, help="输出前缀 (目录会自动建)")
    ap.add_argument("--padding", type=float, default=10.0, help="水盒 padding (Å)")
    ap.add_argument("--ionic", type=float, default=0.15, help="离子强度 (M NaCl)")
    ap.add_argument("--ligand-ff", default="gaff-2.11",
                    help="配体力场: gaff-2.11 | openff-2.1.0")
    ap.add_argument("--protein-ff", default="amber14-all.xml",
                    help="蛋白力场 (amber14-all.xml=ff14SB; openmmforcefields 此版本未带 ff19SB)")
    ap.add_argument("--water-ff", default="amber14/tip3p.xml")
    args = ap.parse_args()

    # 延迟导入 (AutoDL gamd 环境)
    try:
        from openmm import app, unit, XmlSerializer
        from openmm.app import Modeller, PDBFile
        from pdbfixer import PDBFixer
        from openmmforcefields.generators import SystemGenerator
    except ImportError as e:
        sys.exit(f"需在 gamd 环境运行 (openmm/pdbfixer/openmmforcefields): {e}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # 1) 受体清洗
    log(f"PDBFixer 清洗受体: {args.receptor}")
    fixer = PDBFixer(filename=args.receptor)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(keepWater=False)   # 去水/缓冲剂 (金属另行处理)
    fixer.addMissingHydrogens(7.4)
    n_missing = sum(len(v) for v in fixer.missingResidues.values())
    log(f"  补缺失残基段: {len(fixer.missingResidues)} (共 {n_missing} 残基), 加氢 pH7.4")

    # 2) 配体 (可选)
    molecules = []
    lig_mol = None
    if not args.apo and args.ligand:
        from openff.toolkit import Molecule
        log(f"读配体: {args.ligand}")
        lig_mol = Molecule.from_file(args.ligand)
        if isinstance(lig_mol, list):
            lig_mol = lig_mol[0]
        molecules.append(lig_mol)

    # 3) SystemGenerator
    log(f"SystemGenerator: 蛋白={args.protein_ff} 配体={args.ligand_ff} 水={args.water_ff}")
    ff_kwargs = {"constraints": app.HBonds, "rigidWater": True,
                 "nonbondedCutoff": 1.0 * unit.nanometer, "hydrogenMass": 1.5 * unit.amu}
    sysgen = SystemGenerator(
        forcefields=[args.protein_ff, args.water_ff],
        small_molecule_forcefield=args.ligand_ff,
        molecules=molecules if molecules else None,
        forcefield_kwargs=ff_kwargs,
        cache=str(out.parent / "ligand_ff_cache.json"),
    )

    # 4) 合并 + 溶剂化
    modeller = Modeller(fixer.topology, fixer.positions)
    if lig_mol is not None:
        lig_top = lig_mol.to_topology().to_openmm()
        lig_pos = lig_mol.conformers[0].to_openmm()
        modeller.add(lig_top, lig_pos)
        log("  已并入配体 pose")
    log(f"溶剂化: padding={args.padding}Å, {args.ionic}M NaCl")
    modeller.addSolvent(
        sysgen.forcefield, model="tip3p",
        padding=args.padding * unit.angstrom,
        ionicStrength=args.ionic * unit.molar,
        neutralize=True,
    )
    n_atoms = modeller.topology.getNumAtoms()
    log(f"  体系原子数: {n_atoms}")

    # 5) 构建 System + 序列化 (openmm xml + amber prmtop/inpcrd 供 gamdRunner)
    system = sysgen.create_system(modeller.topology)
    out_xml = out.with_suffix(".system.xml")
    out_pdb = out.with_suffix(".pdb")
    out_xml.write_text(XmlSerializer.serialize(system))
    with out_pdb.open("w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh)
    # amber 格式 (gamd-openmm gamdRunner 输入)
    # 坑: 主 system 用 HBonds 约束 + rigidWater -> 含氢/水键被转成 constraint,
    #     从 HarmonicBondForce 消失, parmed 导 AMBER 时这些键 type=None -> 'used' 报错。
    # AMBER prmtop 须含"无约束"完整键参数; 约束(SHAKE/HBonds) 由 equilibrate/gamdRunner
    # 在 createSystem 时再施加。故为导出单独建 constraints=None, rigidWater=False 的 system。
    try:
        import parmed
        export_system = sysgen.forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=None, rigidWater=False, removeCMMotion=True)
        st = parmed.openmm.load_topology(modeller.topology, export_system,
                                         xyz=modeller.positions)
        st.save(str(out.with_suffix(".prmtop")), overwrite=True)
        st.save(str(out.with_suffix(".inpcrd")), format="rst7", overwrite=True)
        log(f"  amber: {out.with_suffix('.prmtop').name} + .inpcrd "
            f"(无约束完整键参数; SHAKE 运行时施加)")
    except Exception as e:
        log(f"  [warn] amber 导出失败 ({e}); 仍可用 openmm xml 路径")
    log(f"[OK] -> {out_xml.name}  {out_pdb.name}  (atoms={n_atoms})")
    log("下一步: python equilibrate.py --prmtop %s --inpcrd %s ..."
        % (out.with_suffix('.prmtop').name, out.with_suffix('.inpcrd').name))
    return 0


# NOTE 体系2 (4DQY, 含 Zn x4 + DNA 双链):
#   - DNA: SystemGenerator forcefields 加 'amber/DNA.OL15.xml' (或 bsc1)
#   - Zn:  removeHeterogens 会去掉金属; 需保留并用 (a) 非键 12-6-4 LJ 模型 (Li/Merz)
#          或 (b) ZAFF 键合模型. 建议 tleap+MCPB.py 单独参数化锌指, 再转 OpenMM.
#   - 因 PARP1+DNA+抑制剂 无实验三元复合物, 起始构象 = 4DQY 域排布 + 对接抑制剂 (见 docking/),
#     全长 linker 用 ColabFold 模型叠合补全 (见 plan v3.0 §4b).

if __name__ == "__main__":
    raise SystemExit(main())
