#!/usr/bin/env python3
"""
equilibrate.py — GaMD 前预平衡 (OpenMM, 标准 min -> NVT 升温 -> NPT)

读 amber prmtop/inpcrd -> 能量最小化 -> 限制性 NVT 逐步升温(10->300K)
-> 放松限制 NPT 密度平衡 -> 输出平衡后 rst7 (喂给 gamdRunner) + pdb。

用法:
  python equilibrate.py --prmtop sys.prmtop --inpcrd sys.inpcrd --out sys_equil \
      --temp 300 --npt-ns 1.0
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path


def log(m: str) -> None:
    print(f"[equil] {m}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--inpcrd", required=True)
    ap.add_argument("--out", required=True, help="输出前缀")
    ap.add_argument("--temp", type=float, default=300.0, help="目标温度 K")
    ap.add_argument("--nvt-ns", type=float, default=0.2, help="NVT 升温时长 ns")
    ap.add_argument("--npt-ns", type=float, default=1.0, help="NPT 平衡时长 ns")
    ap.add_argument("--restraint-k", type=float, default=10.0,
                    help="蛋白重原子位置限制 kcal/mol/Å²")
    ap.add_argument("--platform", default="CUDA")
    args = ap.parse_args()

    try:
        from openmm import (app, unit, LangevinMiddleIntegrator,
                            MonteCarloBarostat, Platform, CustomExternalForce)
        import parmed
    except ImportError as e:
        sys.exit(f"需 gamd 环境 (openmm/parmed): {e}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    K = unit.kelvin
    ps = unit.picoseconds

    log(f"读 {args.prmtop} / {args.inpcrd}")
    prm = app.AmberPrmtopFile(args.prmtop)
    crd = app.AmberInpcrdFile(args.inpcrd)
    system = prm.createSystem(nonbondedMethod=app.PME,
                              nonbondedCutoff=1.0 * unit.nanometer,
                              constraints=app.HBonds, rigidWater=True,
                              hydrogenMass=1.5 * unit.amu)

    # 蛋白重原子位置限制 (升温阶段)
    restraint = CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    restraint.addGlobalParameter("k", args.restraint_k *
                                 unit.kilocalories_per_mole / unit.angstrom**2)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    prot_res = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HID","HIE",
                "HIP","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
    n_rest = 0
    for atom in prm.topology.atoms():
        if atom.residue.name in prot_res and atom.element is not None \
                and atom.element.symbol != "H":
            pv = crd.positions[atom.index].value_in_unit(unit.nanometer)
            restraint.addParticle(atom.index, [pv.x, pv.y, pv.z])
            n_rest += 1
    restraint.setForceGroup(31)
    restraint_idx = system.addForce(restraint)
    log(f"  位置限制蛋白重原子: {n_rest}")

    # 平衡阶段用 1fs (docked 体系起步最脆弱, 小步长 + HBonds 约束最稳, 避免 NVT NaN)
    DT_PS = 0.001
    integ = LangevinMiddleIntegrator(args.temp * K, 1.0 / ps, DT_PS * ps)
    plat = Platform.getPlatformByName(args.platform)
    sim = app.Simulation(prm.topology, system, integ, plat)
    sim.context.setPositions(crd.positions)
    if crd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*crd.boxVectors)

    def pe():
        st_ = sim.context.getState(getEnergy=True)
        return st_.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    e_rest = sim.context.getState(getEnergy=True, groups={31}
                                  ).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    log(f"位置限制力初始能量 = {e_rest:.1f} kJ/mol (应≈0)")

    # 1) 最小化 (两段: 先粗解冲突, 再收敛到力阈, 压掉 docked pose 残余张力)
    log(f"最小化前 PE = {pe():.1f} kJ/mol")
    sim.minimizeEnergy(maxIterations=2000)
    log(f"  粗最小化后 PE = {pe():.1f} kJ/mol")
    sim.minimizeEnergy(tolerance=5.0 * unit.kilojoule_per_mole / unit.nanometer,
                       maxIterations=0)
    log(f"  收敛最小化后 PE = {pe():.1f} kJ/mol")

    # 2) NVT 限制性升温 10->target
    nvt_steps = int(args.nvt_ns * 1000 / DT_PS)
    log(f"NVT 升温 ({args.nvt_ns} ns, {nvt_steps} steps)...")
    n_stage = 10
    for i in range(1, n_stage + 1):
        T = 10 + (args.temp - 10) * i / n_stage
        integ.setTemperature(T * K)
        sim.context.setVelocitiesToTemperature(T * K)
        sim.step(max(1, nvt_steps // n_stage))

    # 3) 撤限制, NPT 密度平衡
    log(f"NPT 平衡 ({args.npt_ns} ns), 逐步撤位置限制...")
    barostat = MonteCarloBarostat(1.0 * unit.bar, args.temp * K, 25)
    system.addForce(barostat)
    sim.context.reinitialize(preserveState=True)
    npt_steps = int(args.npt_ns * 1000 / DT_PS)
    # 分 5 段把 k 从初值降到 0
    for i in range(5, -1, -1):
        sim.context.setParameter("k", (args.restraint_k * i / 5) *
                                 unit.kilocalories_per_mole / unit.angstrom**2)
        sim.step(max(1, npt_steps // 6))

    # 4) 输出平衡态 rst7 + pdb
    state = sim.context.getState(getPositions=True, getVelocities=True,
                                 enforcePeriodicBox=True)
    st = parmed.openmm.load_topology(prm.topology, system,
                                     xyz=state.getPositions())
    import numpy as np
    bv = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
    lengths = np.linalg.norm(bv, axis=1)
    st.box = [lengths[0], lengths[1], lengths[2], 90.0, 90.0, 90.0]
    rst7 = out.with_suffix(".rst7")
    st.save(str(rst7), format="rst7", overwrite=True)
    with out.with_suffix(".pdb").open("w") as fh:
        app.PDBFile.writeFile(prm.topology, state.getPositions(), fh)
    log(f"[OK] 平衡完成 -> {rst7.name} + {out.with_suffix('.pdb').name}")
    log("下一步: 用 gen_gamd_config.py 生成 config -> gamdRunner xml config.xml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
