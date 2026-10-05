#!/usr/bin/env python3
"""
p04_equilibrate.py — pre-GaMD equilibration (OpenMM: minimisation -> NVT -> NPT)

Purpose:   Read an Amber prmtop/inpcrd -> energy minimisation -> restrained NVT
           heating (10 -> 300 K) -> unrestrained NPT density equilibration ->
           write the equilibrated rst7 (fed to gamdRunner) + pdb.
Author:    Tao Zhu (tjogzt@gmail.com)
Created:   2026-09-15
Inputs:    --prmtop/--inpcrd  Amber topology + coordinates (from p01)
Outputs:   <out>.rst7 + <out>.pdb (equilibrated state)
Depends:   openmm, parmed (conda env 'gamd')
Run:       python3 code/pipeline/p04_equilibrate.py --prmtop sys.prmtop \
               --inpcrd sys.inpcrd --out sys_equil --temp 300 --npt-ns 1.0
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
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--temp", type=float, default=300.0, help="target temperature (K)")
    ap.add_argument("--nvt-ns", type=float, default=0.2, help="NVT heating duration (ns)")
    ap.add_argument("--npt-ns", type=float, default=1.0, help="NPT equilibration duration (ns)")
    ap.add_argument("--restraint-k", type=float, default=10.0,
                    help="protein heavy-atom positional restraint (kcal/mol/A^2)")
    ap.add_argument("--platform", default="CUDA")
    args = ap.parse_args()

    try:
        from openmm import (app, unit, LangevinMiddleIntegrator,
                            MonteCarloBarostat, Platform, CustomExternalForce)
        import parmed
    except ImportError as e:
        sys.exit(f"gamd environment required (openmm/parmed): {e}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    K = unit.kelvin
    ps = unit.picoseconds

    log(f"reading {args.prmtop} / {args.inpcrd}")
    prm = app.AmberPrmtopFile(args.prmtop)
    crd = app.AmberInpcrdFile(args.inpcrd)
    system = prm.createSystem(nonbondedMethod=app.PME,
                              nonbondedCutoff=1.0 * unit.nanometer,
                              constraints=app.HBonds, rigidWater=True,
                              hydrogenMass=1.5 * unit.amu)

    # Protein heavy-atom positional restraint (heating stage)
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
    log(f"  restrained protein heavy atoms: {n_rest}")

    # 1 fs during equilibration (docked systems are most fragile at the start;
    # small steps + HBonds constraints are the most stable, avoiding NVT NaNs)
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
    log(f"restraint force initial energy = {e_rest:.1f} kJ/mol (should be ~0)")

    # 1) Minimisation (two stages: coarse clash relief, then converged to the
    # force tolerance, releasing residual docked-pose strain)
    log(f"PE before minimisation = {pe():.1f} kJ/mol")
    sim.minimizeEnergy(maxIterations=2000)
    log(f"  PE after coarse minimisation = {pe():.1f} kJ/mol")
    sim.minimizeEnergy(tolerance=5.0 * unit.kilojoule_per_mole / unit.nanometer,
                       maxIterations=0)
    log(f"  PE after converged minimisation = {pe():.1f} kJ/mol")

    # 2) Restrained NVT heating 10 -> target
    nvt_steps = int(args.nvt_ns * 1000 / DT_PS)
    log(f"NVT heating ({args.nvt_ns} ns, {nvt_steps} steps)...")
    n_stage = 10
    for i in range(1, n_stage + 1):
        T = 10 + (args.temp - 10) * i / n_stage
        integ.setTemperature(T * K)
        sim.context.setVelocitiesToTemperature(T * K)
        sim.step(max(1, nvt_steps // n_stage))

    # 3) Remove restraints, NPT density equilibration
    log(f"NPT equilibration ({args.npt_ns} ns), releasing restraints stepwise...")
    barostat = MonteCarloBarostat(1.0 * unit.bar, args.temp * K, 25)
    system.addForce(barostat)
    sim.context.reinitialize(preserveState=True)
    npt_steps = int(args.npt_ns * 1000 / DT_PS)
    # Ramp k from its initial value down to 0 in five stages
    for i in range(5, -1, -1):
        sim.context.setParameter("k", (args.restraint_k * i / 5) *
                                 unit.kilocalories_per_mole / unit.angstrom**2)
        sim.step(max(1, npt_steps // 6))

    # 4) Write the equilibrated rst7 + pdb
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
    log(f"[OK] equilibration complete -> {rst7.name} + {out.with_suffix('.pdb').name}")
    log("Next: p05_gen_gamd_config.py to build the config -> gamdRunner xml config.xml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
