#!/usr/bin/env python3
"""Dihedral-group energy via OpenMM (exact gamdRunner definition).

The gamdRunner "dihedral" boost group = PeriodicTorsionForce + CMAPTorsionForce
(tools/gamd-openmm/gamd/integrator_factory.py:40). This script reproduces that
energy per frame from a stripped prmtop + trajectory, for the P0-2 S2 DBE
reconstruction of systems whose gamd.log was not archived.

Validation: on APO/talazoparib the OpenMM dihedral-group energy must match the
archived gamd.log column 4 (Unboosted-Dihedral-Energy) at the aligned frames.

Usage:
  python3 scripts/dihedral_group_energy.py <prmtop> <trajectory> <out.csv> [stride]

Trajectory formats: anything MDAnalysis reads (DCD/NetCDF).
"""
import sys
from pathlib import Path

import numpy as np

try:
    from openmm import app
    import openmm
    import MDAnalysis as mda
except ImportError as e:
    sys.exit(f"missing dependency: {e}")

FF_XMLS = [
    "amber14/protein.ff14SB.xml",
    "amber14/DNA.bsc1.xml",
    "amber14/tip3p.xml",
]


def build_system(prmtop):
    prm = app.AmberPrmtopFile(prmtop)
    system = prm.createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        rigidWater=False,
    )
    # dihedral group = PeriodicTorsionForce + CMAPTorsionForce (group 2)
    for force in system.getForces():
        if isinstance(force, (openmm.PeriodicTorsionForce,
                              openmm.CMAPTorsionForce)):
            force.setForceGroup(2)
    return prm, system


def main():
    prmtop = sys.argv[1]
    traj = sys.argv[2]
    out_csv = sys.argv[3]
    stride = int(sys.argv[4]) if len(sys.argv) > 4 else 1

    u = mda.Universe(prmtop, traj)
    prm, system = build_system(prmtop)
    integrator = openmm.VerletIntegrator(1.0)
    context = openmm.Context(system, integrator)
    n = len(u.trajectory)
    with open(out_csv, "w") as fh:
        fh.write("frame,index,E_dihedral_kcal\n")
        for i, ts in enumerate(u.trajectory[::stride]):
            context.setPositions(ts.positions)
            e = context.getState(getEnergy=True, groups={2}).getPotentialEnergy()
            e_kcal = e.value_in_unit(openmm.unit.kilocalorie_per_mole)
            fh.write(f"{i},{i * stride},{e_kcal:.6f}\n")
            if (i + 1) % 100 == 0:
                print(f"frame {i+1}/{n // stride}", flush=True)
    print(f"done: {n // stride} frames -> {out_csv}")


if __name__ == "__main__":
    main()
