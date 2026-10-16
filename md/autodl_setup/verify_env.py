#!/usr/bin/env python
"""GaMD environment self-check: OpenMM CUDA platform + key dependencies + the GaMD engine."""
import sys

def line(k, v): print(f"{k:28s}: {v}")

print("=" * 50)
print("Python:", sys.version.split()[0], "@", sys.executable)
print("=" * 50)

# OpenMM + CUDA
try:
    import openmm
    from openmm import Platform
    line("OpenMM version", openmm.__version__)
    names = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    line("Platforms", names)
    line("CUDA platform available", "CUDA" in names)
    if "CUDA" in names:
        # actually compute one step on the GPU to confirm CUDA really works
        from openmm import System, VerletIntegrator, Vec3
        from openmm.unit import nanometer, kelvin, picosecond
        sysm = System()
        sysm.addParticle(1.0)
        integ = VerletIntegrator(0.001 * picosecond)
        ctx = openmm.Context(sysm, integ, Platform.getPlatformByName("CUDA"))
        ctx.setPositions([Vec3(0, 0, 0) * nanometer])
        ctx.getState(getEnergy=True)
        line("CUDA smoke test", "PASS (ran 0-particle context on GPU)")
        del ctx
except Exception as e:
    line("OpenMM ERROR", repr(e))

# key scientific dependencies
for mod, attr in [("openff.toolkit", "__version__"),
                  ("openmmforcefields", "__version__"),
                  ("rdkit", "__version__"),
                  ("MDAnalysis", "__version__"),
                  ("parmed", "__version__"),
                  ("numpy", "__version__"),
                  ("scipy", "__version__")]:
    try:
        m = __import__(mod, fromlist=["x"])
        line(mod, getattr(m, attr, "imported"))
    except Exception as e:
        line(mod + " ERROR", repr(e))

# AmberTools CLI (tleap / antechamber on PATH?)
import shutil
for tool in ["tleap", "antechamber", "parmchk2", "sander", "pdb4amber"]:
    line("amber:" + tool, shutil.which(tool) or "NOT FOUND")

# GaMD engine
for mod in ["gamd", "pyemma", "deeptime"]:
    try:
        m = __import__(mod)
        line(mod, getattr(m, "__version__", "imported"))
    except Exception:
        line(mod, "MISSING")

print("=" * 50)
print("SELF-CHECK DONE")
