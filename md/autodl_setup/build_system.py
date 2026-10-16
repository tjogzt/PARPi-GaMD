#!/usr/bin/env python3
"""
# build_system.py — GaMD system builder (protein-ligand-water)

# Flow: pdbfixer completes missing residues/atoms + adds hydrogens (pH 7.4) + removes waters -> read docked ligand (SDF, pose/bond orders kept)
# -> openmmforcefields.SystemGenerator (protein ff19SB + ligand GAFF2/OpenFF + TIP3P water)
# -> Modeller merge + solvation (padding) + 0.15 M NaCl neutralisation -> serialised system.xml + topology/coordinate PDB.

# Configured for system 1 (6VKK etc. CAT domain + inhibitor; no metal/DNA).
# System 2 (4DQY: with Zn4 + DNA): see NOTE at end; needs a Zn bonding/nonbonded model + amber DNA (OL15) force field; handled separately.

# Usage:
  python build_system.py --receptor rec.pdb --ligand pose.sdf --out md/systems/sys1_talazoparib
  python build_system.py --receptor rec.pdb --apo --out md/systems/apo   # ligand-free control
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path


def log(msg: str) -> None:
    print(f"[build] {msg}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, help="receptor PDB")
    ap.add_argument("--ligand", help="docked ligand SDF (pose kept); omit = apo")
    ap.add_argument("--apo", action="store_true", help="ligand-free system")
    ap.add_argument("--out", required=True, help="output prefix (directories auto-created)")
    ap.add_argument("--padding", type=float, default=10.0, help="water-box padding (A)")
    ap.add_argument("--ionic", type=float, default=0.15, help="ionic strength (M NaCl)")
    ap.add_argument("--ligand-ff", default="gaff-2.11",
                    help="ligand force field: gaff-2.11 | openff-2.1.0")
    ap.add_argument("--protein-ff", default="amber14-all.xml",
                    help="protein force field (amber14-all.xml = ff14SB; this openmmforcefields build lacks ff19SB)")
    ap.add_argument("--water-ff", default="amber14/tip3p.xml")
    args = ap.parse_args()

    # deferred imports (AutoDL gamd environment)
    try:
        from openmm import app, unit, XmlSerializer
        from openmm.app import Modeller, PDBFile
        from pdbfixer import PDBFixer
        from openmmforcefields.generators import SystemGenerator
    except ImportError as e:
        sys.exit(f"requires the gamd environment (openmm/pdbfixer/openmmforcefields): {e}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # 1) receptor cleaning
    log(f"PDBFixer cleaning receptor: {args.receptor}")
    fixer = PDBFixer(filename=args.receptor)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(keepWater=False)   # remove waters/buffer (metals handled separately)
    fixer.addMissingHydrogens(7.4)
    n_missing = sum(len(v) for v in fixer.missingResidues.values())
    log(f"  completing missing-residue segments: {len(fixer.missingResidues)} ({n_missing} residues), adding H at pH 7.4")

    # 2) ligand (optional)
    molecules = []
    lig_mol = None
    if not args.apo and args.ligand:
        from openff.toolkit import Molecule
        log(f"reading ligand: {args.ligand}")
        lig_mol = Molecule.from_file(args.ligand)
        if isinstance(lig_mol, list):
            lig_mol = lig_mol[0]
        # user-supplied charge injection point: if <ligand>.bcc_charges.txt exists beside it, use it directly
        # (for ligands where sqm SCF does not converge; computed offline via antechamber -ek scfconv)
        chg_file = Path(args.ligand).with_suffix(".bcc_charges.txt")
        if chg_file.exists():
            from openff.units import unit as off_unit
            vals = [float(x) for x in chg_file.read_text().split() if x.strip()]
            if len(vals) == lig_mol.n_atoms:
                lig_mol.partial_charges = off_unit.Quantity(
                    vals, off_unit.elementary_charge)
                log(f"  injected preset BCC charges ({len(vals)} atoms from {chg_file.name})")
            else:
                log(f"  [warn] charge-file atom count mismatch ({len(vals)} vs {lig_mol.n_atoms}); ignored")
        molecules.append(lig_mol)

    # 3) SystemGenerator
    log(f"SystemGenerator: protein={args.protein_ff} ligand={args.ligand_ff} water={args.water_ff}")
    ff_kwargs = {"constraints": app.HBonds, "rigidWater": True,
                 "nonbondedCutoff": 1.0 * unit.nanometer, "hydrogenMass": 1.5 * unit.amu}
    sysgen = SystemGenerator(
        forcefields=[args.protein_ff, args.water_ff],
        small_molecule_forcefield=args.ligand_ff,
        molecules=molecules if molecules else None,
        forcefield_kwargs=ff_kwargs,
        cache=str(out.parent / "ligand_ff_cache.json"),
    )

    # 4) merge + solvation
    modeller = Modeller(fixer.topology, fixer.positions)
    if lig_mol is not None:
        lig_top = lig_mol.to_topology().to_openmm()
        lig_pos = lig_mol.conformers[0].to_openmm()
        modeller.add(lig_top, lig_pos)
        log("  ligand pose merged")
    log(f"solvation: padding={args.padding} A, {args.ionic} M NaCl")
    modeller.addSolvent(
        sysgen.forcefield, model="tip3p",
        padding=args.padding * unit.angstrom,
        ionicStrength=args.ionic * unit.molar,
        neutralize=True,
    )
    n_atoms = modeller.topology.getNumAtoms()
    log(f"  system atom count: {n_atoms}")

    # 5) build System + serialise (openmm xml + amber prmtop/inpcrd for gamdRunner)
    system = sysgen.create_system(modeller.topology)
    out_xml = out.with_suffix(".system.xml")
    out_pdb = out.with_suffix(".pdb")
    out_xml.write_text(XmlSerializer.serialize(system))
    with out_pdb.open("w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh)
    # amber format (gamd-openmm gamdRunner input)
    # pitfall: the main system uses HBonds constraints + rigidWater -> H/water bonds become constraints
    #      and disappear from HarmonicBondForce; parmed then exports these bonds with type=None -> 'used' error.
    #      The AMBER prmtop must contain the full UNconstrained bond parameters; constraints (SHAKE/HBonds) are
    #      applied later by equilibrate/gamdRunner in createSystem. Hence export uses a separate constraints=None, rigidWater=False system.
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
            f"(full unconstrained bond parameters; SHAKE applied at runtime)")
    except Exception as e:
        log(f"  [warn] amber export failed ({e}); the openmm xml path can still be used")
    log(f"[OK] -> {out_xml.name}  {out_pdb.name}  (atoms={n_atoms})")
    log("next: python equilibrate.py --prmtop %s --inpcrd %s ..."
        % (out.with_suffix('.prmtop').name, out.with_suffix('.inpcrd').name))
    return 0


# NOTE system 2 (4DQY, Zn x4 + DNA duplex):
#   - DNA: add 'amber/DNA.OL15.xml' (or bsc1) to the SystemGenerator force fields
#   - Zn:  removeHeterogens strips metals; keep it and use (a) the nonbonded 12-6-4 LJ model (Li/Merz)
#          or (b) the ZAFF bonded model. Recommended: parameterise the zinc finger separately with tleap + MCPB.py, then convert to OpenMM.
#   - since PARP1+DNA+inhibitor has no experimental ternary complex, the starting conformation = 4DQY domain arrangement + docked inhibitor (see docking/),
#     with the full-length linker completed by ColabFold model superposition (see plan v3.0 section 4b).

if __name__ == "__main__":
    raise SystemExit(main())
