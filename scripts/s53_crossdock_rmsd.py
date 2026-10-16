#!/usr/bin/env python3
"""s53_crossdock_rmsd.py — cross-docking RMSD validation (Table S16 provenance).

Purpose:  For each S1 ligand with a PARP1 CAT-domain co-crystal, compare the
          AutoDock Vina docked pose (docked into 6VKK chain A) against the
          crystallographic ligand: (1) superpose the co-crystal protein onto
          6VKK chain A (Kabsch on common CA atoms; 2RD6 uses its own numbering,
          handled by a +661 UniProt offset); (2) transform the crystal ligand
          into the docking frame; (3) compute the symmetry-aware heavy-atom
          RMSD against the pose (RDKit substructure-mapped atom pairs, minimum
          over mapped orderings).
Inputs:   data/01_modeling/docking_system1/poses/<ligand>_docked.sdf (+ .pdbqt provenance)
          data/01_modeling/docking_system1/ligands/<ligand>.sdf
          data/01_modeling/docking_system1/crystals/{6VKK,5DS3,4R6E,4UND,2RD6}.pdb
Outputs:  data/analysis/crossdock_rmsd_results.csv
Depends:  MDAnalysis, RDKit, numpy; Open Babel (obabel/obrms) on PATH
Run:      python3 scripts/s53_crossdock_rmsd.py   (from the repository root)
"""
import subprocess
import tempfile
from pathlib import Path

import numpy as np

WORK = Path(__file__).resolve().parents[1]
DOCK = WORK / "data" / "01_modeling" / "docking_system1"
OUT = WORK / "data" / "analysis" / "crossdock_rmsd_results.csv"

PAIRS = [
    # ligand, crystal, ligand resname, note
    ("rucaparib", "6VKK", "RPB", "in-frame redocking gate"),
    ("olaparib", "5DS3", "09L", ""),
    ("niraparib", "4R6E", "3JD", ""),
    ("talazoparib", "4UND", "2YQ", ""),
    ("veliparib", "2RD6", "78P", ""),
]

import MDAnalysis as mda  # noqa: E402


def kabsch(P, Q):
    """Least-squares superposition P -> Q (row-vector convention)."""
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cQ - R @ cP


def main():
    u_ref = mda.Universe(str(DOCK / "crystals" / "6VKK.pdb"))
    ref_map = {a.resid: a.position.copy() for a in u_ref.select_atoms("segid A and name CA")}

    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds

    tmpdir = Path(tempfile.mkdtemp())
    rows = []
    for lig, pdb, rn, note in PAIRS:
        u = mda.Universe(str(DOCK / "crystals" / f"{pdb}.pdb"))
        mob_map = {a.resid: a.position.copy() for a in u.select_atoms("segid A and name CA")}
        offset = 661 if pdb == "2RD6" else 0
        common = [r for r in sorted(ref_map) if (r - offset) in mob_map and 790 <= r <= 1005]
        P = np.array([mob_map[r - offset] for r in common])
        Q = np.array([ref_map[r] for r in common])
        R, t = kabsch(P, Q)
        lig_at = u.select_atoms(f"resname {rn}")
        sel = lig_at[lig_at.segids == "A"] if "A" in list(lig_at.segids) else lig_at
        frm = tmpdir / f"{pdb}_{rn}_final.pdb"
        with open(frm, "w") as f:
            for a in sel:
                x = R @ a.position + t
                f.write(f"HETATM{a.index:5d} {a.name:<4} {rn:>3} A   1    "
                        f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00          {a.name[0]:>2}\n")
            f.write("END\n")
        subprocess.run(["obabel", "-ipdb", str(frm), "-osdf", "-O",
                        str(frm).replace(".pdb", ".sdf")], capture_output=True)

        tmpl = Chem.MolFromMolFile(str(DOCK / "ligands" / f"{lig}.sdf"), removeHs=True)
        pose = Chem.MolFromMolFile(str(DOCK / "poses" / f"{lig}_docked.sdf"), removeHs=True)
        xtal = Chem.MolFromMolFile(str(frm).replace(".pdb", ".sdf"), removeHs=True)

        # talazoparib quirk: one F mis-recorded in 4UND -> compare the 26 non-F
        # heavy atoms; 6VKK's RPB extract needs bond-order repair before mapping.
        def strip_F(m):
            rw = Chem.RWMol(m)
            for i in sorted([a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() == "F"],
                            reverse=True):
                rw.RemoveAtom(i)
            m2 = rw.GetMol()
            Chem.SanitizeMol(m2)
            return m2

        if lig == "talazoparib":
            xtal_main = max(Chem.GetMolFrags(xtal, asMols=True, sanitizeFrags=True),
                            key=lambda m: m.GetNumAtoms())
            tmpl_, pose_, xtal_ = strip_F(tmpl), strip_F(pose), strip_F(xtal_main)
            rdDetermineBonds.DetermineConnectivity(xtal_, useHueckel=False)
            xtal_ = AllChem.AssignBondOrdersFromTemplate(tmpl_, xtal_)
        else:
            tmpl_, pose_, xtal_ = tmpl, pose, xtal
            rdDetermineBonds.DetermineConnectivity(xtal_, useHueckel=False)
            xtal_ = AllChem.AssignBondOrdersFromTemplate(tmpl_, xtal_)

        pm = pose_.GetSubstructMatches(tmpl_, uniquify=False)
        cm = xtal_.GetSubstructMatches(tmpl_, uniquify=False)
        if not cm:
            # fallback: direct obrms on the framed crystal ligand vs the pose
            rr = subprocess.run(["obrms", str(frm).replace(".pdb", ".sdf"),
                                 str(DOCK / "poses" / f"{lig}_docked.sdf")],
                                capture_output=True, text=True)
            val = rr.stdout.strip().split()[-1]
            best = float(val)
            print(f"{lig:<12}{pdb}/{rn}  RMSD(obrms) = {best}")
        else:
            pc, xc = pose_.GetConformer(), xtal_.GetConformer()
            best = None
            for p in pm[:80]:
                pa = np.array([list(pc.GetAtomPosition(i)) for i in p])
                for c in cm[:80]:
                    cc = np.array([list(xc.GetAtomPosition(i)) for i in c])
                    r = float(np.sqrt(((pa - cc) ** 2).sum(1).mean()))
                    best = r if best is None else min(best, r)
        rows.append((lig, pdb, rn, round(best, 2), note))
        print(f"{lig:<12}{pdb}/{rn}  RMSD = {best:.2f} A  {note}")

    with open(OUT, "w") as f:
        f.write("ligand,crystal,ligand_resname,rmsd_A,note\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    print("saved", OUT)


if __name__ == "__main__":
    main()
