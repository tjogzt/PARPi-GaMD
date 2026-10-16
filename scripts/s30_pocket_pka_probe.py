#!/usr/bin/env python3
#
# Purpose:  pocket-pKa probe (direction 1, quick version):
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s30_pocket_pka_probe.py   (from the repository root)

"""pocket-pKa probe (direction 1, quick version):
PROPKA3 pKa shifts of titratable residues in the S2 PARP1 pocket, apo vs
niraparib(protonated)-bound, as a first-order estimate of binding-coupled
protonation effects around the catalytic site.

Inputs:
  data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb   (S2 receptor)
  docking/results_system2/mol2/niraparib_best_prot.mol2      (protonated ligand, +1)
Outputs: results/analysis/pocket_pka_probe.txt
Note: PROPKA treats the ligand only as background electrostatics unless its
atom types are recognized; the residue-level shifts are the interpretable
output. Ligand pKa itself is NOT computed here (that needs FEP/constant-pH;
see future work).
"""
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "analysis"
OUT.mkdir(parents=True, exist_ok=True)

REC = ROOT / "data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb"
LIG = ROOT / "docking/results_system2/niraparib_best_prot.sdf"

# mol2 -> PDB (HETATM) with RDKit
from rdkit import Chem
m = Chem.MolFromMolFile(str(LIG), removeHs=False)
assert m is not None, "ligand mol2 parse failed"
conf = m.GetConformer()
lines = ["HETATM%5d  C1  UNL A 501    %8.3f%8.3f%8.3f  1.00  0.00           C  " %
         (i + 1, conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
          conf.GetAtomPosition(i).z) for i in range(m.GetNumAtoms())]
lig_pdb = "\n".join(lines) + "\nTER\nEND\n"

# apo complex: receptor only (protein+DNA+Zn) for PROPKA needs only protein;
# rewrite with MDAnalysis to guarantee standard PDB column widths (PROPKA is
# strict about fixed columns).
def write_protein_only(src, dst, include_lig=False, lig_txt=None):
    import MDAnalysis as mda
    # PARPi_full_DNA_Zn.pdb has 3 non-standard ZN HETATM lines; filter and
    # rewrite only well-formed coordinate lines to a clean PDB.
    clean = []
    for line in src.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            clean.append(line)
    tmp = dst.with_suffix(".raw.pdb")
    tmp.write_text("\n".join(clean) + "\nEND\n")
    u = mda.Universe(str(tmp))
    sel = u.select_atoms("protein or (resname ZN and name ZN)")
    sel.write(str(dst))
    tmp.unlink()
    if include_lig and lig_txt:
        with open(dst, "a") as f:
            f.write(lig_txt)

write_protein_only(REC, OUT / "pocket_apo.pdb")
write_protein_only(REC, OUT / "pocket_holo.pdb", True, lig_pdb)


def run_propka(pdb, tag):
    subprocess.run(["propka3", str(pdb), "-q"],
                   capture_output=True, text=True, timeout=600)
    cand = sorted(Path(str(pdb) + ".pka").parent.glob(f"{Path(pdb).stem}.pka"))
    return cand[0] if cand else None


apo_pka = run_propka(OUT / "pocket_apo.pdb", "apo")
holo_pka = run_propka(OUT / "pocket_holo.pdb", "holo")

report = []
if apo_pka and holo_pka:
    def parse(p):
        d = {}
        for line in Path(p).read_text().splitlines():
            if ":" in line and len(line.split()) > 4:
                toks = line.split()
                try:
                    res, val = toks[1], float(toks[2])
                    d[res] = val
                except (ValueError, IndexError):
                    continue
        return d
    apo, holo = parse(apo_pka), parse(holo_pka)
    # pocket-relevant residues (catalytic + HD interface + D766/D770 alpha-F)
    pocket = ["D766", "D770", "H862", "E988", "D770", "H826", "D740",
              "E763", "D787", "H780", "E791"]
    report.append(f"{'res':>6s} {'apo':>7s} {'holo':>7s} {'shift':>7s}")
    for r in sorted(set(pocket) & set(apo) & set(holo)):
        sh = holo[r] - apo[r]
        report.append(f"{r:>6s} {apo[r]:7.2f} {holo[r]:7.2f} {sh:+7.2f}")
    big = [l for l in report[1:] if abs(float(l.split()[-1])) >= 0.5]
    report.append("")
    report.append(f"residues with |shift| >= 0.5: {len(big)}")
    for l in big:
        report.append("  " + l)
else:
    report.append("PROPKA run incomplete; see *.pka files")

txt = "\n".join(report)
print(txt)
(OUT / "pocket_pka_probe.txt").write_text(txt + "\n")
print("WROTE", OUT / "pocket_pka_probe.txt")
