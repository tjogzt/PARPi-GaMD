#!/usr/bin/env python3
#
# Purpose:  GaMD starting-structure diagnostics (pre-cleaning)
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   (see docstring)
# Outputs:  data/00_raw/pdb/structure_report.md
# Depends:  Bio
# Run:      python3 scripts/s13_analyze_structures.py

"""
analyze_structures.py -- GaMD starting-structure diagnostics (pre-cleaning)

Parses each PDB: chain composition (protein/nucleic-acid/residue ranges), non-water ligands, missing residues (REMARK 465).
Output: data/00_raw/pdb/structure_report.md
Purpose: decide the system-construction strategy + whether ColabFold loop completion is needed; feeds the Methods report.
"""
from __future__ import annotations
import re
from pathlib import Path
from collections import defaultdict, Counter

from Bio.PDB import PDBParser, is_aa
from Bio.PDB.Polypeptide import PPBuilder

PDB_DIR = Path(__file__).resolve().parents[1] / "data" / "00_raw" / "pdb"
OUT = PDB_DIR / "structure_report.md"
PDBS = ["4DQY", "6VKK", "6VKO", "6VKQ", "6NTU", "6BHV"]

NUCLEIC = {"DA", "DT", "DG", "DC", "DU", "A", "U", "G", "C"}
IGNORE_HET = {"HOH", "WAT"}
COMMON_IONS = {"ZN", "MG", "NA", "CL", "K", "CA", "MN", "SO4", "PO4",
               "GOL", "EDO", "ACT", "PEG", "DMS", "FMT", "CO3"}


def parse_remark465(pdb_path: Path) -> dict[str, int]:
    """Count missing residues per chain (REMARK 465)."""
    missing = defaultdict(int)
    in_block = False
    for line in pdb_path.read_text(errors="ignore").splitlines():
        if line.startswith("REMARK 465"):
            body = line[10:].strip()
            # skip header rows
            m = re.match(r"^([A-Z0-9]{1,3})\s+([A-Za-z0-9])\s+(-?\d+)", body)
            if m and "RES" not in body and "SSSEQI" not in body:
                missing[m.group(2)] += 1
        elif in_block and not line.startswith("REMARK"):
            break
    return dict(missing)


def classify_chain(chain) -> str:
    aa = na = 0
    for res in chain:
        rn = res.get_resname().strip()
        if is_aa(res, standard=False):
            aa += 1
        elif rn in NUCLEIC:
            na += 1
    if aa and aa >= na:
        return "protein"
    if na:
        return "nucleic"
    return "other"


def main() -> None:
    parser = PDBParser(QUIET=True)
    lines = ["# GaMD starting-structure diagnostics report\n",
             "_Auto-generated (Bio.PDB). Ligands including metals/buffer classified. Missing residues from REMARK 465._\n"]
    for pid in PDBS:
        path = PDB_DIR / f"{pid}.pdb"
        if not path.exists():
            lines.append(f"\n## {pid}\n- **file missing**\n")
            continue
        s = parser.get_structure(pid, str(path))
        model = next(iter(s))
        missing = parse_remark465(path)
        lines.append(f"\n## {pid}\n")
        # chains
        lines.append("| Chain | Type | Residues | Range | Missing (465) |")
        lines.append("|----|------|:------:|:--------:|:---------:|")
        het_counter = Counter()
        for chain in model:
            ctype = classify_chain(chain)
            poly = [r for r in chain if r.id[0] == " "]
            het = [r for r in chain if r.id[0] not in (" ", "W")]
            for r in het:
                rn = r.get_resname().strip()
                if rn not in IGNORE_HET:
                    het_counter[rn] += 1
            if poly:
                nums = [r.id[1] for r in poly]
                rng = f"{min(nums)}-{max(nums)}"
            else:
                rng = "-"
            miss = missing.get(chain.id, 0)
            if poly or het:
                lines.append(f"| {chain.id} | {ctype} | {len(poly)} | {rng} | {miss or '-'} |")
        # ligand classification
        ligands = {k: v for k, v in het_counter.items() if k not in COMMON_IONS}
        ions = {k: v for k, v in het_counter.items() if k in COMMON_IONS}
        lines.append("")
        lines.append(f"- **Ligands/small molecules**: {ligands or 'none'}")
        lines.append(f"- **Metals/buffer**: {ions or 'none'}")
        total_missing = sum(missing.values())
        lines.append(f"- **Missing residues total**: {total_missing}"
                     + ("  WARNING loop completion needed" if total_missing > 0 else "  (none missing)"))
    OUT.write_text("\n".join(lines) + "\n")
    print(f"[OK] -> {OUT}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
