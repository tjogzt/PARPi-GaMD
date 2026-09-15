"""Ligand panel metadata (single source: common/ligands.csv).

Columns: ligand, label, class, color, lty, shape, trap_potency
(trap_potency = trapping relative to olaparib; None = not determined).
"""
import csv
from pathlib import Path

CSV = Path(__file__).with_name("ligands.csv")


def load():
    with CSV.open(newline="") as fh:
        rows = []
        for r in csv.DictReader(fh):
            r["trap_potency"] = float(r["trap_potency"]) if r["trap_potency"] else None
            r["shape"] = int(r["shape"])
            rows.append(r)
        return rows


LIGANDS = {r["ligand"]: r for r in load()}


def colors():
    return {lig: r["color"] for lig, r in LIGANDS.items()}
