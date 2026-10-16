#!/usr/bin/env python3
#
# Purpose:  system-1 pocket docking + redocking validation
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   Usage: python run_docking.py <receptor.pdbqt> <box_file> <ref_sdf> [exhaustiveness]
# Depends:  Python 3 standard library
# Run:      python3 scripts/s12_run_docking.py

"""
run_docking.py -- system-1 pocket docking + redocking validation

Docks every ligand in docking/ligands/*.pdbqt into the specified receptor pocket with AutoDock Vina,
aggregates the top binding affinity; for rucaparib also runs a redocking RMSD check
(top pose vs the crystallographic RPB, obrms symmetry-aware) as a docking-protocol quality gate.

Usage: python run_docking.py <receptor.pdbqt> <box_file> <ref_sdf> [exhaustiveness]
"""
from __future__ import annotations
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIG_DIR = ROOT / "docking" / "ligands"
RES_DIR = ROOT / "docking" / "results"
SEED = 42


def parse_box(box: Path) -> dict[str, float]:
    d = {}
    for line in box.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=")
            d[k.strip()] = float(v)
    return d


def vina_dock(rec: Path, lig: Path, box: dict, out: Path, exh: int) -> float | None:
    cmd = ["vina", "--receptor", str(rec), "--ligand", str(lig),
           "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
           "--center_z", str(box["center_z"]), "--size_x", str(box["size_x"]),
           "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
           "--exhaustiveness", str(exh), "--num_modes", "9", "--seed", str(SEED),
           "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [vina FAIL] {lig.name}: {r.stderr.strip()[:200]}", file=sys.stderr)
        return None
    # top affinity from the first REMARK VINA RESULT in the output pdbqt
    for line in out.read_text().splitlines():
        m = re.search(r"VINA RESULT:\s*(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    return None


def top_pose_sdf(pdbqt: Path) -> Path:
    sdf = pdbqt.with_suffix(".top.sdf")
    subprocess.run(["obabel", str(pdbqt), "-O", str(sdf), "-f", "1", "-l", "1"],
                   capture_output=True, text=True)
    return sdf


def redock_rmsd(ref_sdf: Path, pose_sdf: Path) -> float | None:
    r = subprocess.run(["obrms", str(ref_sdf), str(pose_sdf)],
                       capture_output=True, text=True)
    m = re.search(r"(-?\d+\.\d+)", r.stdout)
    return float(m.group(1)) if m else None


def main() -> int:
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    rec = Path(sys.argv[1])
    box = parse_box(Path(sys.argv[2]))
    ref_sdf = Path(sys.argv[3])
    exh = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    RES_DIR.mkdir(parents=True, exist_ok=True)

    ligs = sorted(LIG_DIR.glob("*.pdbqt"))
    print(f"docking {len(ligs)} ligands -> {rec.name}  (exhaustiveness={exh})\n")
    results = []
    for lig in ligs:
        name = lig.stem
        out = RES_DIR / f"{name}_docked.pdbqt"
        score = vina_dock(rec, lig, box, out, exh)
        rmsd = ""
        if name == "rucaparib" and score is not None:
            pose = top_pose_sdf(out)
            val = redock_rmsd(ref_sdf, pose)
            rmsd = f"  redock RMSD={val:.2f} A" if val is not None else "  RMSD=?"
        print(f"  {name:12s} affinity = {score if score is not None else '?'} kcal/mol{rmsd}")
        results.append((name, score, rmsd.strip()))

    # aggregate csv
    out_csv = RES_DIR / "docking_summary.csv"
    with out_csv.open("w") as fh:
        fh.write("name,top_affinity_kcal_mol,note\n")
        for name, score, note in sorted(results, key=lambda x: (x[1] is None, x[1] or 0)):
            fh.write(f"{name},{score if score is not None else ''},{note}\n")
    print(f"\n[DONE] -> {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
