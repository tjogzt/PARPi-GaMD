#!/usr/bin/env python3
"""
build_seed_dataset.py — PARP1 trapping seed dataset builder

Purpose:   Validate the curated PARP1-inhibitor seed ligands and compute RDKit
           physicochemical descriptors for the trapping ML seed dataset.
Inputs:    data/01_curated/seed_ligands.csv
           (compound identity + trapping labels; all data externalized)
Outputs:   data/01_curated/trapping_seed_dataset.csv
Depends:   rdkit

Notes:
- trapping_class (low/medium/high) is a field-consensus qualitative label
  (Murai 2012/2014, Hopkins 2015, and others; see quantitative_trapping_data.md).
  Exact quantitative values (fold-trapping / EC50) remain to be extracted from
  full texts.
- time_split: train = first disclosure <= 2021; test = 2022+ (blind test of
  next-generation selective PARP1 inhibitors).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    from rdkit import RDLogger
except ImportError:
    sys.exit("RDKit not found. `mamba install -c conda-forge rdkit`")

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "data" / "01_curated" / "seed_ligands.csv"
OUT = ROOT / "data" / "01_curated" / "trapping_seed_dataset.csv"

COLS = ["name", "synonyms", "pubchem_cid", "molecular_formula", "isomeric_smiles",
        "inchikey_pubchem", "trapping_class", "trapping_rank", "parp1_selective",
        "first_disclosure_year", "time_split", "primary_pmid",
        "zandarashvili_type", "catalytic_ic50_nM_murai", "notes",
        "rdkit_canonical_smiles", "inchikey_rdkit", "inchikey_match",
        "mw", "clogp", "hbd", "hba", "tpsa", "rotatable_bonds",
        "aromatic_rings", "stereo_centers", "n_heavy", "n_rings"]


def main() -> int:
    with IN.open(newline="") as fh:
        compounds = list(csv.DictReader(fh))

    required = ["compound_id", "isomeric_smiles", "inchikey_pubchem",
                "pubchem_cid", "trapping_class", "source_verified"]
    for c in compounds:
        for key in required:
            if not c.get(key, "").strip():
                print(f"[FAIL] {c.get('compound_id', '?')}: missing {key}", file=sys.stderr)
                return 1

    rows = []
    n_fail = 0
    for c in compounds:
        name = c["compound_id"]
        smi = c["isomeric_smiles"]
        ikey_pc = c["inchikey_pubchem"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"[FAIL] {name}: RDKit cannot parse SMILES", file=sys.stderr)
            n_fail += 1
            continue
        can = Chem.MolToSmiles(mol)
        ikey_rd = Chem.MolToInchiKey(mol)
        match = (ikey_rd == ikey_pc)
        if not match:
            print(f"[WARN] {name}: InChIKey mismatch PubChem={ikey_pc} RDKit={ikey_rd}",
                  file=sys.stderr)
        stereo = len(Chem.FindMolChiralCenters(mol, useLegacyImplementation=False,
                                                includeUnassigned=True))
        year = int(c["first_disclosure_year"])
        time_split = "test" if year >= 2022 else "train"

        def to_float(v):
            return None if not v.strip() else round(float(v), 2)

        rows.append({
            "name": name, "synonyms": c["synonyms"],
            "pubchem_cid": int(c["pubchem_cid"]),
            "molecular_formula": c["molecular_formula"],
            "isomeric_smiles": smi, "inchikey_pubchem": ikey_pc,
            "trapping_class": c["trapping_class"],
            "trapping_rank": int(c["trapping_rank"]),
            "parp1_selective": c["parp1_selective"] == "TRUE",
            "first_disclosure_year": year, "time_split": time_split,
            "primary_pmid": c["primary_pmid"],
            "zandarashvili_type": c["zandarashvili_type"],
            "catalytic_ic50_nM_murai": to_float(c["catalytic_ic50_nM_murai"]),
            "notes": c["notes"],
            "rdkit_canonical_smiles": can, "inchikey_rdkit": ikey_rd,
            "inchikey_match": match,
            "mw": round(Descriptors.MolWt(mol), 2),
            "clogp": round(Crippen.MolLogP(mol), 2),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "n_heavy": Descriptors.HeavyAtomCount(mol),
            "n_rings": rdMolDescriptors.CalcNumRings(mol),
            "stereo_centers": stereo,
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)

    print(f"[OK] wrote {len(rows)} compounds -> {OUT}")
    n_train = sum(r["time_split"] == "train" for r in rows)
    n_test = sum(r["time_split"] == "test" for r in rows)
    n_sel = sum(r["parp1_selective"] for r in rows)
    print(f"     time-split: train={n_train} test={n_test} | PARP1-selective={n_sel}")
    print(f"     InChIKey match: {sum(r['inchikey_match'] for r in rows)}/{len(rows)}")
    if n_fail:
        print(f"[WARN] {n_fail} compound(s) failed to parse", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
