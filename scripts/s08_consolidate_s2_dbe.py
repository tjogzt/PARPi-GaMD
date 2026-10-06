#!/usr/bin/env python3
"""Consolidate the final 10-system S2 DBE + AAI table (data-layer version).

Reads (archived data layer; no hard-coded well depths)
------------------------------------------------------
- data/analysis/s2_dbe_well_depths.npy      7 exact-log systems (CV1/CV2 C1-C3)
- data/01_curated/nira_prot_dbe.csv         niraparib S2 rows + S1 C3, protonated re-run
- data/01_curated/rerun_cumulant_wells.csv  olaparib/rucaparib S2 rows, exact re-runs
- data/analysis/s1_dbe_unified_wells.csv    unified S1 C3, seven original systems
- data/01_curated/extension_s1_wells.csv    extension trio S1 C3

Writes
------
- data/analysis/s2_dbe_final.csv  (ligand, metric, C1, C2, C3, source)
  CV1/CV2 rows: well depths in the archived display forms (1-decimal from the
  npy panel; source strings passed through from the curated CSVs).
  AAI rows: S1 C3 / S2 CV2 C3 at 6 decimals.
  Source tags: exact_log_dbe | nira_prot_dbe | rerun_exact_dbe.

Notes
-----
- The niraparib rows inside rerun_cumulant_wells.csv are superseded drafts;
  the final niraparib values come from nira_prot_dbe.csv.
- Anchors asserted below match the manuscript display values
  (AAI talazoparib ~0.83, AZD5305 ~1.19; niraparib CV1 62.0; olaparib CV2 48.0).

Purpose:  Consolidate the final 10-system S2 DBE + AAI table from the archived data layer.
Created:  2026-09-17 (header standardised 2026-10-05; data-layer rewrite 2026-10-06)
Changelog: scripts/CHANGELOG.md — mandatory for any further change (this table is the manuscript's numerical anchor).
Depends:  common.paths, numpy
Run:      python3 scripts/s08_consolidate_s2_dbe.py   (from the repository root)
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir, REPO_ROOT

A = analysis_dir()
CUR = REPO_ROOT / "data" / "01_curated"


def read_rows(p):
    """Read a small curated CSV into a list of dicts (exact display strings)."""
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


wells = np.load(str(A / "s2_dbe_well_depths.npy"), allow_pickle=True).item()
nira = read_rows(CUR / "nira_prot_dbe.csv")
rerun = read_rows(CUR / "rerun_cumulant_wells.csv")
unified = read_rows(A / "s1_dbe_unified_wells.csv")
ext = read_rows(CUR / "extension_s1_wells.csv")


def d1(v):
    """1-decimal display form used by the exact-log panel."""
    return f"{float(v):.1f}"


# --- CV rows ------------------------------------------------------------------
EXACT = ["APO", "talazoparib", "AZD5305", "veliparib",
         "fluzoparib", "pamiparib", "senaparib"]
rows_cv = {}
for lig in EXACT:
    rows_cv[lig] = {cvn: [d1(x) for x in wells[f"{lig}_{cvn}"]] for cvn in ("CV1", "CV2")}
for lig, key, src in [("niraparib", "sys2_niraparib_prot", nira),
                      ("olaparib", "sys2_olaparib", rerun),
                      ("rucaparib", "sys2_rucaparib", rerun)]:
    # column-name case differs between the two curated sources (c1/c2/c3 vs C1/C2/C3)
    rows_cv[lig] = {r["cv"]: [r.get("C1", r.get("c1")), r.get("C2", r.get("c2")),
                              r.get("C3", r.get("c3"))]
                    for r in src if r["system"] == key}

# --- S1 C3 numerators (display forms preserved from the sources) --------------
s1 = {r["ligand"].lower(): r["C3"] for r in unified}
s1["niraparib"] = next(r["c3"] for r in nira if r["system"] == "sys1_niraparib_prot")
for r in ext:
    s1[r["ligand"].lower()] = d1(r["C3"])

# --- assemble -----------------------------------------------------------------
ORDER = ["APO", "talazoparib", "AZD5305", "veliparib", "fluzoparib",
         "pamiparib", "senaparib", "niraparib", "olaparib", "rucaparib"]
TAG = {l: "exact_log_dbe" for l in EXACT}
TAG.update({"niraparib": "nira_prot_dbe", "olaparib": "rerun_exact_dbe",
            "rucaparib": "rerun_exact_dbe"})
rows = []
for lig in ORDER:
    for cvn in ("CV1", "CV2"):
        rows.append([lig, cvn, *rows_cv[lig][cvn], TAG[lig]])
    s1v = float(s1[lig.lower()])
    s2v = float(rows_cv[lig]["CV2"][2])
    rows.append([lig, "AAI", s1[lig.lower()], rows_cv[lig]["CV2"][2],
                 f"{s1v / s2v:.6f}", TAG[lig]])

# --- manuscript anchors -------------------------------------------------------
anch = {(r[0], r[1]): r for r in rows}
assert len(rows) == 30
assert anch[("niraparib", "CV1")][4] == "62.0"
assert anch[("niraparib", "CV2")][4] == "50.5"
assert anch[("olaparib", "CV2")][4] == "48.0"
assert anch[("rucaparib", "CV1")][4] == "59.3"
assert anch[("talazoparib", "AAI")][4].startswith("0.830")
assert anch[("AZD5305", "AAI")][4].startswith("1.186")

with open(A / "s2_dbe_final.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["ligand", "metric", "C1", "C2", "C3", "source"])
    w.writerows(rows)
print(f"wrote {A/'s2_dbe_final.csv'} ({len(rows)} rows)")
