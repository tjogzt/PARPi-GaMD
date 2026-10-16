#!/usr/bin/env python3
"""Consolidate the frozen S2 DBE + AAI table from the corrected-span data layer.

Reads (archived data layer; no hard-coded well depths)
------------------------------------------------------
- data/analysis/r6/r6c_s2_spanfix.csv        CV1 C1-C3 spans, cutoff-2 convention (R6)
- data/analysis/r6/r6c_s2_cv2_spanfix.csv    CV2 C1-C3 spans (R6)
- data/analysis/r6/r6_panel_values.csv       veliparib-protonated S1 C3 (R9)
- data/analysis/r6/r6c_s1_spanfix.csv        niraparib-protonated S1 C3 (R6)
- data/analysis/s1_dbe_unified_wells.csv     unified S1 C3, seven original systems
- data/01_curated/extension_s1_wells.csv    extension trio S1 C3
- data/analysis/veli_prot/veli_prot_s2_summary.json  veliparib-protonated S2 rows (R9)

Writes
------
- data/analysis/s2_dbe_final.csv  (ligand, metric, C1, C2, C3, source)
  CV1/CV2 rows: occupied-support spans at 1-decimal display precision.
  AAI rows: S1 C3 / S2 CV2 C3 at 6 decimals.
  Source tags: exact_log_dbe | nira_prot_dbe | rerun_exact_dbe | veli_prot_dbe.

Notes
-----
- Asserts pin the frozen manuscript state; if the frozen table is present, the
  rebuilt table must match it cell-for-cell (golden-file check) or the script
  fails loudly (this table is the manuscript's numerical anchor).

Purpose:  Consolidate the frozen 11-system S2 DBE + AAI table from the corrected-span data layer.
Created:  2026-09-16 (header standardised 2026-10-05; corrected-span rewrite 2026-10-11)
Changelog: scripts/CHANGELOG.md — mandatory for any further change (this table is the manuscript's numerical anchor).
Depends:  csv, json
Run:      python3 scripts/s08_consolidate_s2_dbe.py   (from the repository root)
"""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(".")
A = ROOT / "data/analysis"
C = ROOT / "data/01_curated"


def load_rows(p, key):
    with open(p, newline="") as fh:
        return {r[key]: r for r in csv.DictReader(fh)}


def d1(x):
    """1-decimal display convention used by the frozen table."""
    return f"{float(x):.1f}"


def pick(d, *keys):
    for k in keys:
        if k in d:
            return d[k]
    raise KeyError(f"none of {keys} in {sorted(d)}")


s2cv1 = load_rows(A / "r6/r6c_s2_spanfix.csv", "system")
s2cv2 = load_rows(A / "r6/r6c_s2_cv2_spanfix.csv", "system")
s1sp = load_rows(A / "r6/r6c_s1_spanfix.csv", "system")
panel = load_rows(A / "r6/r6_panel_values.csv", "ligand")
uni = {r["ligand"].lower(): r for r in csv.DictReader(open(A / "s1_dbe_unified_wells.csv"))}
ext = load_rows(C / "extension_s1_wells.csv", "ligand")
v2 = json.loads((A / "veli_prot/veli_prot_s2_summary.json").read_text())

# entity -> (display label, source tag) — frozen order (veliparib-prot last)
ENT = [("APO", "APO", "exact_log_dbe"),
       ("talazoparib", "talazoparib", "exact_log_dbe"),
       ("AZD5305", "AZD5305", "exact_log_dbe"),
       ("veliparib", "veliparib", "exact_log_dbe"),
       ("fluzoparib", "fluzoparib", "exact_log_dbe"),
       ("pamiparib", "pamiparib", "exact_log_dbe"),
       ("senaparib", "senaparib", "exact_log_dbe"),
       ("niraparib_prot", "niraparib", "nira_prot_dbe"),
       ("olaparib", "olaparib", "rerun_exact_dbe"),
       ("rucaparib", "rucaparib", "rerun_exact_dbe")]

rows = []
for ent, lab, tag in ENT:
    cv1 = [lab, "CV1", d1(s2cv1[ent]["c1_valid"]), d1(s2cv1[ent]["c2_valid"]),
           d1(s2cv1[ent]["c3_valid"]), tag]
    cv2 = [lab, "CV2", d1(s2cv2[ent]["c1"]), d1(s2cv2[ent]["c2"]),
           d1(s2cv2[ent]["c3"]), tag]
    s1 = d1(uni[lab.lower()]["C3"]) if lab.lower() in ("apo", "olaparib", "rucaparib") \
        else {"talazoparib": d1(uni["talazoparib"]["C3"]),
              "azd5305": d1(uni["azd5305"]["C3"]),
              "veliparib": d1(uni["veliparib"]["C3"]),
              "fluzoparib": d1(ext["fluzoparib"]["C3"]),
              "pamiparib": d1(ext["pamiparib"]["C3"]),
              "senaparib": d1(ext["senaparib"]["C3"]),
              "niraparib": d1(s1sp["niraparib_prot"]["pmf_c3_v2"])}[lab.lower()]
    rows += [cv1, cv2,
             [lab, "AAI", s1, cv2[4], f"{float(s1) / float(cv2[4]):.6f}", tag]]

# veliparib-protonated rows from the archived R9 summary (display label: veliparib)
vp_s1 = panel["veliparib_prot"]["S1_C3"]
for cvn in ("CV1", "CV2"):
    blk = v2[cvn]
    rows.append(["veliparib", cvn, d1(pick(blk, "C1", "c1")), d1(pick(blk, "C2", "c2")),
                 d1(pick(blk, "C3", "c3")), "veli_prot_dbe"])
cv2c3 = rows[-1][4]
rows.append(["veliparib", "AAI", vp_s1, cv2c3, f"{float(vp_s1) / float(cv2c3):.6f}", "veli_prot_dbe"])

# ---- frozen-state assertions ----
idx = {(r[0], r[1]): r for r in rows if r[5] == "exact_log_dbe" or r[5] == "nira_prot_dbe"
       or r[5] == "rerun_exact_dbe"}
assert len(rows) == 33, len(rows)
assert idx[("APO", "CV1")][2:5] == ["6.4", "11.3", "19.5"]
assert idx[("talazoparib", "CV1")][2:5] == ["4.7", "12.2", "22.9"]
assert idx[("AZD5305", "CV1")][4] == "49.5"
assert idx[("niraparib", "CV1")][2:5] == ["3.9", "12.1", "34.1"]
assert idx[("olaparib", "CV2")][4] == "16.6"
assert idx[("rucaparib", "CV1")][4] == "28.7"
assert idx[("talazoparib", "AAI")][4].startswith("1.1625")
assert idx[("AZD5305", "AAI")][4].startswith("1.8423")
assert idx[("APO", "AAI")][4].startswith("4.2888")
assert rows[-1][5] == "veli_prot_dbe" and rows[-1][4].startswith("1.410")
assert rows[-2][4] == "16.8" and rows[-1][:5] == ["veliparib", "AAI", "23.7", "16.8", "1.410714"]

# ---- golden-file check / write ----
dst = A / "s2_dbe_final.csv"
if dst.exists():
    with open(dst, newline="") as fh:
        cur = [r for r in csv.reader(fh)][1:]
    if len(cur) != len(rows) or any(a != b for a, b in zip(cur, rows)):
        for a, b in zip(cur, rows):
            if a != b:
                print("MISMATCH:", a, "vs", b, file=sys.stderr)
        sys.exit("golden-file check failed")
    print(f"golden check passed ({len(rows)} rows identical to the frozen table)")
else:
    with open(dst, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand", "metric", "C1", "C2", "C3", "source"])
        w.writerows(rows)
    print(f"wrote {dst} ({len(rows)} rows)")
