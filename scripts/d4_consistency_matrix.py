#!/usr/bin/env python3
"""D4-2: estimator x metric consistency matrix (rho grid) + heatmap data.

Panel: niraparib-corrected (S2 curves are the protonated-niraparib xvg files;
S1 niraparib uses the rebuilt protonated curve, see S1_TAGS). Tab rows are
the corrected-panel table-convention values from
scripts/audit_r1_table_convention.py (table_convention_rows.csv).

Rows: metric (S2 CV1 C3, S2 CV2 C3, S1 C3, AAI)
Cols: estimator (min-max, 5-95% quantile, 95% truncated, effective support)
Cells: tie-aware exact-permutation rho (p) of the n=5 span values vs trapping rank.

AAI cells use matched estimators in numerator (S1) and denominator (S2 CV2).

Purpose:  Estimator x metric consistency matrix (rho grid) for the corrected panel; feeds the consistency-matrix figure (code/29).
Created:  packaged 2026-10-06
Depends:  numpy, scipy
Run:      python3 scripts/d4_consistency_matrix.py   (from the repository root)
"""
import numpy as np
from pathlib import Path
from itertools import permutations
from scipy.stats import rankdata

ks = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
trap_rank = np.array([5.0, 4.0, 2.5, 2.5, 1.0])
A = Path("results/analysis")
PERMS = np.array([p for p in permutations(range(5))])

S2_TAGS = {"talazoparib": "sys2_talazoparib", "niraparib": "sys2_niraparib",
           "olaparib": "sys2_olaparib", "rucaparib": "sys2_rucaparib",
           "veliparib": "sys2_veliparib"}
S1_TAGS = {"talazoparib": "sys1_talazoparib_pmf_c3", "niraparib": "sys1_niraparib_pmf_c3_prot",
           "olaparib": "sys1_olaparib_pmf_c3", "rucaparib": "sys1_rucaparib_pmf_c3",
           "veliparib": "sys1_veliparib_pmf_c3"}
# NOTE (panel audit 2026-09-29): the S2 xvg files on disk already carry the
# niraparib-protonated curves (gen_nira_prot_xvgs.py output), so S2 rows are
# the corrected panel; the S1 niraparib curve is the rebuilt protonated one
# (scripts/rebuild_s1_nira_prot.py, archived per-frame data), so all rows
# now render the CORRECTED panel to match the main-text diagnostic (2).


def read_pmf(path):
    lines = [ln for ln in open(path) if ln.strip() and not ln.startswith(('#', '@'))]
    arr = np.array([[float(x) for x in ln.split()] for ln in lines])
    return arr[:, 0], arr[:, 1]  # RC coordinate, PMF (kcal/mol)


def estimators(rc, pmf):
    mm = pmf.max() - pmf.min()
    q90 = np.percentile(pmf, 95) - np.percentile(pmf, 5)
    t95 = np.percentile(pmf, 97.5) - np.percentile(pmf, 2.5)
    sel = pmf <= pmf.min() + 3.0
    eff = rc[sel].max() - rc[sel].min()  # RC width of the PMF<=min+3 region
    return {"mm": mm, "q90": q90, "trunc95": t95, "eff": eff}


def exact_p(values):
    v = np.array(values)
    rank_v = rankdata(v, method='average')
    rank_y = rankdata(trap_rank, method='average')
    rho_obs = np.corrcoef(rank_v, rank_y)[0, 1]
    rho_perm = np.array([np.corrcoef(rankdata(v[p], method='average'), rank_y)[0, 1]
                         for p in PERMS])
    p = (np.abs(rho_perm) >= abs(rho_obs) - 1e-12).mean()
    return rho_obs, p


# collect per-system PMFs
s2cv1, s2cv2, s1 = {}, {}, {}
for lig in ks:
    s2cv1[lig] = read_pmf(A / f"pmf-c3-{S2_TAGS[lig]}_CV1_cv.dat.xvg")
    s2cv2[lig] = read_pmf(A / f"pmf-c3-{S2_TAGS[lig]}_CV2_cv.dat.xvg")
    s1[lig] = read_pmf(A / f"{S1_TAGS[lig]}.xvg")

metrics = {"S2_CV1_C3": s2cv1, "S2_CV2_C3": s2cv2, "S1_C3": s1}
est_names = ["mm", "q90", "trunc95", "eff"]

rows = []
print(f"{'metric':12s}" + "".join(f"{e:>14s}" for e in est_names))
for mname, pmfs in metrics.items():
    line = f"{mname:12s}"
    row = {"metric": mname}
    for e in est_names:
        spans = [estimators(*pmfs[lig])[e] for lig in ks]
        rho, p = exact_p(spans)
        row[e + "_rho"] = round(rho, 3)
        row[e + "_p"] = round(p, 4)
        line += f"{rho:>9.3f}/{p:<4.2f}"
    rows.append(row)
    print(line)

# AAI: matched estimator S1/S2CV2
line = f"{'AAI':12s}"
row = {"metric": "AAI"}
for e in est_names:
    aai_vals = [estimators(*s1[lig])[e] / estimators(*s2cv2[lig])[e] for lig in ks]
    rho, p = exact_p(aai_vals)
    row[e + "_rho"] = round(rho, 3)
    row[e + "_p"] = round(p, 4)
    line += f"{rho:>9.3f}/{p:<4.2f}"
rows.append(row)
print(line)

# Table-convention fixed rows: written by scripts/audit_r1_table_convention.py
# (manuscript Table 1 / Table S1-S5 C3 values, tie-aware exact permutation).
# Read from CSV rather than embedded, so the matrix renders both the
# mean-curve convention (computed above) and the table convention side by side.
import csv as _csv
TAB_ROWS = []
with open(A / "table_convention_rows.csv") as f:
    for rec in _csv.DictReader(f):
        if rec["panel"] != "corrected":
            continue  # figure renders the corrected panel (main-text diagnostic (2))
        vals = {}
        for e in est_names:
            if e == "mm":
                vals[f"{e}_rho"] = float(rec[f"mm_rho"])
                vals[f"{e}_p"] = float(rec[f"mm_p"])
            else:
                vals[f"{e}_rho"] = float("nan")
                vals[f"{e}_p"] = float("nan")
        row = {"metric": rec["metric"], **vals}
        TAB_ROWS.append(row)
        print(f"{rec['metric']:12s} table-convention row (corrected)")
for row in TAB_ROWS:
    rows.append(row)

with open(A / "d4_consistency_matrix.csv", "w") as f:
    f.write("metric," + ",".join(f"{e}_rho,{e}_p" for e in est_names) + "\n")
    for row in rows:
        f.write(row["metric"] + "," + ",".join(
            f"{row[e+'_rho']},{row[e+'_p']}" for e in est_names) + "\n")
print("\nWROTE results/analysis/d4_consistency_matrix.csv")
