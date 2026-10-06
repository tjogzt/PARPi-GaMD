#!/usr/bin/env python3
"""R1 audit: recompute all associations from TABLE C3 values (single source).
Table S1 C3 (S1 span), Table S5 C3 (S2 CV1 / CV2). Neutral + protonated niraparib.
Tie-aware midranks, exact permutation p over 5! = 120.

Purpose:  R1-audit table-convention recomputation of the associations from table C3 values (co-produces the *_tab rows of data/analysis/d4_consistency_matrix.csv).
Created:  packaged 2026-10-06
Depends:  numpy
Run:      python3 scripts/audit_r1_table_convention.py   (from the repository root)
"""
import itertools
import numpy as np

def midrank(v):
    order = np.argsort(v, kind="stable")
    r = np.empty(len(v))
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        r[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return r

def spearman_exact(xs, ys):
    rx = midrank(np.array(xs, float)); ry = midrank(np.array(ys, float))
    rho = float(np.corrcoef(rx, ry)[0, 1]) if rx.std() and ry.std() else 0.0
    cnt = tot = 0
    for perm in itertools.permutations(xs):
        tot += 1
        rp = float(np.corrcoef(midrank(np.array(perm, float)), ry)[0, 1])
        if abs(rp) >= abs(rho) - 1e-12:
            cnt += 1
    return rho, cnt / tot

ks = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
trap = {"talazoparib": 5, "niraparib": 4, "olaparib": 2.5, "rucaparib": 2.5, "veliparib": 1}
y = [trap[k] for k in ks]

# Table S5 C3: S2 CV1 (neutral nira 56.6; protonated 62.0), S2 CV2 (neutral 53.4; prot 50.5)
s2cv1_neu = {"talazoparib": 53.1, "niraparib": 56.6, "olaparib": 58.5, "rucaparib": 59.3, "veliparib": 58.9}
s2cv1_pro = dict(s2cv1_neu, niraparib=62.0)
s2cv2_neu = {"talazoparib": 61.9, "niraparib": 53.4, "olaparib": 48.0, "rucaparib": 52.6, "veliparib": 72.9}
s2cv2_pro = dict(s2cv2_neu, niraparib=50.5)
# Table S1 C3: S1 span
s1_neu = {"talazoparib": 51.4, "niraparib": 49.4, "olaparib": 51.6, "rucaparib": 48.7, "veliparib": 62.9}
s1_pro = dict(s1_neu, niraparib=59.3)

def aai(s1, s2):
    return {k: s1[k] / s2[k] for k in ks}

print(f"{'metric':<18}{'neutral':>24}{'protonated':>24}")
for name, neu, pro in [
    ("S2 CV1 span", [s2cv1_neu[k] for k in ks], [s2cv1_pro[k] for k in ks]),
    ("S1 span",      [s1_neu[k] for k in ks],    [s1_pro[k] for k in ks]),
    ("AAI = S1C3/S2CV2C3", [aai(s1_neu, s2cv2_neu)[k] for k in ks],
                          [aai(s1_pro, s2cv2_pro)[k] for k in ks]),
]:
    rn, pn = spearman_exact(neu, y)
    rp, pp = spearman_exact(pro, y)
    print(f"{name:<18}{rn:+.3f}/p={pn:.3f}{rp:+.3f}/p={pp:.3f}")
    if name == "S1 span":
        s1_tab = (rn, pn)
        s1_tab_pro = (rp, pp)
    if name == "AAI = S1C3/S2CV2C3":
        aai_tab = (rn, pn)
        aai_tab_pro = (rp, pp)

print()
print("AAI values (table convention):")
print("  neutral:   ", {k: round(aai(s1_neu, s2cv2_neu)[k], 3) for k in ks})
print("  protonated:", {k: round(aai(s1_pro, s2cv2_pro)[k], 3) for k in ks})
print()
print("S1 protonated spans:", {k: s1_pro[k] for k in ks})

# Write the table-convention rows consumed by scripts/d4_consistency_matrix.py.
# Both panels: neutral (original parameterization) and corrected (niraparib
# protonated) mm cells; the consistency matrix's *_tab rows use the corrected
# panel to match the main-text diagnostic (2) description.
with open("results/analysis/table_convention_rows.csv", "w") as f:
    f.write("panel,metric,mm_rho,mm_p\n")
    f.write(f"neutral,S1_C3_tab,{s1_tab[0]:.3f},{s1_tab[1]:.4f}\n")
    f.write(f"neutral,AAI_tab,{aai_tab[0]:.3f},{aai_tab[1]:.4f}\n")
    f.write(f"corrected,S1_C3_tab,{s1_tab_pro[0]:.3f},{s1_tab_pro[1]:.4f}\n")
    f.write(f"corrected,AAI_tab,{aai_tab_pro[0]:.3f},{aai_tab_pro[1]:.4f}\n")
print("WROTE results/analysis/table_convention_rows.csv")
