#!/usr/bin/env python3
#
# Purpose:  Three-state panel: exact permutation Spearman under table convention.
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s41_three_state_audit.py   (from the repository root)

"""Three-state panel: exact permutation Spearman under table convention.
States: (A) original neutral; (B) niraparib protonated; (C) rucaparib protonated.
All values = table C3 convention (Table S1 / Table S5).
"""
import itertools
import numpy as np

def spearman_exact(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    rx = np.argsort(np.argsort(x)) + 1
    ry = np.argsort(np.argsort(y)) + 1
    # tie-aware midranks
    def midrank(v):
        r = np.argsort(np.argsort(v)) + 1.0
        for val in np.unique(v):
            idx = np.where(v == val)[0]
            if len(idx) > 1:
                r[idx] = r[idx].mean()
        return r
    rx = midrank(x); ry = midrank(y)
    rho = np.corrcoef(rx, ry)[0, 1]
    obs = abs(rho)
    cnt = 0; tot = 0
    for perm in itertools.permutations(range(n)):
        tot += 1
        rp = np.corrcoef(rx[list(perm)], ry)[0, 1]
        if abs(rp) >= obs - 1e-12:
            cnt += 1
    return rho, cnt / tot

# trapping ranks: tala 5, nira 4, ola 2.5, ruca 2.5, veli 1
y = [5, 4, 2.5, 2.5, 1]

# Table values
s2_neu = [53.1, 56.6, 58.5, 59.3, 58.9]     # tala, nira, ola, ruca, veli
s2_nira_prot = [53.1, 62.0, 58.5, 59.3, 58.9]
s2_ruca_prot = [53.1, 62.0, 58.5, 79.7, 58.9]  # ruca 59.3 -> 79.7

s1_neu = [51.4, 49.4, 51.6, 48.7, 62.9]
s1_nira_prot = [51.4, 59.3, 51.6, 48.7, 62.9]
s1_ruca_prot = [51.4, 59.3, 51.6, 61.0, 62.9]  # ruca 48.7 -> 61.0

aai_neu = [51.4/61.9, 49.4/53.4, 51.6/48.0, 48.7/52.6, 62.9/72.9]
aai_nira_prot = [51.4/61.9, 59.3/50.5, 51.6/48.0, 48.7/52.6, 62.9/72.9]
aai_ruca_prot = [51.4/61.9, 59.3/50.5, 51.6/48.0, 61.0/75.4, 62.9/72.9]

print("metric | state           | rho    | p")
print("-" * 55)
for name, vals in [
    ("S2 CV1", s2_neu), ("S2 CV1", s2_nira_prot), ("S2 CV1", s2_ruca_prot),
    ("S1     ", s1_neu), ("S1     ", s1_nira_prot), ("S1     ", s1_ruca_prot),
    ("AAI    ", aai_neu), ("AAI    ", aai_nira_prot), ("AAI    ", aai_ruca_prot),
]:
    state = ["neutral  ", "nira prot", "ruca prot"][[0, 1, 2][[s2_neu, s2_nira_prot, s2_ruca_prot, s1_neu, s1_nira_prot, s1_ruca_prot, aai_neu, aai_nira_prot, aai_ruca_prot].index(vals) // 3]] if False else ["neutral", "nira_prot", "ruca_prot"][([s2_neu, s2_nira_prot, s2_ruca_prot, s1_neu, s1_nira_prot, s1_ruca_prot, aai_neu, aai_nira_prot, aai_ruca_prot].index(vals)) % 3]
    rho, p = spearman_exact(vals, y)
    print(f"{name} | {state:9s} | {rho:+.3f} | {p:.3f}")

# LOO for ruca-prot panel, S2 CV1
print("\nLOO (ruca-prot, S2 CV1):")
for i, lbl in enumerate(["tala", "nira", "ola", "ruca", "veli"]):
    v = [x for j, x in enumerate(s2_ruca_prot) if j != i]
    yy = [x for j, x in enumerate(y) if j != i]
    rho, p = spearman_exact(v, yy)
    print(f"  excl {lbl:5s}: rho={rho:+.3f} p={p:.3f}")

print("\nwindow noise check (S1): 61.0 vs neutral 48.7 -> +12.3 kcal shift")
print("S2 CV1 shift: 59.3 -> 79.7 (+20.4); S2 CV2 shift: 52.6 -> 75.4 (+22.8)")
