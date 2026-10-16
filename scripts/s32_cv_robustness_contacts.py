#!/usr/bin/env python3
#
# Purpose:  CV-robustness supplement: HD-ART CA contact number and min-distance as
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s32_cv_robustness_contacts.py   (from the repository root)

"""CV-robustness supplement: HD-ART CA contact number and min-distance as
alternative (unreweighted, geometric) collective variables.

Scope: seven original-panel S2 CA trajectories (md_analysis/sys2_*_ca.dcd,
protein CA only) + the six rerun and four protonated trajectories (full DCD).
For each system the script computes two per-frame CVs:
  cv_contact = number of HD(662-787)-ART(788-1011) CA pairs within 8 A
  cv_mindist = minimum HD-ART CA distance
and reports the production-segment span (max-min) and 5-95% quantile range,
plus the n=5 rank correlation against the trapping rank order for the
neutral panel (tala/veli/nira/ola/ruca) — a descriptive CV-choice check that
does NOT enter the manuscript's table convention (which is C3-cumulant).
Output: results/analysis/cv_robustness_contacts.csv
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root

DATA = data_root()
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis"
CUTOFF = 8.0  # Angstrom, CA-CA

# HD/ART residue indices (1-based) for the S2 numbering
HD = slice(661, 787)   # 662-787
ART = slice(787, 1011)  # 788-1011


def calc(u, hd_idx, art_idx, prod_frac=0.5):
    """Contact number and min distance per frame over the last prod_frac of frames."""
    n = len(u.trajectory)
    start = int(n * (1 - prod_frac))
    contacts, mindists = [], []
    for i, ts in enumerate(u.trajectory):
        if i < start:
            continue
        x = u.atoms.positions
        d = np.linalg.norm(x[hd_idx][:, None, :] - x[art_idx][None, :, :], axis=2)
        contacts.append((d < CUTOFF).sum())
        mindists.append(d.min())
    return np.array(contacts), np.array(mindists)


def spearman_exact(x, y):
    from itertools import permutations
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    def midrank(v):
        r = np.argsort(np.argsort(v)) + 1.0
        for val in np.unique(v):
            idx = np.where(v == val)[0]
            if len(idx) > 1:
                r[idx] = r[idx].mean()
        return r
    rx, ry = midrank(x), midrank(y)
    rho = np.corrcoef(rx, ry)[0, 1]
    cnt = tot = 0
    for perm in permutations(range(n)):
        tot += 1
        if abs(np.corrcoef(rx[list(perm)], ry)[0, 1]) >= abs(rho) - 1e-12:
            cnt += 1
    return rho, cnt / tot


def main():
    import MDAnalysis as mda

    trap_rank = {"talazoparib": 5.0, "niraparib": 4.0, "olaparib": 2.5,
                 "rucaparib": 2.5, "veliparib": 1.0}
    rows = []
    print(f"{'system':22s} {'n_frames':>8s} {'cont_span':>9s} {'cont_q90':>8s} "
          f"{'mind_span':>9s} {'mind_q90':>8s}", flush=True)

    # 1) original-panel CA trajectories (protein CA only; residue order from prmtop)
    for lig in ["APO", "AZD5305", "talazoparib", "veliparib", "niraparib",
                "olaparib", "rucaparib"]:
        dcd = DATA / "md_analysis" / f"sys2_{lig}_ca.dcd"
        prm = next((DATA / f"sys2_{lig}").glob("*.prmtop"), None)
        if not dcd.exists() or prm is None:
            print(f"{lig:22s} CA dcd/prmtop missing", flush=True)
            continue
        u = mda.Universe(str(dcd))
        ref = mda.Universe(str(prm))
        resid = ref.select_atoms("protein and name CA").resids
        assert len(resid) == u.atoms.n_atoms, \
            f"{lig}: CA count mismatch {len(resid)} vs {u.atoms.n_atoms}"
        hd_idx = np.where((resid >= 662) & (resid <= 787))[0]
        art_idx = np.where((resid >= 788) & (resid <= 1011))[0]
        assert len(hd_idx) > 0 and len(art_idx) > 0, f"{lig}: empty domain slices"
        c, m = calc(u, hd_idx, art_idx)
        row = [f"orig_{lig}", len(c), c.max() - c.min(),
               np.percentile(c, 95) - np.percentile(c, 5),
               m.max() - m.min(), np.percentile(m, 95) - np.percentile(m, 5)]
        rows.append(row)
        print(f"{row[0]:22s} {row[1]:8d} {row[2]:9.1f} {row[3]:8.1f} "
              f"{row[4]:9.2f} {row[5]:8.2f}", flush=True)

    # 2) rerun + protonated trajectories (full DCD; HD/ART via resname numbering)
    import gzip
    def calc_full(prmtop, dcd, prod_step_frac=0.5):
        u = mda.Universe(str(prmtop), str(dcd))
        sel = u.select_atoms("protein and name CA")
        resid = sel.resids
        hd_idx = np.where((resid >= 662) & (resid <= 787))[0]
        art_idx = np.where((resid >= 788) & (resid <= 1011))[0]
        n = len(u.trajectory)
        start = int(n * (1 - prod_step_frac))
        contacts, mindists = [], []
        for i, ts in enumerate(u.trajectory):
            if i < start:
                continue
            x = sel.positions
            d = np.linalg.norm(x[hd_idx][:, None, :] - x[art_idx][None, :, :], axis=2)
            contacts.append((d < CUTOFF).sum())
            mindists.append(d.min())
        return np.array(contacts), np.array(mindists)

    rerun = data_root() / "rerun_202609"
    for lig in ["niraparib", "olaparib", "rucaparib"]:
        d = rerun / "runs" / f"sys2_{lig}"
        prmtop = next(d.glob("*.prmtop"), None)
        dcd = d / "gamd_out" / "output.dcd"
        if not (prmtop and dcd.exists()):
            continue
        c, m = calc_full(prmtop, dcd)
        row = [f"rerun_{lig}", len(c), c.max() - c.min(),
               np.percentile(c, 95) - np.percentile(c, 5),
               m.max() - m.min(), np.percentile(m, 95) - np.percentile(m, 5)]
        rows.append(row)
        print(f"{row[0]:22s} {row[1]:8d} {row[2]:9.1f} {row[3]:8.1f} "
              f"{row[4]:9.2f} {row[5]:8.2f}", flush=True)
    for tag, sub in [("nira_prot", "runs_nira_prot"), ("ruca_prot", "runs_ruca_prot")]:
        d = rerun / sub / "sys2"
        prmtop = next(d.glob("*.prmtop"), None)
        dcd = d / "output.dcd"
        if prmtop and dcd.exists():
            c, m = calc_full(prmtop, dcd)
            row = [f"{tag}", len(c), c.max() - c.min(),
                   np.percentile(c, 95) - np.percentile(c, 5),
                   m.max() - m.min(), np.percentile(m, 95) - np.percentile(m, 5)]
            rows.append(row)
            print(f"{row[0]:22s} {row[1]:8d} {row[2]:9.1f} {row[3]:8.1f} "
                  f"{row[4]:9.2f} {row[5]:8.2f}", flush=True)

    # 3) rank correlation (neutral panel, contact-number span)
    names = [r[0] for r in rows]
    def get(name, col):
        return rows[names.index(name)][col]
    neutral = {"talazoparib": "orig_talazoparib", "veliparib": "orig_veliparib",
               "niraparib": "rerun_niraparib", "olaparib": "rerun_olaparib",
               "rucaparib": "rerun_rucaparib"}
    if all(v in names for v in neutral.values()):
        x = [get(neutral[k], 2) for k in ["talazoparib", "niraparib", "olaparib",
                                          "rucaparib", "veliparib"]]
        y = [trap_rank[k] for k in ["talazoparib", "niraparib", "olaparib",
                                    "rucaparib", "veliparib"]]
        rho, p = spearman_exact(x, y)
        print(f"\nNeutral panel, contact-number span vs trapping rank: "
              f"rho={rho:+.3f}, exact p={p:.3f}", flush=True)
        x2 = [get(neutral[k], 4) for k in ["talazoparib", "niraparib", "olaparib",
                                           "rucaparib", "veliparib"]]
        rho2, p2 = spearman_exact(x2, y)
        print(f"Neutral panel, min-distance span vs trapping rank:    "
              f"rho={rho2:+.3f}, exact p={p2:.3f}", flush=True)

    out = OUT / "cv_robustness_contacts.csv"
    with open(out, "w") as f:
        f.write("system,n_frames,contact_span,contact_q90,min_dist_span,min_dist_q90\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    print(f"WROTE {out}", flush=True)


if __name__ == "__main__":
    main()
