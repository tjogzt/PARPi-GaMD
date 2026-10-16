#!/usr/bin/env python3
#
# Purpose:  Windowed block-bootstrap convergence diagnostics (route 1 of TECH_ROADMAP.md).
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s33_convergence_diagnostics.py   (from the repository root)

"""Windowed block-bootstrap convergence diagnostics (route 1 of TECH_ROADMAP.md).

For each key system, slice the production trajectory into 2 ns windows
(40 frames at 50 ps), reweight each window with PyReweighting C3 (DBE),
take the per-window span, then block-bootstrap the window spans to get
confidence intervals and cumulative-window convergence curves.

Systems: S2 rerun (niraparib/olaparib/rucaparib, CV1 and CV2, 13 windows)
and S1 talazoparib/veliparib 3-seed replicates (12 windows each).
Outputs results/analysis/convergence_diag.csv + *_window_spans.npz."""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir, data_root, pyrew

BETA = 1.0 / (0.001987 * 300.0)
A = Path(analysis_dir())
D = Path(data_root())
RERUN = D / 'rerun_202609'
PROD_STEP = 6.0e6
WIN = 40          # 2 ns at 50 ps/frame
NBOOT = 10000
RNG = np.random.default_rng(42)

# reuse run_and_r from collect_s1_table1 (writes weights, runs PyReweighting-1D)
import importlib.util
spec = importlib.util.spec_from_file_location("s1t1", "scripts/collect_s1_table1.py")
s1t1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s1t1)


def window_spans(cv, w, win=WIN):
    """Per-window C3 spans. w is the (beta*dV, 0, dV) weight array aligned to cv."""
    spans = []
    n = min(len(cv), len(w))
    for i in range(0, n - win + 1, win):
        d1, d2, d3, r = s1t1.run_and_r(cv[i:i + win], w[i:i + win])
        spans.append(d3)
    return np.array(spans)


def bootstrap(spans, nboot=NBOOT, rng=RNG):
    s = np.asarray(spans)
    boots = np.array([np.mean(rng.choice(s, size=len(s), replace=True))
                      for _ in range(nboot)])
    return np.percentile(boots, [2.5, 97.5])


def s2_rerun(tag):
    import MDAnalysis as mda
    lig = tag.split('_')[1]
    prmtop = D / f"sys2_{lig}" / f"sys2_{lig}.prmtop"
    u = mda.Universe(str(prmtop),
                     str(RERUN / "runs" / tag / "gamd_out" / "output.dcd"))
    hd = u.select_atoms("resid 662:787 and name CA")
    art = u.select_atoms("resid 788:1011 and name CA")
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    cv1, cv2 = [], []
    for ts in u.trajectory:
        cv1.append(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
        cv2.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    log = np.loadtxt(RERUN / "rerun_logs" / f"{tag}_gamd.log", comments="#")
    prod = log[:, 1] >= PROD_STEP
    dv = log[:, 7][prod]
    w = np.column_stack([BETA * dv, np.zeros(len(dv)), dv])
    return (np.array(cv1)[prod], np.array(cv2)[prod], w)


def s1_repl(tag, seed):
    d = D / "review_runs" / f"sys1_{tag}_{seed}"
    log = np.loadtxt(d / "gamd.log", comments="#")
    prod = log[:, 1] >= PROD_STEP
    dv = log[:, 7][prod][:480]
    cv_full = np.loadtxt(A / "replicates" / f"sys1_{tag}_{seed}" / "cv.dat")
    cv = cv_full[prod][:480]
    w = np.column_stack([BETA * dv, np.zeros(len(dv)), dv])
    return cv, w


def main():
    rows = []
    npz = {}
    # S2 rerun: CV1/CV2
    for lig in ["niraparib", "olaparib", "rucaparib"]:
        tag = f"sys2_{lig}"
        cv1, cv2, w = s2_rerun(tag)
        for cvn, cv in [("CV1", cv1), ("CV2", cv2)]:
            n = min(len(cv), len(w))
            sp = window_spans(cv[:n], w[:n])
            lo, hi = bootstrap(sp)
            npz[f"s2_{lig}_{cvn}"] = sp
            rows.append(dict(system=f"s2_{lig}", cv=cvn, nwin=len(sp),
                             span_mean=sp.mean(), span_sd=sp.std(),
                             ci_lo=lo, ci_hi=hi, ci_width=hi - lo,
                             full_span=s1t1.run_and_r(cv[:n], w[:n])[2]))
            print(f"s2_{lig} {cvn}: {len(sp)} win, "
                  f"span {sp.mean():.1f} +/- {sp.std():.1f}, "
                  f"CI [{lo:.1f}, {hi:.1f}]", flush=True)
    # S1 replicates
    for tag in ["talazoparib", "veliparib"]:
        per_seed = []
        for seed in ["s7", "s49", "s123"]:
            cv, w = s1_repl(tag, seed)
            n = min(len(cv), len(w))
            sp = window_spans(cv[:n], w[:n])
            per_seed.append(sp)
            npz[f"s1_{tag}_{seed}"] = sp
            print(f"s1_{tag}_{seed}: {len(sp)} win, span {sp.mean():.1f} +/- {sp.std():.1f}",
                  flush=True)
        allsp = np.concatenate(per_seed)
        lo, hi = bootstrap(allsp)
        seed_means = [s.mean() for s in per_seed]
        rows.append(dict(system=f"s1_{tag}", cv="HD-ART", nwin=len(allsp),
                         span_mean=allsp.mean(), span_sd=allsp.std(),
                         ci_lo=lo, ci_hi=hi, ci_width=hi - lo,
                         full_span=np.mean(seed_means)))
        print(f"s1_{tag} 3-seed: span {allsp.mean():.1f} +/- {allsp.std():.1f}, "
              f"CI [{lo:.1f}, {hi:.1f}], seed means {['%.1f' % s for s in seed_means]}",
              flush=True)
    import csv
    with open(A / "convergence_diag.csv", "w", newline="") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(rows)
    np.savez(A / "convergence_diag_window_spans.npz", **npz)
    print("DONE ->", A / "convergence_diag.csv")


if __name__ == "__main__":
    main()
