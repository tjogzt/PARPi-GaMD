#!/usr/bin/env python3
#
# Purpose:  Robust well-depth estimators on the C3 PMF curves (corrected definitions).
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s43_robust_span_estimators.py   (from the repository root)

"""Robust well-depth estimators on the C3 PMF curves (corrected definitions).

Manuscript convention: well depth = pmf.max() - pmf.min() over the 2-col xvg
(RC, PMF kcal/mol) — the y-range of the reweighted C3 curve.
Robust alternatives on the same curves:
  span_mm      = pmf.max() - pmf.min()                    (current convention)
  span_q90     = 5%-95% interquantile range of PMF values
  span_trunc95 = 2.5%-97.5% interquantile range of PMF values
  span_eff     = RC width (Angstrom) of the region with PMF <= min + 3.0
Then re-run the n=5 trapping association (Spearman, exact permutation p,
leave-one-out, tie sensitivity) for each estimator on S2 CV1.
Outputs results/analysis/robust_spans.csv and robust_association.csv."""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import analysis_dir

A = Path(analysis_dir())
N5 = ["talazoparib", "niraparib", "olaparib", "rucaparib", "veliparib"]
TRAP = {"talazoparib": 5.0, "niraparib": 4.0, "olaparib": 2.5,
        "rucaparib": 2.5, "veliparib": 1.0}


def load_pmf(path):
    lines = [ln for ln in open(path) if ln.strip() and not ln.startswith(("#", "@"))]
    arr = np.array([[float(x) for x in ln.split()] for ln in lines])
    return arr[:, 0], arr[:, 1]


def spans(x, f):
    return dict(
        span_mm=float(f.max() - f.min()),
        span_q90=float(np.percentile(f, 95) - np.percentile(f, 5)),
        span_trunc95=float(np.percentile(f, 97.5) - np.percentile(f, 2.5)),
        span_eff=float(x[f <= f.min() + 3.0].max() - x[f <= f.min() + 3.0].min()),
    )


def s2_path(lig, cvn):
    return A / f"pmf-c3-sys2_{lig}_{cvn}_cv.dat.xvg"


def s1_path(lig):
    return A / f"sys1_{lig}_pmf_c3.xvg"


def exact_p(x, y):
    """tie-aware exact permutation p (matches 21_spearman_correlation.R)"""
    from itertools import permutations
    xr = rankdata(x)
    yr = rankdata(y)
    obs = np.corrcoef(xr, yr)[0, 1]
    cnt = tot = 0
    for perm in set(permutations(range(len(xr)))):
        r = np.corrcoef(rankdata(xr[np.array(perm)]), yr)[0, 1]
        tot += 1
        if abs(r) >= abs(obs) - 1e-12:
            cnt += 1
    return obs, cnt / tot


def main():
    import csv
    rows = []
    for lig in ["APO", "AZD5305", "talazoparib", "veliparib",
                "niraparib", "olaparib", "rucaparib"]:
        for cvn in ["CV1", "CV2"]:
            x, f = load_pmf(s2_path(lig, cvn))
            rows.append(dict(panel="S2", lig=lig, cv=cvn, **spans(x, f)))
        x, f = load_pmf(s1_path(lig))
        rows.append(dict(panel="S1", lig=lig, cv="HD-ART", **spans(x, f)))
    with open(A / "robust_spans.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("S2 CV1 table:")
    for r in rows:
        if r["panel"] == "S2" and r["cv"] == "CV1":
            print(f"  {r['lig']:12s} mm={r['span_mm']:6.1f} q90={r['span_q90']:6.1f} "
                  f"t95={r['span_trunc95']:6.1f} eff={r['span_eff']:6.1f} A")

    d = {r["lig"]: r for r in rows if r["panel"] == "S2" and r["cv"] == "CV1"}
    assoc = []
    for est in ["span_mm", "span_q90", "span_trunc95", "span_eff"]:
        x = np.array([d[l][est] for l in N5])
        y = np.array([TRAP[l] for l in N5])
        rho, p = exact_p(x, y)
        loos = [np.corrcoef(rankdata(np.delete(x, k)),
                            rankdata(np.delete(y, k)))[0, 1] for k in range(5)]
        y2 = np.array([5.0, 4.0, 2.0, 1.0, 0.5])  # rucaparib below olaparib (manuscript convention)
        rho2, p2 = exact_p(x, y2)
        assoc.append(dict(estimator=est, rho=rho, p=p,
                          loo_min=min(loos), loo_max=max(loos),
                          tie_rho=rho2, tie_p=p2))
        print(f"\n{est}: rho={rho:+.3f} p={p:.4f} | LOO [{min(loos):+.3f}, {max(loos):+.3f}]"
              f" | tie rho={rho2:+.3f} p={p2:.4f}")
        print(f"   values: " + ", ".join(f"{l}:{d[l][est]:.1f}" for l in N5))
    with open(A / "robust_association.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(assoc[0].keys()))
        w.writeheader()
        w.writerows(assoc)
    print("\nDONE ->", A / "robust_association.csv")


if __name__ == "__main__":
    main()
