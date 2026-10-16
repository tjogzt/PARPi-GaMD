#!/usr/bin/env python3
#
# Purpose:  main-pipeline analysis of the REV-1/REV-2 new replicate trajectories (tala 3x200ns, veli 2x400ns).
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   <DATA_ROOT> new REP-1/REP-2 replicate trajectories (see docstring)
# Outputs:  Output: results/analysis/replicates/<TAG>_s<seed>/{cv.dat, weights.dat, pmf_c{1,2,3}.xvg, summary.json}
# Depends:  numpy, common, MDAnalysis
# Run:      python3 scripts/s21_new_replicate_dbe.py

"""analyze_new_replicates.py -- main-pipeline analysis of the REV-1/REV-2 new replicate trajectories (tala 3x200ns, veli 2x400ns).

Convention identical to Table S1:
  CV   = HD (internal 0-118) vs ART (internal 119-344) CA center-of-mass distance (Angstrom)
  Weight = canonical DBE: w1 = beta*dV_D (gamd.log col 8, kcal/mol); 3 columns (beta*dV, 0, dV)
  PMF  = PyReweighting-1D cumulant C1/C2/C3 -> well depth + location

Usage:
  python scripts/analyze_new_replicates.py <run_dir>...
  Each run_dir contains: TAG.prmtop, output.dcd, gamd.log
Output: results/analysis/replicates/<TAG>_s<seed>/{cv.dat, weights.dat, pmf_c{1,2,3}.xvg, summary.json}
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root  # noqa: E402

try:
    import MDAnalysis as mda
except ImportError:
    mda = None

KB = 0.0019872041  # kcal/(mol K)
TEMP = 300.0
BETA = 1.0 / (KB * TEMP)

OUT = Path(__file__).resolve().parents[1] / "results" / "analysis" / "replicates"


def extract_cv(dcd, prmtop):
    """HD 0-118 vs ART 119-344 CA COM distance (internal numbering, same convention as Table S1)."""
    u = mda.Universe(str(prmtop), str(dcd))
    hd = u.select_atoms("resid 1:119 and name CA")   # MDAnalysis 1-based
    art = u.select_atoms("resid 120:345 and name CA")
    assert len(hd) == 119 and len(art) == 226, f"unexpected CA selection: {len(hd)}/{len(art)}"
    cv = []
    for ts in u.trajectory:
        cv.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    return np.array(cv)


def parse_gamd_log(log):
    """Parse gamd.log cols 7/8, dV_T (Total-BE) / dV_D (Dih-BE) (kcal/mol).
    Dual-boost canonical weight = beta*(dV_T + dV_D)."""
    rows = []
    with open(log) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 8:
                try:
                    rows.append((float(parts[6]), float(parts[7])))
                except ValueError:
                    continue
    return np.array(rows)


def well_depth(pmf):
    """Well depth = max(PMF) - min(PMF) over the sampled RC range (same convention as S1_PMF_c3_stats.csv)."""
    if len(pmf) < 2:
        return None, None
    y = pmf[:, 1]
    return float(y.max() - y.min()), float(pmf[int(np.argmin(y)), 0])


def analyze(run_dir, prmtop=None, tag=None, seed=None):
    from common.pmf import run_pyrew

    run_dir = Path(run_dir)
    # fetch layout: review_runs/sys1_talazoparib_s49/{output.dcd,gamd.log}, prmtop in the parent directory
    name = run_dir.name
    if tag is None:
        tag = name.rsplit("_s", 1)[0] if "_s" in name else name
    if seed is None:
        seed = name.rsplit("_s", 1)[1] if "_s" in name else "na"
    if prmtop is None:
        prmtop = run_dir.parent / f"{tag}.prmtop"
    dcd = run_dir / "output.dcd"
    log = run_dir / "gamd.log"
    if not (dcd.exists() and Path(prmtop).exists() and log.exists()):
        print(f"[skip] {run_dir} missing dcd/prmtop/log")
        return
    outd = OUT / f"{tag}_s{seed}"
    outd.mkdir(parents=True, exist_ok=True)

    # cv is reusable; weights re-derived from the log each run (dual-boost canonical: beta*(dV_T+dV_D))
    if (outd / "cv.dat").exists():
        cv = np.loadtxt(outd / "cv.dat")
        print(f"{tag}_s{seed}: reuse existing cv ({len(cv)} frames)")
    else:
        cv = extract_cv(dcd, prmtop)
        np.savetxt(outd / "cv.dat", cv, fmt="%.4f")

    dv = parse_gamd_log(log)
    n = min(len(cv), len(dv))
    if len(dv) < len(cv):
        print(f"  WARNING {tag}_s{seed}: gamd.log weight gap {len(cv)-len(dv)} frames "
              f"({(len(cv)-len(dv))*50} ps) -- reconstruction from trajectory required")
    cv, dv = cv[:n], dv[:n]
    dvtot = dv[:, 0] + dv[:, 1]
    print(f"{tag}_s{seed}: frames={n}, cv={cv.min():.2f}-{cv.max():.2f} A, "
          f"dV_T={dv[:,0].mean():.2f} dV_D={dv[:,1].mean():.2f} (dual-boost total {dvtot.mean():.2f})")

    w = np.column_stack([BETA * dvtot, np.zeros(n), dvtot])
    np.savetxt(outd / "weights.dat", w, fmt="%.4f")

    summary = {"tag": tag, "seed": seed, "frames": n,
               "cv_mean": float(cv.mean()), "cv_sd": float(cv.std()),
               "dvt_mean": float(dv[:, 0].mean()), "dvd_mean": float(dv[:, 1].mean())}
    for order in ["c1", "c2", "c3"]:
        r = run_pyrew(cv, w, wfmt="%.6f", cvfmt="%.6f")
        src = r.tmpdir / f"pmf-{order}-cv.dat.xvg"
        dst = outd / f"pmf_{order}.xvg"
        dst.write_text(src.read_text())
        pmf = np.loadtxt(dst, comments=["#", "@"])
        depth, rcmin = well_depth(pmf)
        summary[f"{order}_depth"] = depth
        summary[f"{order}_rcmin"] = rcmin
    (outd / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  -> {outd}: " + " ".join(f"{k.split('_')[0].upper()}={summary.get(k)}"
          for k in ["c1_depth", "c2_depth", "c3_depth"]))
    return summary


if __name__ == "__main__":
    for rd in sys.argv[1:]:
        analyze(rd)
