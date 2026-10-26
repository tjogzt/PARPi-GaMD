#!/usr/bin/env python3
"""s73_tala26_structure_qc.py — structural QC for the independent 26-ns
talazoparib seed-2026 trajectory.

Production frames only (step >= 6.0e6; frames 240-759 of the 760-frame DCD,
520 frames at 50 ps). Computes:
  - CA backbone RMSD after Kabsch superposition to the first production frame;
  - Zn-SG minimum distances and coordination numbers for Zn1-Zn3;
  - HD/ART domain self-aligned per-residue RMSF (same protocol as
    code/36_recompute_s2_profiles.py);
  - HD x ART residue-level DCCM block.

Outputs:
  data/analysis/s73_tala26_qc_frames.csv
  data/analysis/s73_tala26_rmsf_residue.csv
  data/analysis/s73_tala26_dccm_hd_art.csv
  data/analysis/s73_tala26_structure_qc.csv
"""
import csv
import sys
import warnings
from pathlib import Path

import MDAnalysis as mda
import numpy as np

warnings.filterwarnings("ignore")

ROOT = Path("/Users/taozhu/clacky_workspace/PARPi_design/PARPi-GaMD")
sys.path.insert(0, str(ROOT))
from common.kabsch import kabsch_align  # noqa: E402

PRMTOP = Path("/Volumes/tjogzt4T/PARPi_data/sys2_talazoparib/sys2_talazoparib.prmtop")
DCD = Path("/Volumes/tjogzt4T/PARPi_data/extension/tala_s2_26ns/output.dcd")
OUTDIR = ROOT / "data/analysis"
HD = (662, 787)
ART = (788, 1014)
PROD_STEP = 6_000_000
N_PROD = 520


def kabsch_rmsd(P, Q):
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _s, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((Pc @ R.T - Qc) ** 2, axis=1))))


def main():
    u = mda.Universe(str(PRMTOP), str(DCD))
    ca = u.select_atoms("protein and name CA")
    zn = u.select_atoms("name ZN")
    sg = u.select_atoms("name SG")
    hd_sel = u.select_atoms(f"resid {HD[0]}:{HD[1]} and name CA")
    art_sel = u.select_atoms(f"resid {ART[0]}:{ART[1]} and name CA")
    resids = ca.residues.resids
    print(f"frames={len(u.trajectory)} CA={len(ca)} Zn={len(zn)} SG={len(sg)} "
          f"HD={len(hd_sel)} ART={len(art_sel)}", flush=True)
    assert len(zn) == 3 and len(sg) > 0

    # DCD frames k=1..760; production starts at frame 240 (step 6.0e6).
    prod = list(range(240, 240 + N_PROD))
    u.trajectory[prod[0]]
    ref_ca = ca.positions.copy()
    ref_all = ref_ca - ref_ca.mean(axis=0)
    ref_hd = hd_sel.positions.copy()
    ref_hd = ref_hd - ref_hd.mean(axis=0)
    ref_art = art_sel.positions.copy()
    ref_art = ref_art - ref_art.mean(axis=0)

    sp_hd = np.zeros((len(hd_sel), 3))
    sp2_hd = np.zeros_like(sp_hd)
    sp_art = np.zeros((len(art_sel), 3))
    sp2_art = np.zeros_like(sp_art)
    xyz = np.zeros((N_PROD, len(ca), 3))
    rows = []
    for k, fi in enumerate(prod):
        u.trajectory[fi]
        rmsd = kabsch_rmsd(ca.positions, ref_ca)
        zpos = zn.positions.copy()
        spos = sg.positions.copy()
        d = np.linalg.norm(zpos[:, None, :] - spos[None, :, :], axis=2)
        mins = d.min(axis=1)
        ncoord = (d <= 2.8).sum(axis=1)
        rows.append((fi, (fi - 239) * 0.05, rmsd, *mins, *ncoord))
        xyz[k] = kabsch_align(ca.positions.copy(), ref_all)
        ah = kabsch_align(hd_sel.positions.copy(), ref_hd)
        sp_hd += ah
        sp2_hd += ah ** 2
        aa = kabsch_align(art_sel.positions.copy(), ref_art)
        sp_art += aa
        sp2_art += aa ** 2
        if (k + 1) % 50 == 0:
            print(f"  {k+1}/{N_PROD} frames", flush=True)

    mp_hd = sp_hd / N_PROD
    rmsf_hd = np.sqrt((sp2_hd / N_PROD - mp_hd ** 2).sum(axis=1))
    mp_art = sp_art / N_PROD
    rmsf_art = np.sqrt((sp2_art / N_PROD - mp_art ** 2).sum(axis=1))

    hd_m = (resids >= HD[0]) & (resids <= HD[1])
    art_m = (resids >= ART[0]) & (resids <= ART[1])

    # Whole-protein-aligned DCCM, residue-level displacement Pearson.
    X = xyz  # frames x CA x 3
    Xm = X - X.mean(axis=0, keepdims=True)
    cov = np.einsum("tia,tja->ij", Xm, Xm)
    sd = np.sqrt(np.diag(cov))
    R = cov / np.outer(sd, sd)
    block = R[np.ix_(hd_m, art_m)]

    OUTDIR.mkdir(parents=True, exist_ok=True)
    with open(OUTDIR / "s73_tala26_qc_frames.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frame_index", "time_ns", "rmsd_ca", "zn1_min", "zn2_min",
                    "zn3_min", "zn1_ncoord", "zn2_ncoord", "zn3_ncoord"])
        for r in rows:
            w.writerow([r[0], f"{r[1]:.1f}", f"{r[2]:.3f}",
                        *[f"{x:.3f}" for x in r[3:6]], *r[6:9]])

    with open(OUTDIR / "s73_tala26_rmsf_residue.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["resid", "rmsf_hd", "rmsf_art"])
        hd_ids = resids[hd_m]
        art_ids = resids[art_m]
        for rid, v in zip(hd_ids, rmsf_hd):
            w.writerow([int(rid), f"{v:.4f}", ""])
        for rid, v in zip(art_ids, rmsf_art):
            w.writerow([int(rid), "", f"{v:.4f}"])

    np.savetxt(OUTDIR / "s73_tala26_dccm_hd_art.csv", block, delimiter=",",
               fmt="%.5f", header="HD rows x ART columns; resid order "
               f"HD {HD[0]}-{HD[1]}, ART {ART[0]}-{ART[1]}")

    summary = {
        "system": "talazoparib_seed2026_26ns",
        "n_frames": N_PROD,
        "rmsd_first": rows[0][2],
        "rmsd_last": rows[-1][2],
        "rmsd_max": max(r[2] for r in rows),
        "zn_min": min(min(r[3:6]) for r in rows),
        "zn_max": max(max(r[3:6]) for r in rows),
        "zn1_ncoord_min": min(r[6] for r in rows),
        "zn2_ncoord_min": min(r[7] for r in rows),
        "zn3_ncoord_min": min(r[8] for r in rows),
        "HD_rmsf": float(rmsf_hd.mean()),
        "HD_sd": float(rmsf_hd.std()),
        "ART_rmsf": float(rmsf_art.mean()),
        "ART_sd": float(rmsf_art.std()),
        "dccm_abs": float(np.abs(block).mean()),
        "dccm_mean": float(block.mean()),
        "dccm_pos_pct": float((block > 0.3).mean() * 100),
        "dccm_neg_pct": float((block < -0.3).mean() * 100),
    }
    with open(OUTDIR / "s73_tala26_structure_qc.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary))
        w.writeheader()
        w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v)
                    for k, v in summary.items()})
    print("\nsummary")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"\nwrote {OUTDIR}/s73_tala26_*.csv")


if __name__ == "__main__":
    main()
