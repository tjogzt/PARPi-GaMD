#!/usr/bin/env python3
"""s55_s10_archive_frame_cvs.py — per-frame CV + boost archive for the S10 cross-seed runs.

Purpose:  Builds the per-frame analysis records for the two additional seeds
          (7 and 123) and two chemical states of the niraparib S2 system
          (26-ns production each; identical protocol), used for the cross-seed
          stability analysis of the primary S2 CV1 span (SI Table S17).
          Pairing is step-based, not row-based: DCD frame k <-> simulation step
          25,000*k, with gamd.log rows looked up by their step column. This is
          required for the seed-123 neutral run, whose gamd.log lost ~51
          buffered rows (steps 16.35-17.60 M, 2.55 ns) when the pre-migration
          process was stopped; its DCD is complete (760 frames = 38 ns).
Conventions (identical to s31_archive_frame_cvs.py):
          CV1 = |COM(protein CA) - COM(DNA heavy)|;
          CV2 = |COM(HD 662:787 CA) - COM(ART 788:1011 CA)|;
          dV_D = gamd.log column 8 (index 7); dV_T = index 6 (reported);
          production mask: step >= 6.0e6 (production-start-step.txt);
          frame values are written at 50-ps cadence.
Inputs:   runs/sys2_niraparib/{gamd_out_s7_r14, gamd_out_s123_r14}/
          (output.dcd, gamd.log) and runs/sys2_niraparib/sys2_niraparib.prmtop;
          runs_nira_prot/sys2_niraparib_prot/{gamd_out_s7_r14,
          gamd_out_s123_r14}/ and build/niraparib_best_prot.prmtop
Outputs:  results/analysis/s10/{tag}_frames.csv.gz (production frames only);
          curated copies: data/per_frame/{tag}_frames.csv.gz
Depends:  numpy, MDAnalysis
Run:      python3 scripts/s55_s10_archive_frame_cvs.py   (from the repository root)
"""
import gzip
import hashlib
from pathlib import Path

import numpy as np

try:
    import MDAnalysis as mda
except ImportError:  # pragma: no cover - environment check
    raise SystemExit("MDAnalysis is required (conda/pip); see environment_gamd_min.yml")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "analysis" / "s10"
OUT.mkdir(parents=True, exist_ok=True)
PROD_STEP = 6.0e6
STEP_PER_FRAME = 25_000  # 50 ps at 2-fs integration

RUNS = {
    "s10_neu_s7": dict(
        prmtop=ROOT / "runs/sys2_niraparib/sys2_niraparib.prmtop",
        rundir=ROOT / "runs/sys2_niraparib/gamd_out_s7_r14"),
    "s10_prot_s7": dict(
        prmtop=ROOT / "runs_nira_prot/sys2_niraparib_prot/build/niraparib_best_prot.prmtop",
        rundir=ROOT / "runs_nira_prot/sys2_niraparib_prot/gamd_out_s7_r14"),
    "s10_neu_s123": dict(
        prmtop=ROOT / "runs/sys2_niraparib/sys2_niraparib.prmtop",
        rundir=ROOT / "runs/sys2_niraparib/gamd_out_s123_r14"),
    "s10_prot_s123": dict(
        prmtop=ROOT / "runs_nira_prot/sys2_niraparib_prot/build/niraparib_best_prot.prmtop",
        rundir=ROOT / "runs_nira_prot/sys2_niraparib_prot/gamd_out_s123_r14"),
}


def compute_cv(prmtop, dcd):
    """CV1/CV2 COMD series (identical selections to s31_archive_frame_cvs.py)."""
    u = mda.Universe(str(prmtop), str(dcd))
    hd = u.select_atoms("resid 662:787 and name CA")
    art = u.select_atoms("resid 788:1011 and name CA")
    prot = u.select_atoms("protein and name CA")
    dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                         "resname DA3 DT3 DG3 DC3) and not element H")
    assert len(hd) > 0 and len(art) > 0 and len(prot) > 0 and len(dna) > 0, \
        (len(hd), len(art), len(prot), len(dna))
    cv1, cv2 = [], []
    for _ in u.trajectory:
        cv1.append(np.linalg.norm(prot.center_of_mass() - dna.center_of_mass()))
        cv2.append(np.linalg.norm(hd.center_of_mass() - art.center_of_mass()))
    return np.asarray(cv1), np.asarray(cv2)


def main():
    manifest = []
    for tag, cfg in RUNS.items():
        dcd = cfg["rundir"] / "output.dcd"
        logp = cfg["rundir"] / "gamd.log"
        log = np.loadtxt(logp, comments="#")
        dvd_map = {}
        for s, dv, dvt in zip(log[:, 1], log[:, 7], log[:, 6]):
            dvd_map[int(s)] = (float(dv), float(dvt))
        cv1, cv2 = compute_cv(cfg["prmtop"], dcd)
        rows, missing = [], []
        for k in range(1, len(cv1) + 1):
            s = STEP_PER_FRAME * k
            if s in dvd_map:
                dv, dvt = dvd_map[s]
                rows.append((k, s * 2e-6, dv, dvt, cv1[k - 1], cv2[k - 1]))  # time in ns
            else:
                missing.append(s)
        arr = np.asarray(rows)
        arrp = arr[arr[:, 1] >= PROD_STEP * 2e-6]     # production starts at 12 ns
        outp = OUT / f"{tag}_frames.csv.gz"
        with gzip.open(outp, "wt") as f:
            f.write("frame,time_ns,dV_D,dV_T,cv1_prot_dna_com,cv2_hd_art_com\n")
            np.savetxt(f, arrp, fmt="%.6f")
        h = hashlib.md5(outp.read_bytes()).hexdigest()[:12]
        print(f"[{tag}] log rows={len(log)} missing={len(missing)} "
              f"({len(missing) * 0.05:.2f} ns) paired={len(arr)} "
              f"prod frames={len(arrp)} cv1 mean={arrp[:, 4].mean():.1f} md5={h}",
              flush=True)
        if missing:
            print(f"        missing step range: {min(missing):,}..{max(missing):,}")
        manifest.append(dict(tag=tag, n_prod=len(arrp), n_missing=len(missing),
                             t_first=round(float(arrp[0, 1]), 3),
                             t_last=round(float(arrp[-1, 1]), 3),
                             cv1_mean=round(float(arrp[:, 4].mean()), 3),
                             cv1_min=round(float(arrp[:, 4].min()), 3),
                             cv1_max=round(float(arrp[:, 4].max()), 3),
                             cv2_mean=round(float(arrp[:, 5].mean()), 3), md5=h))
    import csv as _csv
    with open(OUT / "s10_run_manifest.csv", "w", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=list(manifest[0].keys()))
        w.writeheader()
        w.writerows(manifest)
    print("DONE ->", OUT)


if __name__ == "__main__":
    main()
