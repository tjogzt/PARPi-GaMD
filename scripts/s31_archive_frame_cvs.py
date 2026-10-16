#!/usr/bin/env python3
#
# Purpose:  export per-frame CV + boost energies for every
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s31_archive_frame_cvs.py   (from the repository root)

"""archive_frame_cvs.py — export per-frame CV + boost energies for every
production trajectory of the current pipeline, as the machine-readable input
of the manuscript's descriptors (one authoritative source per trajectory).

For each trajectory the script writes:
  results/archive/per_frame/{tag}_frames.csv.gz
    columns: frame, step, dV_D, dV_T, cv (S1), or cv1,cv2 (S2)
  results/archive/manifest.csv
    one row per trajectory: tag, n_prod_frames, cv_mean, cv_min, cv_max, hash

Conventions (identical to s18_rerun_exact_dbe.py / analyze_{nira,ruca}_prot.py):
  - production mask: step >= 6.0e6 (production-start-step.txt)
  - dV_D = gamd.log column 8 (index 7); dV_T = column 7 (index 6) where present
  - S1 CV = |COM(HD 1:119 CA) - COM(ART 120:345 CA)|
  - S2 CV1 = |COM(protein CA) - COM(DNA heavy)|; CV2 = |COM(HD 662:787) - COM(ART 788:1011)|

Run: python scripts/archive_frame_cvs.py [--all | --rerun-only]
"""
import argparse
import gzip
import hashlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "archive" / "per_frame"
OUT.mkdir(parents=True, exist_ok=True)
PROD_STEP = 6.0e6
DATA = data_root()
RERUN = DATA / "rerun_202609"


def md5_file(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()[:12]


def load_log(log_path):
    log = np.loadtxt(log_path, comments="#")
    prod = log[:, 1] >= PROD_STEP
    dv_d = log[:, 7]
    dv_t = log[:, 6] if log.shape[1] > 7 else np.zeros(len(log))
    return dv_d, dv_t, prod


def compute_cv(prmtop, dcd, which):
    import MDAnalysis as mda
    u = mda.Universe(str(prmtop), str(dcd))
    if which == "S1":
        hd = u.select_atoms("resid 1:119 and name CA")
        art = u.select_atoms("resid 120:345 and name CA")
        assert len(hd) > 0 and len(art) > 0
        cvs = np.array([np.linalg.norm(hd.center_of_mass() - art.center_of_mass())
                        for ts in u.trajectory])
        return {"cv": cvs}
    else:
        hd = u.select_atoms("resid 662:787 and name CA")
        art = u.select_atoms("resid 788:1011 and name CA")
        prot = u.select_atoms("protein and name CA")
        dna = u.select_atoms("(resname DA DT DG DC or resname DA5 DT5 DG5 DC5 or "
                             "resname DA3 DT3 DG3 DC3) and not element H")
        assert len(hd) > 0 and len(art) > 0 and len(dna) > 0
        cv1 = np.array([np.linalg.norm(prot.center_of_mass() - dna.center_of_mass())
                        for ts in u.trajectory])
        cv2 = np.array([np.linalg.norm(hd.center_of_mass() - art.center_of_mass())
                        for ts in u.trajectory])
        return {"cv1": cv1, "cv2": cv2}


def archive(tag, prmtop, dcd, log, which, manifest_rows):
    print(f"[{tag}] log={log.name} dcd={dcd.name}", flush=True)
    dv_d_all, dv_t_all, prod_mask = load_log(log)
    cvs = compute_cv(prmtop, dcd, which)
    # align DCD frames to log rows (both span the full run), then take production
    n_full = min(len(dv_d_all), len(next(iter(cvs.values()))))
    for k in cvs:
        cvs[k] = cvs[k][:n_full][prod_mask[:n_full]]
    dv_d = dv_d_all[:n_full][prod_mask[:n_full]]
    dv_t = dv_t_all[:n_full][prod_mask[:n_full]]
    n = len(dv_d)
    steps = np.arange(n) + 1
    cols = [np.arange(n) + 1, steps, dv_d]
    header = "frame,step_ps,dV_D"
    if which == "S1" and np.abs(dv_t).max() > 0:
        cols.append(dv_t)
        header += ",dV_T"
    if which == "S1":
        cols.append(cvs["cv"])
        header += ",cv_hd_art_com"
    else:
        cols += [cvs["cv1"], cvs["cv2"]]
        header += ",cv1_prot_dna_com,cv2_hd_art_com"
    arr = np.column_stack(cols)
    out_path = OUT / f"{tag}_frames.csv.gz"
    with gzip.open(out_path, "wt") as f:
        f.write(header + "\n")
        np.savetxt(f, arr, fmt="%.6f")
    h = md5_file(out_path)
    cv_col = arr[:, -1] if which == "S1" else arr[:, 3]
    manifest_rows.append({
        "tag": tag, "construct": which, "n_prod_frames": int(n),
        "cv_mean": round(float(cv_col.mean()), 4),
        "cv_min": round(float(cv_col.min()), 4),
        "cv_max": round(float(cv_col.max()), 4),
        "archive_hash": h,
    })
    print(f"    -> {out_path.name}: {n} frames, cv {cv_col.min():.1f}..{cv_col.max():.1f}",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rerun-only", action="store_true")
    args = ap.parse_args()
    manifest = []

    # 1) six rerun trajectories (nira/ola/ruca x S1/S2)
    for lig in ["niraparib", "olaparib", "rucaparib"]:
        for which, tag in [("S1", f"sys1_{lig}"), ("S2", f"sys2_{lig}")]:
            if which == "S1":
                prmtop = RERUN / "sys_files" / f"{tag}.prmtop"
                dcd = RERUN / "runs" / tag / "gamd_out" / "output.dcd"
                log = RERUN / "rerun_logs" / f"{tag}_gamd.log"
            else:
                d = RERUN / "runs" / tag
                prmtop = next(d.glob("*.prmtop"), None)
                dcd = d / "gamd_out" / "output.dcd"
                log = RERUN / "rerun_logs" / f"{tag}_gamd.log"
            if prmtop and dcd.exists() and log.exists():
                archive(tag, prmtop, dcd, log, which, manifest)

    # 2) protonated re-simulations
    for run_dir, tag_base in [(RERUN / "runs_nira_prot", "nira_prot"),
                              (RERUN / "runs_ruca_prot", "ruca_prot")]:
        for which, sysn in [("S1", "sys1"), ("S2", "sys2")]:
            d = run_dir / sysn
            log = d / "gamd.log"
            dcd = d / "output.dcd"
            prmtop = next(d.glob("*.prmtop"), None)
            if prmtop and dcd.exists() and log.exists():
                archive(f"{sysn}_{tag_base}", prmtop, dcd, log, which, manifest)

    # 3) original-panel systems (skip with --rerun-only)
    if not args.rerun_only:
        orig = {
            ("S1", "APO"): (DATA / "sys1_APO", None),
            ("S2", "APO"): (DATA / "sys2_APO", None),
            ("S1", "AZD5305"): (DATA / "sys1_AZD5305", None),
            ("S2", "AZD5305"): (DATA / "sys2_AZD5305", None),
            ("S2", "talazoparib"): (DATA / "sys2_talazoparib", None),
            ("S2", "veliparib"): (DATA / "sys2_veliparib", None),
            ("S2", "niraparib"): (DATA / "sys2_niraparib", None),
            ("S2", "olaparib"): (DATA / "sys2_olaparib", None),
            ("S2", "rucaparib"): (DATA / "sys2_rucaparib", None),
        }
        for (which, lig), (d, _) in orig.items():
            if not d.exists():
                continue
            prmtop = next(d.glob("*.prmtop"), None)
            dcd = next(d.glob("**/output.dcd"), None)
            log = next(d.glob("**/gamd.log"), None)
            if prmtop and dcd and log:
                try:
                    archive(f"orig_{which.lower()}_{lig}", prmtop, dcd, log,
                            which, manifest)
                except Exception as e:
                    print(f"[orig_{which}_{lig}] SKIP: {e}", flush=True)

    import csv
    mpath = ROOT / "results" / "archive" / "manifest.csv"
    with open(mpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["tag", "construct", "n_prod_frames",
                                          "cv_mean", "cv_min", "cv_max",
                                          "archive_hash"])
        w.writeheader()
        for row in manifest:
            w.writerow(row)
    print(f"WROTE {mpath} ({len(manifest)} trajectories)")


if __name__ == "__main__":
    main()
