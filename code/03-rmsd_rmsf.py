#!/usr/bin/env python3
"""
03-rmsd_rmsf.py — RMSD and RMSF analysis (in-memory, no temp files).
Stride=100. Domain-level decomposition: HD, ART, DNA.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def log(msg: str) -> None:
    print(f"[rmsd] {msg}", flush=True)

def align_and_get_rmsd(mobile, ref_coords, mobile_coords):
    """Kabsch superposition + RMSD between mobile_coords and ref_coords.
    Both N×3 arrays. Returns (rmsd, rotated_mobile)."""
    from scipy.spatial.transform import Rotation as R
    # Center
    ref_c = ref_coords - ref_coords.mean(axis=0)
    mob_c = mobile_coords - mobile_coords.mean(axis=0)
    # Covariance matrix
    C = mob_c.T @ ref_c
    V, S, Wt = np.linalg.svd(C)
    # Rotation matrix
    rot = Wt.T @ V.T
    # Ensure proper rotation (det=1)
    if np.linalg.det(rot) < 0:
        S[-1] *= -1
        rot = Wt.T @ np.diag(S) @ V.T
    # Apply rotation
    aligned = mob_c @ rot
    rmsd = np.sqrt(np.mean((aligned - ref_c) ** 2))
    # Return aligned coordinates in original frame's center + ref's center
    return rmsd, aligned + ref_coords.mean(axis=0)

def main():
    import MDAnalysis as mda

    DATA = Path("/Volumes/tjogzt4T/PARPi_data")
    OUT = Path("results/analysis")
    OUT.mkdir(parents=True, exist_ok=True)

    systems = ["sys2_APO", "sys2_talazoparib"]

    for name in systems:
        top = DATA / name / f"{name}.prmtop"
        traj = DATA / name / "output.dcd"
        if not top.exists() or not traj.exists():
            log(f"SKIP {name}: missing files")
            continue

        log(f"Loading {name}...")
        u = mda.Universe(str(top), str(traj))
        n_frames = len(u.trajectory)
        log(f"  frames={n_frames} atoms={len(u.atoms)}")

        # Selection groups
        ca_all = u.select_atoms("protein and name CA")
        hd_ca = u.select_atoms("resid 1:375 and name CA")   # HD domain
        art_ca = u.select_atoms("resid 390:525 and name CA") # ART domain 
        dna_p = u.select_atoms("nucleic and name P")

        n_ca = {"all": len(ca_all), "HD": len(hd_ca),
                "ART": len(art_ca), "DNA": len(dna_p)}
        for k, v in n_ca.items():
            log(f"  {k} atoms: {v}")

        stride = 100
        n_analyze = n_frames // stride

        # Pre-allocate
        rmsd = {k: np.zeros(n_analyze) for k in n_ca}
        aligned_sum = {k: np.zeros((n_ca[k], 3)) for k in n_ca}
        aligned_sum2 = {k: np.zeros((n_ca[k], 3)) for k in n_ca}

        # Reference: first frame coordinates
        ref = {
            "all": ca_all.positions.copy(),
            "HD": hd_ca.positions.copy(),
            "ART": art_ca.positions.copy(),
            "DNA": dna_p.positions.copy() if n_ca["DNA"] > 0 else None,
        }

        log(f"  Processing {n_analyze} frames (stride={stride})...")
        for frame_idx, ts in enumerate(u.trajectory[::stride]):
            # All-protein-CA alignment (overall superposition)
            rmsd_all, aligned_all = align_and_get_rmsd(
                ca_all, ref["all"], ca_all.positions.copy())
            rmsd["all"][frame_idx] = rmsd_all
            aligned_sum["all"] += aligned_all
            aligned_sum2["all"] += aligned_all ** 2

            # HD domain RMSD (using same overall rotation applied to CA, then subset)
            # For domain RMSD: align on domain CA, compute within-domain RMSD
            rmsd_hd, aligned_hd = align_and_get_rmsd(
                hd_ca, ref["HD"], hd_ca.positions.copy())
            rmsd["HD"][frame_idx] = rmsd_hd
            aligned_sum["HD"] += aligned_hd
            aligned_sum2["HD"] += aligned_hd ** 2

            rmsd_art, aligned_art = align_and_get_rmsd(
                art_ca, ref["ART"], art_ca.positions.copy())
            rmsd["ART"][frame_idx] = rmsd_art
            aligned_sum["ART"] += aligned_art
            aligned_sum2["ART"] += aligned_art ** 2

            if n_ca["DNA"] > 0:
                rmsd_dna, aligned_dna = align_and_get_rmsd(
                    dna_p, ref["DNA"], dna_p.positions.copy())
                rmsd["DNA"][frame_idx] = rmsd_dna
                aligned_sum["DNA"] += aligned_dna
                aligned_sum2["DNA"] += aligned_dna ** 2

            if frame_idx % 50 == 0:
                parts = [f"All={rmsd_all:.2f}", f"HD={rmsd_hd:.2f}",
                         f"ART={rmsd_art:.2f}"]
                if n_ca["DNA"] > 0:
                    parts.append(f"DNA={rmsd_dna:.2f}")
                log(f"  frame {frame_idx*stride}: " + " ".join(parts))

        # RMSF from aligned coordinates
        log("  Computing RMSF...")
        rmsf = {}
        for k in n_ca:
            mean_pos = aligned_sum[k] / n_analyze
            rmsf_raw = np.sqrt(np.sum(aligned_sum2[k] / n_analyze - mean_pos**2, axis=1))
            rmsf[k] = rmsf_raw

        # Save numerical results
        for domain in n_ca:
            np.savetxt(
                str(OUT / f"{name}_rmsd_{domain}.dat"),
                np.column_stack([np.arange(n_analyze) * stride, rmsd[domain]]),
                fmt="%d %.4f", header="frame rmsd(A)")

            # For RMSF, map to residue IDs
            if domain == "all":
                res_ids = ca_all.residues.resids
            elif domain == "HD":
                res_ids = hd_ca.residues.resids
            elif domain == "ART":
                res_ids = art_ca.residues.resids
            else:  # DNA
                res_ids = dna_p.residues.resids if n_ca["DNA"] > 0 else np.arange(len(rmsf[domain]))

            np.savetxt(
                str(OUT / f"{name}_rmsf_{domain}.dat"),
                np.column_stack([res_ids, rmsf[domain]]),
                fmt="%d %.4f", header="resid rmsf(A)")

        # Summary
        log(f"\n--- {name} Summary ---")
        for domain in n_ca:
            r = rmsd[domain]
            f = rmsf[domain]
            log(f"  {domain} RMSD: mean={r.mean():.2f} ± {r.std():.2f} Å, "
                f"max={r.max():.2f}")
            log(f"  {domain} RMSF: mean={f.mean():.2f} ± {f.std():.2f} Å, "
                f"max={f.max():.2f}")

    log("\n[DONE]")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
