#!/usr/bin/env python3
"""
19-rmsf_all_systems.py — Per-residue RMSF for all 7 S2 systems.
Correct PARP1 domain boundaries (UniProt P09874):
  ZnF1: 1-96, ZnF2: 97-214, ZnF3: 215-383
  BRCT: 384-517, WGR: 518-661
  HD: 662-787, ART: 788-1014
Aligns on protein CA, computes per-residue RMSF per domain.
Uses APO prmtop + stripped DCDs for all systems.
"""
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import warnings
import sys
warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.kabsch import kabsch_rmsd
from common.paths import data_root

DATA = data_root()
OUT = Path("results/analysis")
OUT.mkdir(parents=True, exist_ok=True)
APO_TOP = DATA / "sys2_APO" / "sys2_APO.prmtop"

# Domain definitions (residue IDs in S2 prmtop, 1-indexed)
DOMAINS = {
    "HD":    (662, 787),
    "ART":   (788, 1014),
    "ZnF":   (1, 383),    # ZnF1+ZnF2+ZnF3
    "HDART": (662, 1014), # HD+ART combined (catalytic region)
}

# Systems: (name, dcd_path, stride_for_output)
# APO and talazoparib have their own prmtop + full DCDs
# Others use APO prmtop + stripped DCDs
systems = {}
for name in ["sys2_APO", "sys2_talazoparib"]:
    top = DATA / name / f"{name}.prmtop"
    dcd = DATA / name / "output.dcd"
    if top.exists() and dcd.exists():
        systems[name] = (str(top), str(dcd))

for name in ["sys2_AZD5305", "sys2_veliparib", "sys2_niraparib",
             "sys2_olaparib", "sys2_rucaparib"]:
    dcd = DATA / name / "output_stripped.dcd"
    if dcd.exists():
        systems[name] = (str(APO_TOP), str(dcd))

print(f"Systems to process: {len(systems)}")
for name in systems:
    print(f"  {name}")

# Align function: common/kabsch.kabsch_rmsd (imported above)

# Process each system
for sys_name, (top_path, dcd_path) in systems.items():
    out_prefix = OUT / sys_name
    rmsf_file = OUT / f"{sys_name}_rmsf_residue.dat"

    if rmsf_file.exists():
        print(f"[SKIP] {sys_name}: RMSF data exists")
        continue

    print(f"\n{'='*50}")
    print(f"Processing {sys_name}...")

    u = mda.Universe(top_path, dcd_path)
    n_frames = len(u.trajectory)
    print(f"  Frames: {n_frames}")

    # Selection: all protein CA
    ca_all = u.select_atoms("protein and name CA")
    n_ca = len(ca_all)
    ca_indices = ca_all.indices  # 0-based atom indices
    ca_resids = ca_all.residues.resids  # residue IDs

    print(f"  CA atoms: {n_ca}, residues: {len(ca_resids)}")

    # Reference: first frame
    u.trajectory[0]
    ref_coords = ca_all.positions.copy()

    # Accumulators
    n_used = 0
    rmsd_all = np.zeros(n_frames)
    sum_pos = np.zeros((n_ca, 3))
    sum_pos2 = np.zeros((n_ca, 3))

    print(f"  Computing RMSD + accumulating for RMSF...")
    for i, ts in enumerate(u.trajectory):
        rmsd_val, aligned = kabsch_rmsd(ca_all.positions.copy(), ref_coords)
        rmsd_all[i] = rmsd_val
        sum_pos += aligned
        sum_pos2 += aligned ** 2
        n_used += 1

        if (i + 1) % 1000 == 0:
            print(f"    frame {i+1}/{n_frames}, RMSD={rmsd_val:.2f} Å")

    # Compute RMSF
    mean_pos = sum_pos / n_used
    rmsf = np.sqrt(sum_pos2 / n_used - mean_pos ** 2)
    rmsf_per_res = np.sqrt(np.sum(rmsf, axis=1))  # per-residue: sqrt of xyz variance

    # Save per-residue RMSF
    with open(rmsf_file, 'w') as f:
        f.write("# resid rmsf(A)\n")
        for resid, val in zip(ca_resids, rmsf_per_res):
            f.write(f"{resid} {val:.4f}\n")
    print(f"  Saved: {rmsf_file}")

    # Domain-level summaries
    for dom_name, (start, end) in DOMAINS.items():
        mask = (ca_resids >= start) & (ca_resids <= end)
        dom_rmsf = rmsf_per_res[mask]
        dom_resids = ca_resids[mask]
        dom_mean = np.mean(dom_rmsf)
        dom_std = np.std(dom_rmsf)

        # Save domain-specific file
        dom_file = OUT / f"{sys_name}_rmsf_{dom_name}.dat"
        with open(dom_file, 'w') as f:
            f.write(f"# {dom_name} domain (resid {start}-{end}), mean={dom_mean:.3f}±{dom_std:.3f}\n")
            f.write("# resid rmsf(A)\n")
            for resid, val in zip(dom_resids, dom_rmsf):
                f.write(f"{resid} {val:.4f}\n")

        print(f"  {dom_name}: mean RMSF={dom_mean:.3f}±{dom_std:.3f} Å")

    # Overall summary
    print(f"  Total: mean RMSF={np.mean(rmsf_per_res):.3f}±{np.std(rmsf_per_res):.3f} Å")
    print(f"  Mean RMSD: {np.mean(rmsd_all):.2f}±{np.std(rmsd_all):.2f} Å")

print("\n" + "="*50)
print("All RMSF computations complete!")
