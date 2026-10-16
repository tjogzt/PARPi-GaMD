#!/usr/bin/env python3
#
# Purpose:  Recompute RMSF (HD/ART) + DCCM (HDxART) for the three re-simulated S2 systems (niraparib, olaparib, rucaparib) from the 2026-09 rerun trajectories (26 ns production, 50 ps frame rate, seed 42) under the identical prot...
# Created:  2026-09 (header standardised 2026-10-11)
# Inputs:   (see docstring)
# Depends:  numpy, MDAnalysis, common
# Run:      python3 scripts/s25_rerun3_rmsf_dccm.py

"""Recompute RMSF (HD/ART) + DCCM (HDxART) for the three re-simulated S2 systems
(niraparib, olaparib, rucaparib) from the 2026-09 rerun trajectories (26 ns
production, 50 ps frame rate, seed 42) under the identical protocol used by
recompute_all_s2.py for the other four systems."""
import warnings
warnings.filterwarnings('ignore')
import sys
from pathlib import Path
import numpy as np
import MDAnalysis as mda

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.kabsch import kabsch_align
from common.paths import data_root

DATA = data_root()
RERUN = DATA / 'rerun_202609'
HD = (662, 787); ART = (788, 1014)

# Re-simulated systems: (name, prmtop, dcd, production frames)
SYS = [
    ('niraparib', DATA/'sys2_niraparib'/'sys2_niraparib.prmtop',
     RERUN/'runs'/'sys2_niraparib'/'gamd_out'/'output.dcd'),
    ('olaparib', DATA/'sys2_olaparib'/'sys2_olaparib.prmtop',
     RERUN/'runs'/'sys2_olaparib'/'gamd_out'/'output.dcd'),
    ('rucaparib', DATA/'sys2_rucaparib'/'sys2_rucaparib.prmtop',
     RERUN/'runs'/'sys2_rucaparib'/'gamd_out'/'output.dcd'),
]
N_PROD = 520  # production segment: last 520 frames = 26 ns at 50 ps/frame


def analyze(name, top, dcd):
    u = mda.Universe(str(top), str(dcd))
    ca = u.select_atoms('protein and name CA')
    hd_sel = u.select_atoms(f'resid {HD[0]}:{HD[1]} and name CA')
    art_sel = u.select_atoms(f'resid {ART[0]}:{ART[1]} and name CA')
    resids = ca.residues.resids
    ntot = len(u.trajectory)
    if N_PROD is None:
        frames = list(range(ntot))
    else:
        frames = list(range(ntot - N_PROD, ntot))
    n = len(frames)
    u.trajectory[frames[0]]
    ref_all = ca.positions.copy(); rc_all = ref_all - ref_all.mean(axis=0)
    ref_hd = hd_sel.positions.copy(); rc_hd = ref_hd - ref_hd.mean(axis=0)
    ref_art = art_sel.positions.copy(); rc_art = ref_art - ref_art.mean(axis=0)
    sp_hd = np.zeros((len(hd_sel), 3)); sp2_hd = np.zeros_like(sp_hd)
    sp_art = np.zeros((len(art_sel), 3)); sp2_art = np.zeros_like(sp_art)
    xyz = np.zeros((n, len(ca), 3))
    for k, fi in enumerate(frames):
        u.trajectory[fi]
        xyz[k] = kabsch_align(ca.positions.copy(), rc_all)
        al_hd = kabsch_align(hd_sel.positions.copy(), rc_hd)
        sp_hd += al_hd; sp2_hd += al_hd**2
        al_art = kabsch_align(art_sel.positions.copy(), rc_art)
        sp_art += al_art; sp2_art += al_art**2
    mp_hd = sp_hd / n
    rmsf_hd = np.sqrt((sp2_hd / n - mp_hd**2).sum(axis=1))
    mp_art = sp_art / n
    rmsf_art = np.sqrt((sp2_art / n - mp_art**2).sum(axis=1))
    hd_m = (resids >= HD[0]) & (resids <= HD[1])
    art_m = (resids >= ART[0]) & (resids <= ART[1])
    Xres = xyz.transpose(1, 0, 2)
    R = np.zeros((len(ca), len(ca)))
    means = Xres.mean(axis=1); ss = ((Xres - means[:, None, :])**2).sum(axis=(1, 2))
    for i in range(len(ca)):
        ci = Xres[i] - means[i]
        for j in range(i+1, len(ca)):
            cj = Xres[j] - means[j]
            R[i, j] = R[j, i] = (ci * cj).sum() / np.sqrt(ss[i] * ss[j])
    block = R[np.ix_(hd_m, art_m)]
    rmsf_full = np.full(len(ca), np.nan)
    rmsf_full[hd_m] = rmsf_hd
    rmsf_full[art_m] = rmsf_art
    return dict(name=name, n=n, hd=rmsf_hd.mean(), hd_sd=rmsf_hd.std(),
                art=rmsf_art.mean(), art_sd=rmsf_art.std(),
                dccm_abs=np.abs(block).mean(), dccm_mean=block.mean(),
                dccm_pos=(block > 0.3).mean() * 100,
                dccm_neg=(block < -0.3).mean() * 100,
                rmsf=rmsf_full, resids=resids, block=block,
                hd_ids=resids[hd_m], art_ids=resids[art_m])


print(f'{"system":<13}{"frames":>5}{"HD RMSF":>9}{"ART RMSF":>10}{"DCCM|r|":>9}{"r>0.3%":>8}{"r<-0.3%":>9}')
results = {}
for name, top, dcd in SYS:
    r = analyze(name, top, dcd)
    results[name] = r
    print(f'{name:<13}{r["n"]:>5}{r["hd"]:>8.2f}\u00b1{r["hd_sd"]:.2f}'
          f'{r["art"]:>9.2f}\u00b1{r["art_sd"]:.2f}{r["dccm_abs"]:>9.3f}'
          f'{r["dccm_pos"]:>8.1f}{r["dccm_neg"]:>9.1f}')

outdir = Path('results/analysis/rmsf_rerun3'); outdir.mkdir(exist_ok=True)
for name, r in results.items():
    np.savetxt(outdir / f'{name}_rmsf_residue.dat',
               np.column_stack([r['resids'], r['rmsf']]),
               fmt='%d %.4f', header='resid rmsf')
    np.save(outdir / f'{name}_dccm_block.npy', r['block'])
    print(f"{name}: signed_mean={r['dccm_mean']:.3f} signed_sd={r['block'].std():.3f}")
print('saved to', outdir)
