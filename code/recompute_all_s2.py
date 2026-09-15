#!/usr/bin/env python3
"""Unified-protocol recomputation of all 7 S2 systems: RMSF (HD/ART) + DCCM (HDxART)
Protocol: production segment, 50 ps frame rate, whole-protein CA Kabsch alignment
of the first frame, per-residue RMSF; DCCM = residue-level CA displacement-vector
Pearson correlation.
Legacy systems (1 ps frame rate): stride 50; rebuilt systems (50 ps): last 517
production frames.
"""
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import warnings
import sys
warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.kabsch import kabsch_align
from common.paths import data_root

DATA = data_root()
NEW = data_root() / 's2_new_drugs'
HD = (662, 787); ART = (788, 1014); CAT = (531, 1011)  # CAT-domain alignment (Methods)

# Legacy systems: (name, prmtop, dcd, stride)
OLD = [
    ('APO', DATA/'sys2_APO'/'sys2_APO.prmtop', DATA/'sys2_APO'/'output.dcd', 50),
    ('talazoparib', DATA/'sys2_talazoparib'/'sys2_talazoparib.prmtop', DATA/'sys2_talazoparib'/'output.dcd', 50),
    ('niraparib', DATA/'sys2_niraparib'/'sys2_niraparib.prmtop', DATA/'sys2_niraparib'/'output.dcd', 50),
    ('olaparib', DATA/'sys2_olaparib'/'sys2_olaparib.prmtop', DATA/'sys2_olaparib'/'output.dcd', 50),
    ('rucaparib', DATA/'sys2_rucaparib'/'sys2_rucaparib.prmtop', DATA/'sys2_rucaparib'/'output.dcd', 50),
]
NEWSYS = [('AZD5305', NEW/'sys2_AZD5305'/'AZD5305_best.prmtop', NEW/'sys2_AZD5305'/'output.dcd'),
          ('veliparib', NEW/'sys2_veliparib'/'veliparib_best.prmtop', NEW/'sys2_veliparib'/'output.dcd')]

def analyze(name, top, dcd, stride=None, n_prod=None, outdir=Path('results/analysis/rmsf_recomp')):
    u = mda.Universe(str(top), str(dcd))
    ca = u.select_atoms('protein and name CA')
    hd_sel = u.select_atoms(f'resid {HD[0]}:{HD[1]} and name CA')
    art_sel = u.select_atoms(f'resid {ART[0]}:{ART[1]} and name CA')
    resids = ca.residues.resids
    ntot = len(u.trajectory)
    if n_prod is None:
        frames = list(range(0, ntot, stride))
    else:
        frames = list(range(ntot - n_prod, ntot))
    n = len(frames)
    # Reference: whole protein / HD / ART (each first frame)
    u.trajectory[frames[0]]
    ref_all = ca.positions.copy(); rc_all = ref_all - ref_all.mean(axis=0)
    ref_hd = hd_sel.positions.copy(); rc_hd = ref_hd - ref_hd.mean(axis=0)
    ref_art = art_sel.positions.copy(); rc_art = ref_art - ref_art.mean(axis=0)

    sp_hd = np.zeros((len(hd_sel), 3)); sp2_hd = np.zeros_like(sp_hd)
    sp_art = np.zeros((len(art_sel), 3)); sp2_art = np.zeros_like(sp_art)
    xyz = np.zeros((n, len(ca), 3))   # whole-protein superposition (DCCM)
    for k, fi in enumerate(frames):
        u.trajectory[fi]
        # Whole-protein superposition -> xyz (DCCM)
        xyz[k] = kabsch_align(ca.positions.copy(), rc_all)
        # HD self-alignment -> HD RMSF
        al_hd = kabsch_align(hd_sel.positions.copy(), rc_hd)
        sp_hd += al_hd; sp2_hd += al_hd**2
        # ART self-alignment -> ART RMSF
        al_art = kabsch_align(art_sel.positions.copy(), rc_art)
        sp_art += al_art; sp2_art += al_art**2
    mp_hd = sp_hd / n
    rmsf_hd = np.sqrt((sp2_hd / n - mp_hd**2).sum(axis=1))
    mp_art = sp_art / n
    rmsf_art = np.sqrt((sp2_art / n - mp_art**2).sum(axis=1))

    hd_m = (resids >= HD[0]) & (resids <= HD[1])
    art_m = (resids >= ART[0]) & (resids <= ART[1])
    # DCCM (residue-level displacement-vector Pearson)
    Xres = xyz.transpose(1, 0, 2)
    R = np.zeros((len(ca), len(ca)))
    means = Xres.mean(axis=1); ss = ((Xres - means[:, None, :])**2).sum(axis=(1, 2))
    for i in range(len(ca)):
        ci = Xres[i] - means[i]
        for j in range(i+1, len(ca)):
            cj = Xres[j] - means[j]
            R[i, j] = R[j, i] = (ci * cj).sum() / np.sqrt(ss[i] * ss[j])
    block = R[np.ix_(hd_m, art_m)]
    # per-residue RMSF: HD segment from HD self-alignment, ART segment from
    # ART self-alignment, others NaN
    rmsf_full = np.full(len(ca), np.nan)
    rmsf_full[hd_m] = rmsf_hd
    rmsf_full[art_m] = rmsf_art
    return dict(name=name, n=n, hd=rmsf_hd.mean(), hd_sd=rmsf_hd.std(),
                art=rmsf_art.mean(), art_sd=rmsf_art.std(),
                dccm_abs=np.abs(block).mean(), dccm_mean=block.mean(),
                dccm_pos=(block > 0.3).mean() * 100,
                dccm_neg=(block < -0.3).mean() * 100, rmsf=rmsf_full, resids=resids,
                block=block, hd_ids=resids[hd_m], art_ids=resids[art_m])

print(f'{"system":<13}{"frames":>5}{"HD RMSF":>9}{"ART RMSF":>10}{"DCCM|r|":>9}{"r>0.3%":>8}{"r<-0.3%":>9}')
results = {}
for name, top, dcd, stride in OLD:
    r = analyze(name, top, dcd, stride=stride)
    results[name] = r
    print(f'{name:<13}{r["n"]:>5}{r["hd"]:>8.2f}±{r["hd_sd"]:.2f}{r["art"]:>9.2f}±{r["art_sd"]:.2f}{r["dccm_abs"]:>9.3f}{r["dccm_pos"]:>8.1f}{r["dccm_neg"]:>9.1f}')
for name, top, dcd in NEWSYS:
    r = analyze(name, top, dcd, n_prod=517)
    results[name] = r
    print(f'{name:<13}{r["n"]:>5}{r["hd"]:>8.2f}±{r["hd_sd"]:.2f}{r["art"]:>9.2f}±{r["art_sd"]:.2f}{r["dccm_abs"]:>9.3f}{r["dccm_pos"]:>8.1f}{r["dccm_neg"]:>9.1f}')

# Save per-residue data (for dRMSF recompute + RMSF figure regeneration)
outdir = Path('results/analysis/rmsf_recomp'); outdir.mkdir(exist_ok=True)
import csv
for name, r in results.items():
    np.savetxt(outdir / f'{name}_rmsf_residue.dat',
               np.column_stack([r['resids'], r['rmsf']]), fmt='%d %.4f',
               header=f'resid rmsf(A)  [{name}, n={r["n"]} frames]')
    # DCCM HDxART block (rows = HD, cols = ART)
    np.savetxt(outdir / f'{name}_dccm_block.csv', r['block'],
               delimiter=',', fmt='%.5f',
               header=','.join(map(str, r['art_ids'])), comments='# ')  # rows = HD resid order
with open(outdir / 's2_rmsf_dccm_uniform.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['system', 'n_frames', 'HD_rmsf', 'HD_sd', 'ART_rmsf', 'ART_sd',
                'dccm_abs', 'dccm_mean', 'dccm_pos_pct', 'dccm_neg_pct'])
    for name, r in results.items():
        w.writerow([name, r['n'], f"{r['hd']:.3f}", f"{r['hd_sd']:.3f}",
                    f"{r['art']:.3f}", f"{r['art_sd']:.3f}",
                    f"{r['dccm_abs']:.4f}", f"{r['dccm_mean']:.4f}",
                    f"{r['dccm_pos']:.2f}", f"{r['dccm_neg']:.2f}"])
    print(f'\nwrote {outdir/"s2_rmsf_dccm_uniform.csv"}')
