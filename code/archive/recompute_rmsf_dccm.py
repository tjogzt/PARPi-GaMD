#!/usr/bin/env python3
"""用补跑 26ns DCD 重算 AZD5305/veliparib 的 S2 RMSF + DCCM
(与原 19-rmsf_all_systems.py / 20-dccm_network.R 协议一致:
 RMSF = Kabsch 对齐首帧后 per-residue; HD=662-787, ART=788-1014
 DCCM = CA 坐标 Pearson 相关 (bio3d dccm 等价), HD×ART 块 mean|r| + 分数)
"""
import MDAnalysis as mda
import numpy as np
from pathlib import Path

DATA = Path('/Volumes/tjogzt4T/PARPi_data/s2_new_drugs')
HD = (662, 787)
ART = (788, 1014)

def kabsch_align(mobile, ref):
    ref_c = ref - ref.mean(axis=0)
    mob_c = mobile - mobile.mean(axis=0)
    H = mob_c.T @ ref_c
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return mob_c @ R

for name in ['AZD5305', 'veliparib']:
    top = DATA / f'sys2_{name}' / f'{name}_best.prmtop'
    dcd = DATA / f'sys2_{name}' / 'output.dcd'
    u = mda.Universe(str(top), str(dcd))
    ca = u.select_atoms('protein and name CA')
    resids = ca.residues.resids
    # 生产帧 = 最后 517 帧 (gamd.log mode=1, step>=6.1M; cMD+equil 在前 243 帧)
    n_prod = 517
    n_frames = len(u.trajectory)
    start = n_frames - n_prod
    print(f'=== {name}: {n_frames} 帧(含cMD+equil), 用生产段后 {n_prod} 帧, {len(ca)} CA ===')

    # --- RMSF (对齐生产段首帧) ---
    u.trajectory[start]
    ref = ca.positions.copy()
    sum_p = np.zeros((len(ca), 3)); sum_p2 = np.zeros((len(ca), 3))
    for ts in u.trajectory[start:]:
        al = kabsch_align(ca.positions.copy(), ref)
        sum_p += al; sum_p2 += al**2
    mean_p = sum_p / n_prod
    rmsf = np.sqrt((sum_p2 / n_prod - mean_p**2).sum(axis=1))
    hd_mask = (resids >= HD[0]) & (resids <= HD[1])
    art_mask = (resids >= ART[0]) & (resids <= ART[1])
    hd_mean = rmsf[hd_mask].mean(); hd_sd = rmsf[hd_mask].std()
    art_mean = rmsf[art_mask].mean(); art_sd = rmsf[art_mask].std()
    print(f'  HD RMSF = {hd_mean:.2f} ± {hd_sd:.2f} Å  ({hd_mask.sum()} 残基)')
    print(f'  ART RMSF = {art_mean:.2f} ± {art_sd:.2f} Å  ({art_mask.sum()} 残基)')

    # --- DCCM (生产段原始 CA 坐标 Pearson, 残基级向量相关) ---
    xyz = np.zeros((n_prod, len(ca), 3))
    for i, ts in enumerate(u.trajectory[start:]):
        xyz[i] = ca.positions
    Xres = xyz.transpose(1, 0, 2)  # 残基 × 帧 × 3
    R = np.zeros((len(ca), len(ca)))
    for i in range(len(ca)):
        Xi = Xres[i]
        ci = Xi - Xi.mean(axis=0)
        ss_i = (ci**2).sum()
        for j in range(i+1, len(ca)):
            cj = Xres[j] - Xres[j].mean(axis=0)
            ss_j = (cj**2).sum()
            R[i, j] = R[j, i] = (ci * cj).sum() / np.sqrt(ss_i * ss_j)
    block = R[np.ix_(hd_mask, art_mask)]
    print(f'  DCCM HD×ART: mean|r| = {np.abs(block).mean():.3f}, '
          f'r>0.3 = {(block>0.3).mean()*100:.1f}%, r<-0.3 = {(block<-0.3).mean()*100:.1f}%')
