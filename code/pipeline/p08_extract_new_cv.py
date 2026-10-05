#!/usr/bin/env python3
"""p08_extract_new_cv.py — extract the HD-ART CV + DBE weights for the extension S1 systems

Purpose:   Chunk-streamed extraction of the HD-ART centre-of-mass distance CV
           and the DBE weights (gamd.log column 8) for the three extension
           inhibitors' S1 runs. Memory footprint < 500 MB (2 GB container limit).
Author:    Tao Zhu (tjogzt@gmail.com)
Created:   2026-09-16
Inputs:    runs/s1_<drug>/{s1_<drug>.prmtop, gamd_out/output.dcd, gamd_out/gamd.log}
Outputs:   runs/s1_<drug>/analysis_cv.dat + analysis_weights.dat
Depends:   mdtraj, numpy
Run:       python3 code/pipeline/p08_extract_new_cv.py   (on the GPU instance)
"""
import mdtraj as md
import numpy as np

beta = 1.0 / 0.596  # 300 K, matches the original protocol
drugs = ['fluzoparib', 'pamiparib', 'senaparib']

for d in drugs:
    base = f'runs/s1_{d}'
    prmtop = f'{base}/s1_{d}.prmtop'
    dcd = f'{base}/gamd_out/output.dcd'
    logpath = f'{base}/gamd_out/gamd.log'

    top = md.load_prmtop(prmtop)
    # New systems use 0-based residues (0-344 of the original 1-345):
    # HD = first 119 residues, ART = last 226 residues
    hd = top.select('resid 0 to 118 and name CA')
    art = top.select('resid 119 to 344 and name CA')
    print(f'[{d}] HD_CA={len(hd)} ART_CA={len(art)}', flush=True)

    cv_chunks = []
    for chunk in md.iterload(dcd, top=top, chunk=500):
        cv = np.linalg.norm(chunk.xyz[:, hd, :].mean(axis=1) - chunk.xyz[:, art, :].mean(axis=1), axis=1) * 10
        cv_chunks.append(cv)
    cv = np.concatenate(cv_chunks)
    print(f'[{d}] DCD frames={len(cv)}, CV mean={cv.mean():.2f}', flush=True)

    log = np.loadtxt(logpath)
    steps = log[:, 1]
    dbe = log[:, 7]
    print(f'[{d}] log rows={len(log)}, first step={steps[0]:.0f}, last step={steps[-1]:.0f}', flush=True)

    # Frame alignment: DCD frame i corresponds to log row i (both recorded
    # every 25000 steps, aligned from the start)
    n = min(len(cv), len(log))
    prod = steps[:n] >= 6000000
    cv_p = cv[:n][prod]
    dbe_p = dbe[:n][prod]
    print(f'[{d}] production frames={prod.sum()} CV mean={cv_p.mean():.2f} range={cv_p.min():.2f}-{cv_p.max():.2f} '
          f'DBE mean={dbe_p.mean():.2f}±{dbe_p.std():.2f}', flush=True)

    np.savetxt(f'{base}/analysis_cv.dat', cv_p, fmt='%.4f')
    np.savetxt(f'{base}/analysis_weights.dat',
               np.column_stack([dbe_p * beta, np.zeros(len(dbe_p)), dbe_p]), fmt='%.4f')
    print(f'[{d}] saved ({len(cv_p)} frames)', flush=True)

print('ALL DONE')
