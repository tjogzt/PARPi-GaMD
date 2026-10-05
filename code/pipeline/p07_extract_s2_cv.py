#!/usr/bin/env python3
"""p07_extract_s2_cv.py — extract CV series + DBE weights for the S2 systems

Purpose:   Chunk-streamed extraction of the two collective variables for the
           full-length S2 systems:
             CV1 = protein CA centre-of-mass - DNA heavy-atom COM distance
             CV2 = HD(662-787) CA COM - ART(788-1011) CA COM distance
           Weights come from gamd.log column 8 (DBE); the production segment
           starts at step >= 6.1M (cMD 1M + prep 100K + equil 5M).
Author:    Tao Zhu (tjogzt@gmail.com)
Created:   2026-09-15
Inputs:    <runs_dir>/sys2_<drug>/{build/<drug>_best.prmtop, gamd_out/output.dcd,
           gamd_out/gamd.log}
Outputs:   analysis_CV1.dat / analysis_CV2.dat / analysis_weights.dat (per system)
Depends:   mdtraj, numpy
Run:       python3 code/pipeline/p07_extract_s2_cv.py <runs_dir> <drug> [<drug> ...]
"""
import mdtraj as md
import numpy as np

beta = 1.0 / 0.596  # 300 K

def main(runs_dir, drug, prod_step=6100000, prod_frames=None):
    base = f'{runs_dir}/sys2_{drug}'
    prmtop = f'{base}/build/{drug}_best.prmtop'
    dcd = f'{base}/gamd_out/output.dcd'
    logpath = f'{base}/gamd_out/gamd.log'

    top = md.load_prmtop(prmtop)
    hd = top.select('resid 662 to 787 and name CA')
    art = top.select('resid 788 to 1011 and name CA')
    prot = top.select('protein and name CA')
    dna = top.select('(resname DA or resname DT or resname DG or resname DC or '
                     'resname DG5 or resname DC3) and (not element H)')
    assert len(hd) == 126, f'HD CA={len(hd)} (expected 126)'
    assert len(art) == 224, f'ART CA={len(art)} (expected 224)'
    assert len(dna) > 900, f'DNA heavy atoms={len(dna)} (unexpected)'
    print(f'[{drug}] HD={len(hd)} ART={len(art)} ProtCA={len(prot)} DNAheavy={len(dna)}', flush=True)

    cv1_chunks, cv2_chunks = [], []
    for chunk in md.iterload(dcd, top=top, chunk=200):
        pc = chunk.xyz[:, prot, :].mean(axis=1)
        dc = chunk.xyz[:, dna, :].mean(axis=1)
        cv1 = np.linalg.norm(pc - dc, axis=1) * 10
        cv2 = np.linalg.norm(chunk.xyz[:, hd, :].mean(axis=1) - chunk.xyz[:, art, :].mean(axis=1), axis=1) * 10
        cv1_chunks.append(cv1); cv2_chunks.append(cv2)
    cv1 = np.concatenate(cv1_chunks); cv2 = np.concatenate(cv2_chunks)
    print(f'[{drug}] DCD={len(cv1)} frames, CV1 mean={cv1.mean():.1f}, CV2 mean={cv2.mean():.1f}', flush=True)

    log = np.loadtxt(logpath)
    steps, mode, dbe = log[:, 1], log[:, 0], log[:, 7]
    n = min(len(cv1), len(log))
    prod = (steps[:n] >= prod_step) & (mode[:n] == 1)
    if prod_frames:
        prod = prod & (np.arange(n) < prod_frames + int(prod.sum() == 0))
    print(f'[{drug}] production frames={prod.sum()} (of {n})', flush=True)
    assert prod.sum() > 200, f'too few production frames: {prod.sum()}'

    cv1p, cv2p, dbep = cv1[:n][prod], cv2[:n][prod], dbe[:n][prod]
    np.savetxt(f'{base}/analysis_CV1.dat', cv1p, fmt='%.4f')
    np.savetxt(f'{base}/analysis_CV2.dat', cv2p, fmt='%.4f')
    np.savetxt(f'{base}/analysis_weights.dat',
               np.column_stack([dbep * beta, np.zeros(len(dbep)), dbep]), fmt='%.4f')
    print(f'[{drug}] saved: CV1={cv1p.mean():.2f}±{cv1p.std():.2f} '
          f'CV2={cv2p.mean():.2f}±{cv2p.std():.2f} DBE={dbep.mean():.2f}±{dbep.std():.2f}', flush=True)

if __name__ == '__main__':
    import sys
    runs_dir = sys.argv[1] if len(sys.argv) > 1 else 'runs_s2'
    for d in sys.argv[2:]:
        main(runs_dir, d)
    print('ALL DONE')
