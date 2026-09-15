#!/usr/bin/env python3
"""提取三个新药 S1 的 HD-ART COM 距离 CV + DBE 权重 (chunk 流式, 内存<500MB, 适配2GB容器限制)"""
import mdtraj as md
import numpy as np

beta = 1.0 / 0.596  # 300K, 与原协议一致
drugs = ['fluzoparib', 'pamiparib', 'senaparib']

for d in drugs:
    base = f'runs/s1_{d}'
    prmtop = f'{base}/s1_{d}.prmtop'
    dcd = f'{base}/gamd_out/output.dcd'
    logpath = f'{base}/gamd_out/gamd.log'

    top = md.load_prmtop(prmtop)
    # 新体系残基 0 索引(0-344, 对应原 1-345): HD=前119残基, ART=后226残基
    hd = top.select('resid 0 to 118 and name CA')
    art = top.select('resid 119 to 344 and name CA')
    print(f'[{d}] HD_CA={len(hd)} ART_CA={len(art)}', flush=True)

    cv_chunks = []
    for chunk in md.iterload(dcd, top=top, chunk=500):
        cv = np.linalg.norm(chunk.xyz[:, hd, :].mean(axis=1) - chunk.xyz[:, art, :].mean(axis=1), axis=1) * 10
        cv_chunks.append(cv)
    cv = np.concatenate(cv_chunks)
    print(f'[{d}] DCD帧数={len(cv)} CV mean={cv.mean():.2f}', flush=True)

    log = np.loadtxt(logpath)
    steps = log[:, 1]
    dbe = log[:, 7]
    print(f'[{d}] log行数={len(log)} 首step={steps[0]:.0f} 尾step={steps[-1]:.0f}', flush=True)

    # 帧对齐: DCD 帧 i 对应 log 行 i (两者均每 25000 步记录, 从头对齐)
    n = min(len(cv), len(log))
    prod = steps[:n] >= 6000000
    cv_p = cv[:n][prod]
    dbe_p = dbe[:n][prod]
    print(f'[{d}] 生产帧={prod.sum()} CV mean={cv_p.mean():.2f} range={cv_p.min():.2f}-{cv_p.max():.2f} '
          f'DBE mean={dbe_p.mean():.2f}±{dbe_p.std():.2f}', flush=True)

    np.savetxt(f'{base}/analysis_cv.dat', cv_p, fmt='%.4f')
    np.savetxt(f'{base}/analysis_weights.dat',
               np.column_stack([dbe_p * beta, np.zeros(len(dbe_p)), dbe_p]), fmt='%.4f')
    print(f'[{d}] 保存完成 ({len(cv_p)} 帧)', flush=True)

print('ALL DONE')
