#!/usr/bin/env python3
"""从 pmf-c1 xvg 文件提取原面板 S2 CV1/CV2 井深 (cumulant 管线)"""
import numpy as np, glob, os

os.chdir('results/analysis')
for cv in ['CV1', 'CV2']:
    print(f'--- {cv} ---')
    for f in sorted(glob.glob(f'pmf-c1-sys2_*_{cv}_cv.dat.xvg')):
        name = f.replace('pmf-c1-sys2_', '').replace(f'_{cv}_cv.dat.xvg', '')
        d = np.loadtxt(f, comments=['#', '@'])
        if d.ndim == 1:
            d = d.reshape(-1, 2)
        rc, pmf = d[:, 0], d[:, 1]
        imin = pmf.argmin(); imax = pmf.argmax()
        print(f'{name:<12} rc范围 {rc.min():.1f}-{rc.max():.1f}  '
              f'井深 {pmf.max()-pmf.min():.2f}  (min@{rc[imin]:.1f}, max@{rc[imax]:.1f})')
