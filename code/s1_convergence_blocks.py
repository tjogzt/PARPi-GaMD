#!/usr/bin/env python3
"""S1 分块收敛: 新药 well depth vs 累计时长 (每 20ns 一块, 验证稿件'80ns内收敛'的断言)"""
import numpy as np
import subprocess, tempfile

PW = '/Users/taozhu/clacky_workspace/PARPi_design/tools/PyReweighting/PyReweighting-1D.py'
DATA = '/Users/taozhu/clacky_workspace/PARPi_design/results/analysis/new_drugs'

def wd_from_pmf(f):
    d = []
    for line in open(f):
        line = line.strip()
        if line and not line.startswith(('#', '@')):
            p = line.split()
            d.append([float(p[0]), float(p[1])])
    pmf = np.array(d)[:, 1]
    return pmf.max() - pmf.min()

def run_pw(cv, w, tmp):
    np.savetxt(f'{tmp}/cv.dat', cv, fmt='%.4f')
    np.savetxt(f'{tmp}/weights.dat', w, fmt='%.4f')
    subprocess.run(['python3', PW, '-input', 'cv.dat', '-T', '300', '-disc', '0.1',
                    '-Emax', '20', '-cutoff', '2', '-job', 'amdweight_CE',
                    '-weight', 'weights.dat'],
                   cwd=tmp, capture_output=True)
    return wd_from_pmf(f'{tmp}/pmf-c3-cv.dat.xvg')

drugs = ['fluzoparib', 'pamiparib', 'senaparib']
print('=== 累计 well depth (C3) ===')
print(f'{"药":<12}' + ''.join([f'{t:>9}' for t in ['20','40','60','80','100','120','140','160','180','200ns']]))
for d in drugs:
    cv = np.loadtxt(f'{DATA}/sys1_{d}_cv.dat')
    w = np.loadtxt(f'{DATA}/sys1_{d}_weights.dat')
    row = []
    for t in range(20, 201, 20):
        n = t * 20  # 50ps/帧 → 20ns = 400帧
        tmp = tempfile.mkdtemp()
        row.append(f'{run_pw(cv[:n], w[:n], tmp):>9.1f}')
    print(f'{d:<12}' + ''.join(row))
