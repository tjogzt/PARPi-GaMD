#!/usr/bin/env python3
"""用补跑新数据重生成 AZD5305/veliparib 的 S2 pmf-c3 文件 (DFW 权重)"""
import numpy as np
import subprocess, tempfile, shutil

PW = '/Users/taozhu/clacky_workspace/PARPi_design/tools/PyReweighting/PyReweighting-1D.py'
A = '/Users/taozhu/clacky_workspace/PARPi_design/results/analysis'
S2 = f'{A}/new_drugs_s2'

for d in ['AZD5305', 'veliparib']:
    w = np.loadtxt(f'{S2}/{d}/analysis_weights_dfw.dat')
    for cvn in ['CV1', 'CV2']:
        cv = np.loadtxt(f'{S2}/{d}/analysis_{cvn}.dat')
        tmp = tempfile.mkdtemp()
        np.savetxt(f'{tmp}/cv.dat', cv, fmt='%.4f')
        np.savetxt(f'{tmp}/weights.dat', w, fmt='%.6f')
        subprocess.run(['python3', PW, '-input', 'cv.dat', '-T', '300', '-disc', '0.1',
                        '-Emax', '20', '-cutoff', '2', '-job', 'amdweight_CE',
                        '-weight', 'weights.dat'], cwd=tmp, capture_output=True)
        dst = f'{A}/pmf-c3-sys2_{d}_{cvn}_cv.dat.xvg'
        shutil.copy(f'{tmp}/pmf-c3-cv.dat.xvg', dst)
        print(f'{d} {cvn} → {dst} (覆盖 2.6ns 旧值)')
print('DONE')
