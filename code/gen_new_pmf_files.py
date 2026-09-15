#!/usr/bin/env python3
"""生成新药 pmf-c3 文件到图脚本所需位置:
- S1: results/analysis/sys1_<lig>_pmf_c3.xvg (DBE 权重, 200ns)
- S2: results/analysis/pmf-c3-sys2_<lig>_CV1/CV2_cv.dat.xvg (DFW 权重)
"""
import numpy as np
import subprocess, tempfile, shutil

PW = '/Users/taozhu/clacky_workspace/PARPi_design/tools/PyReweighting/PyReweighting-1D.py'
A = '/Users/taozhu/clacky_workspace/PARPi_design/results/analysis'
S1 = f'{A}/new_drugs'
S2 = f'{A}/new_drugs_s2'

def run_pw(cv, w, tag):
    tmp = tempfile.mkdtemp()
    np.savetxt(f'{tmp}/cv.dat', cv, fmt='%.4f')
    np.savetxt(f'{tmp}/weights.dat', w, fmt='%.6f')
    subprocess.run(['python3', PW, '-input', 'cv.dat', '-T', '300', '-disc', '0.1',
                    '-Emax', '20', '-cutoff', '2', '-job', 'amdweight_CE',
                    '-weight', 'weights.dat'], cwd=tmp, capture_output=True)
    return tmp

for d in ['fluzoparib', 'pamiparib', 'senaparib']:
    # S1 (DBE)
    cv = np.loadtxt(f'{S1}/sys1_{d}_cv.dat')
    w = np.loadtxt(f'{S1}/sys1_{d}_weights.dat')
    tmp = run_pw(cv, w, d)
    shutil.copy(f'{tmp}/pmf-c3-cv.dat.xvg', f'{A}/sys1_{d}_pmf_c3.xvg')
    # S2 CV1 + CV2 (DFW)
    w2 = np.loadtxt(f'{S2}/{d}/analysis_weights_dfw.dat')
    for cvn in ['CV1', 'CV2']:
        cv2 = np.loadtxt(f'{S2}/{d}/analysis_{cvn}.dat')
        tmp2 = run_pw(cv2, w2, f'{d}_{cvn}')
        shutil.copy(f'{tmp2}/pmf-c3-cv.dat.xvg', f'{A}/pmf-c3-sys2_{d}_{cvn}_cv.dat.xvg')
    print(f'{d}: S1 + S2 CV1/CV2 pmf-c3 已生成')
print('DONE')
