#!/usr/bin/env python3
"""决定性实验: 旧 S2 数据 × 两种权重方案, 复现稿件 C3 值 (talazoparib CV2 = 30.6)
(a) 旧方案: weight = exp(beta*DFW)  [旧权重文件]
(b) 标准方案: weight = exp(beta*DBE) [gamd.log col8]
"""
import numpy as np
import subprocess, tempfile

PW = '/Users/taozhu/clacky_workspace/PARPi_design/tools/PyReweighting/PyReweighting-1D.py'
CV = '/Users/taozhu/clacky_workspace/PARPi_design/results/analysis/sys2_talazoparib_CV2_cv.dat'
LOG = '/Volumes/tjogzt4T/PARPi_data/sys2_talazoparib/gamd.log'
beta = 1.677571

def wd(f):
    d = [l.split() for l in open(f) if l.strip() and not l.startswith(('#', '@'))]
    pmf = np.array([[float(x[0]), float(x[1])] for x in d])[:, 1]
    return pmf.max() - pmf.min()

def run_pw(cv, wcol1, dV, tmp):
    np.savetxt(f'{tmp}/cv.dat', cv, fmt='%.4f')
    np.savetxt(f'{tmp}/weights.dat', np.column_stack([wcol1, np.zeros(len(cv)), dV]), fmt='%.6f')
    subprocess.run(['python3', PW, '-input', 'cv.dat', '-T', '300', '-disc', '0.1',
                    '-Emax', '20', '-cutoff', '2', '-job', 'amdweight_CE',
                    '-weight', 'weights.dat'], cwd=tmp, capture_output=True)
    return [wd(f'{tmp}/pmf-c{i}-cv.dat.xvg') for i in (1, 2, 3)]

cv = np.loadtxt(CV)
log = np.loadtxt(LOG)
n = min(len(cv), len(log))
cv, log = cv[:n], log[:n]
prod = log[:, 0] == 1
cv_p, dfw, dbe = cv[prod], log[prod, 5], log[prod, 7]
print(f'生产帧: {prod.sum()}, DFW {dfw.mean():.3f}±{dfw.std():.3f}, DBE {dbe.mean():.2f}±{dbe.std():.2f}')

tmp = tempfile.mkdtemp()
wa = run_pw(cv_p, beta * dfw, dfw, tmp)     # 旧方案 col1=beta*DFW, dV=DFW
tmp2 = tempfile.mkdtemp()
wb = run_pw(cv_p, beta * dbe, dbe, tmp2)    # 标准方案 col1=beta*DBE, dV=DBE
print(f'(a) DFW 方案: C1={wa[0]:.1f} C2={wa[1]:.1f} C3={wa[2]:.1f}')
print(f'(b) DBE 方案: C1={wb[0]:.1f} C2={wb[1]:.1f} C3={wb[2]:.1f}')
print(f'稿件 talazoparib S2 CV2 C3 = 30.6')
