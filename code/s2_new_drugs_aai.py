#!/usr/bin/env python3
"""新药 AAI = S1 well depth (C3, DBE 权重, 200ns) / S2 CV2 (C3, DFW 权重)"""
import numpy as np
import subprocess, tempfile

PW = '/Users/taozhu/clacky_workspace/PARPi_design/tools/PyReweighting/PyReweighting-1D.py'
S1DIR = '/Users/taozhu/clacky_workspace/PARPi_design/results/analysis/new_drugs'
S2DIR = '/Users/taozhu/clacky_workspace/PARPi_design/results/analysis/new_drugs_s2'

def wd_from_pmf(f):
    d = []
    for line in open(f):
        line = line.strip()
        if line and not line.startswith(('#', '@')):
            p = line.split()
            d.append([float(p[0]), float(p[1])])
    return np.array(d)[:, 1].max() - np.array(d)[:, 1].min()

def s1_cumulants(drug):
    """S1 C1/C2/C3 (DBE 权重, 全 200ns)"""
    cv = np.loadtxt(f'{S1DIR}/sys1_{drug}_cv.dat')
    w = np.loadtxt(f'{S1DIR}/sys1_{drug}_weights.dat')
    tmp = tempfile.mkdtemp()
    np.savetxt(f'{tmp}/cv.dat', cv, fmt='%.4f')
    np.savetxt(f'{tmp}/weights.dat', w, fmt='%.4f')
    subprocess.run(['python3', PW, '-input', 'cv.dat', '-T', '300', '-disc', '0.1',
                    '-Emax', '20', '-cutoff', '2', '-job', 'amdweight_CE',
                    '-weight', 'weights.dat'], cwd=tmp, capture_output=True)
    return [wd_from_pmf(f'{tmp}/pmf-c{i}-cv.dat.xvg') for i in (1, 2, 3)]

def s2_c3(drug):
    cv = np.loadtxt(f'{S2DIR}/{drug}/analysis_CV2.dat')
    w = np.loadtxt(f'{S2DIR}/{drug}/analysis_weights_dfw.dat')
    tmp = tempfile.mkdtemp()
    np.savetxt(f'{tmp}/cv.dat', cv, fmt='%.4f')
    np.savetxt(f'{tmp}/weights.dat', w, fmt='%.6f')
    subprocess.run(['python3', PW, '-input', 'cv.dat', '-T', '300', '-disc', '0.1',
                    '-Emax', '20', '-cutoff', '2', '-job', 'amdweight_CE',
                    '-weight', 'weights.dat'], cwd=tmp, capture_output=True)
    return wd_from_pmf(f'{tmp}/pmf-c3-cv.dat.xvg')

print(f'{"药":<12}{"S1 C1":>8}{"S1 C2":>8}{"S1 C3":>8}{"S1 均值±SD":>12}{"S2 CV2":>8}{"AAI":>8}{"AAI范围":>12}')
for d in ['fluzoparib', 'pamiparib', 'senaparib']:
    c1, c2, c3 = s1_cumulants(d)
    s2 = s2_c3(d)
    mean, sd = np.mean([c1, c2, c3]), np.std([c1, c2, c3])
    aai = mean / s2
    aai_lo, aai_hi = (mean - sd) / s2, (mean + sd) / s2
    print(f'{d:<12}{c1:>8.1f}{c2:>8.1f}{c3:>8.1f}{mean:>8.1f}±{sd:.1f}{s2:>8.1f}{aai:>8.2f}{aai_lo:>7.2f}-{aai_hi:.2f}')

print('\n旧体系 AAI (稿件表): tala 0.99 / AZD5305 0.96 / nira 1.68 / ruca 1.75 / ola 2.27 / veli 3.47 / APO 1.39')
