#!/usr/bin/env python3
"""AZD5305/veliparib 补跑 S2: DFW 权重 well depth (C1-C3)"""
import numpy as np
import subprocess, tempfile

PW = '/Users/taozhu/clacky_workspace/PARPi_design/tools/PyReweighting/PyReweighting-1D.py'
DATA = '/Users/taozhu/clacky_workspace/PARPi_design/results/analysis/new_drugs_s2'

def wd_from_pmf(f):
    d = []
    for line in open(f):
        line = line.strip()
        if line and not line.startswith(('#', '@')):
            p = line.split()
            d.append([float(p[0]), float(p[1])])
    return np.array(d)[:, 1].max() - np.array(d)[:, 1].min()

def run_pw(cv, w, tmp):
    np.savetxt(f'{tmp}/cv.dat', cv, fmt='%.4f')
    np.savetxt(f'{tmp}/weights.dat', w, fmt='%.6f')
    subprocess.run(['python3', PW, '-input', 'cv.dat', '-T', '300', '-disc', '0.1',
                    '-Emax', '20', '-cutoff', '2', '-job', 'amdweight_CE',
                    '-weight', 'weights.dat'], cwd=tmp, capture_output=True)
    return [wd_from_pmf(f'{tmp}/pmf-c{i}-cv.dat.xvg') for i in (1, 2, 3)]

print(f'{"药":<12}{"CV":<6}{"C1":>8}{"C2":>8}{"C3":>8}')
for d in ['AZD5305', 'veliparib']:
    for cvn in ['CV1', 'CV2']:
        cv = np.loadtxt(f'{DATA}/{d}/analysis_{cvn}.dat')
        w = np.loadtxt(f'{DATA}/{d}/analysis_weights_dfw.dat')
        assert len(cv) == len(w), f'{d} {cvn} 长度不一致'
        tmp = tempfile.mkdtemp()
        c1, c2, c3 = run_pw(cv, w, tmp)
        print(f'{d:<12}{cvn:<6}{c1:>8.1f}{c2:>8.1f}{c3:>8.1f}')
