#!/usr/bin/env python3
"""新药 S2: CV1/CV2 well depth — DFW 权重 (与稿件旧 S2 协议一致, 已复现 talazoparib 30.6)"""
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
    pmf = np.array(d)[:, 1]
    return pmf.max() - pmf.min()

def run_pw(cv, w, tmp):
    np.savetxt(f'{tmp}/cv.dat', cv, fmt='%.4f')
    np.savetxt(f'{tmp}/weights.dat', w, fmt='%.6f')
    r = subprocess.run(['python3', PW, '-input', 'cv.dat', '-T', '300', '-disc', '0.1',
                        '-Emax', '20', '-cutoff', '2', '-job', 'amdweight_CE',
                        '-weight', 'weights.dat'],
                       cwd=tmp, capture_output=True, text=True)
    if not __import__('os').path.exists(f'{tmp}/pmf-c3-cv.dat.xvg'):
        raise RuntimeError(f'PyReweighting 失败: {r.stdout[-300:]} {r.stderr[-300:]}')
    return [wd_from_pmf(f'{tmp}/pmf-c{i}-cv.dat.xvg') for i in (1, 2, 3)]

drugs = ['fluzoparib', 'pamiparib', 'senaparib']
print(f'{"药":<12}{"CV":<6}{"C1":>8}{"C2":>8}{"C3":>8}{"均值±SD":>12}')
for d in drugs:
    for cv_name in ['CV1', 'CV2']:
        cv = np.loadtxt(f'{DATA}/{d}/analysis_{cv_name}.dat')
        w = np.loadtxt(f'{DATA}/{d}/analysis_weights_dfw.dat')
        assert len(cv) == len(w), f'{d} {cv_name}: {len(cv)} vs {len(w)}'
        tmp = tempfile.mkdtemp()
        c1, c2, c3 = run_pw(cv, w, tmp)
        mean, sd = np.mean([c1, c2, c3]), np.std([c1, c2, c3])
        print(f'{d:<12}{cv_name:<6}{c1:>8.1f}{c2:>8.1f}{c3:>8.1f}{mean:>8.1f}±{sd:.1f}')
