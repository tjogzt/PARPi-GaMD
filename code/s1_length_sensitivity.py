#!/usr/bin/env python3
"""S1 时长一致性分析: 把新药 S1 (200ns) 截取到与旧药一致的时长重算 well depth
旧药生产时长(帧×1ps): APO 21.5 / AZD5305 25 / olaparib 19 / niraparib 31.2 / rucaparib 31.2 / veliparib 30.8 ns
新药数据: cv.dat + weights.dat (4001 帧 × 50ps = 200ns)
截取点: 19 / 22 / 31 ns (覆盖旧药范围) + 全 200ns
"""
import numpy as np
import subprocess, os, tempfile

PW = '/Users/taozhu/clacky_workspace/PARPi_design/tools/PyReweighting/PyReweighting-1D.py'
DATA = '/Users/taozhu/clacky_workspace/PARPi_design/results/analysis/new_drugs'

def well_depth_from_pmf(pmf_file):
    data = []
    for line in open(pmf_file):
        line = line.strip()
        if line and not line.startswith(('#', '@')):
            p = line.split()
            data.append([float(p[0]), float(p[1])])
    pmf = np.array(data)[:, 1]
    return pmf.max() - pmf.min()

def run_pw(cv, weights, tag, tmpdir):
    """PyReweighting C1-C3, 返回 c3 well depth"""
    np.savetxt(f'{tmpdir}/cv.dat', cv, fmt='%.4f')
    np.savetxt(f'{tmpdir}/weights.dat', weights, fmt='%.4f')
    r = subprocess.run(
        ['python3', PW, '-input', 'cv.dat', '-T', '300', '-disc', '0.1',
         '-Emax', '20', '-cutoff', '2', '-job', 'amdweight_CE', '-weight', 'weights.dat'],
        cwd=tmpdir, capture_output=True, text=True)
    c3min = None
    for line in r.stdout.splitlines():
        if 'pmf_min-c3' in line:
            c3min = float(line.split('=')[1])
    wd = well_depth_from_pmf(f'{tmpdir}/pmf-c3-cv.dat.xvg')
    return wd, c3min

drugs = ['fluzoparib', 'pamiparib', 'senaparib']
# 截取帧数 (50ps/帧): 19ns=380帧, 22ns=440帧, 31ns=620帧
cutoffs = {'19ns': 380, '22ns': 440, '31ns': 620}

print('=== 新药 S1 well depth: 时长截断敏感性 ===')
print(f'{"药":<12}{"200ns(全)":>10}' + ''.join([f'{k:>10}' for k in cutoffs]))
for d in drugs:
    cv = np.loadtxt(f'{DATA}/sys1_{d}_cv.dat')
    w = np.loadtxt(f'{DATA}/sys1_{d}_weights.dat')
    full_wd, _ = run_pw(cv, w, d, f'{DATA}/{d}')  # 全 200ns (已有结果复核)
    results = [f'{full_wd:>10.1f}']
    for k, n in cutoffs.items():
        tmp = tempfile.mkdtemp()
        wd, _ = run_pw(cv[:n], w[:n], f'{d}_{k}', tmp)
        results.append(f'{wd:>10.1f}')
    print(f'{d:<12}' + ''.join(results))
print()
print('参考(旧药 C3): veliparib 101.1 / olaparib 68.9 / rucaparib 53.4 / niraparib 51.6 / APO 42.5 / AZD5305 28.2')
print('(旧药时长 19-31ns; 新药截断后与旧药同协议可比)')
