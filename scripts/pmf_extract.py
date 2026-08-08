#!/usr/bin/env python3
"""Extract PMF from completed GaMD trajectory (v3 - half-split CA)."""
import mdtraj as md, numpy as np, sys, os

tag = sys.argv[1]
rundir = f'/root/autodl-tmp/PARPi_design/runs/{tag}'
outdir = f'{rundir}/analysis'; os.makedirs(outdir, exist_ok=True)

print(f'=== {tag} ===')
prmtop = f'{rundir}/{tag}.prmtop'
dcd = f'{rundir}/gamd_out/output.dcd'
traj = md.load(dcd, top=prmtop)
print(f'Frames: {len(traj)}, atoms: {traj.n_atoms}')

ca = traj.topology.select('name CA')
hd_ca = ca[:len(ca)//2]
art_ca = ca[len(ca)//2:]
print(f'HD CA: {len(hd_ca)}, ART CA: {len(art_ca)}')

hd_com = traj.xyz[:, hd_ca, :].mean(axis=1)
art_com = traj.xyz[:, art_ca, :].mean(axis=1)
dist = np.linalg.norm(hd_com - art_com, axis=1) * 10

weights_data = np.loadtxt(f'{rundir}/gamd_out/gamd.log')
prod_mask = weights_data[:, 0] == 1
boost = weights_data[prod_mask, -2]
n = min(len(dist), len(boost))
hist, edges = np.histogram(dist[:n], bins=50, weights=np.exp(boost[:n]), density=True)
pmf = -np.log(hist + 1e-10); pmf -= pmf.min()
rc = (edges[:-1] + edges[1:]) / 2

np.save(f'{outdir}/rc.npy', rc)
np.save(f'{outdir}/pmf.npy', pmf)
np.save(f'{outdir}/dist_raw.npy', dist)

with open(f'{outdir}/HD_ART_dist_pmf.xvg', 'w') as f:
    f.write('@    title "HD-ART COM Distance PMF"\n')
    f.write('@    xaxis  label "Distance (A)"\n')
    f.write('@    yaxis  label "PMF (kcal/mol)"\n')
    for x, y in zip(rc, pmf):
        f.write(f'{x:.3f} {y:.3f}\n')

print(f'Well depth: {pmf.max():.1f} kcal/mol')
print('DONE')
