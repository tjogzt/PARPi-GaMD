#!/usr/bin/env python3
"""从旧 S2 DCD 提取最后一帧坐标 → .npy (本地, 供云端 parmed 拼装 PDB)"""
import struct, sys, os
import numpy as np

def parse_dcd_last(path):
    with open(path, 'rb') as f:
        hdr = f.read(84)
        n_frames = struct.unpack('<i', hdr[8:12])[0]
        pos = 84
        sz = struct.unpack('<i', f.read(4))[0]; f.read(sz)
        assert struct.unpack('<i', f.read(4))[0] == sz
        natoms = struct.unpack('<i', f.read(4))[0]
        assert struct.unpack('<i', f.read(4))[0] == natoms
        fsz = struct.unpack('<i', f.read(4))[0]
        assert fsz == natoms * 12, f'frame {fsz} != natoms*12 {natoms*12}'
        filesize = os.path.getsize(path)
        f.seek(filesize - fsz - 4)
        last = f.read(fsz)
    xyz = np.frombuffer(last, dtype=np.float32).reshape(natoms, 3)
    return n_frames, natoms, xyz

for drug, dcd, out in [
    ('AZD5305', '/Volumes/tjogzt4T/PARPi_data/sys2_AZD5305/output.dcd', '/tmp/azd5305_last_xyz.npy'),
    ('veliparib', '/Volumes/tjogzt4T/PARPi_data/sys2_veliparib/output.dcd', '/tmp/veliparib_last_xyz.npy'),
]:
    n, natoms, xyz = parse_dcd_last(dcd)
    np.save(out, xyz)
    print(f'{drug}: {n} 帧, 末帧 {natoms} 原子 → {out} ({os.path.getsize(out)} bytes)')
