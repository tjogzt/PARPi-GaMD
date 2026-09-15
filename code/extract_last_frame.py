#!/usr/bin/env python3
"""Extract last-frame coordinates from a legacy S2 DCD to .npy (local; for cloud parmed PDB assembly)"""
import struct, sys, os
import numpy as np
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root

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
    ('AZD5305', str(data_root() / 'sys2_AZD5305' / 'output.dcd'), '/tmp/azd5305_last_xyz.npy'),
    ('veliparib', str(data_root() / 'sys2_veliparib' / 'output.dcd'), '/tmp/veliparib_last_xyz.npy'),
]:
    n, natoms, xyz = parse_dcd_last(dcd)
    np.save(out, xyz)
    print(f'{drug}: {n} frames, last frame {natoms} atoms -> {out} ({os.path.getsize(out)} bytes)')
