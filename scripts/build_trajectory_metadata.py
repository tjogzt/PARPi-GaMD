#!/usr/bin/env python3
"""Build trajectory metadata table (reviewer Section 8 minimum verification package).

For every system directory under DATA_ROOT: topology hash, trajectory frames,
frame interval, boost log fields, restart history, timestamps. Emits
results/analysis/trajectory_metadata.csv.
"""
import hashlib
import os
import sys
import csv
from pathlib import Path

DATA_ROOT = Path(os.environ.get("DATA_ROOT", "/Volumes/tjogzt4T/PARPi_data"))
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis" / "trajectory_metadata.csv"

import MDAnalysis as mda  # noqa: E402

def md5(p, chunk=1 << 20):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()[:12]

rows = []
for d in sorted(DATA_ROOT.iterdir()):
    if not d.is_dir():
        continue
    prmtops = list(d.glob("*.prmtop"))
    if not prmtops:
        continue
    prmtop = prmtops[0]
    dcdfs = list(d.glob("*.dcd")) + list(d.glob("*.nc"))
    logs = list(d.glob("gamd.log"))
    system_id = d.name
    compound = system_id.split("_", 1)[1] if "_" in system_id else "?"
    rec = {
        "system_id": system_id,
        "compound_id": compound,
        "construct_id": "S1" if system_id.startswith("sys1") else "S2",
        "topology": prmtop.name,
        "topology_hash": md5(prmtop),
        "trajectory": dcdfs[0].name if dcdfs else "NONE",
        "frame_interval_ps": "",
        "frames": "",
        "production_ns": "",
        "boost_log": logs[0].name if logs else "NONE",
        "log_fields": "",
        "restart_history": "",
        "random_seed": "unknown",
        "software_commit": "see release env file",
        "analysis_config": "see data_manifest.md",
    }
    # restart-history heuristics from file names
    restr = [f.name for f in d.iterdir() if "restart" in f.name.lower() or f.suffix == ".rst7"]
    rec["restart_history"] = ";".join(sorted(restr))[:200] or "none"
    # trajectory frames (try candidates; some are gzip-compressed copies)
    if dcdfs:
        opened = False
        for cand in dcdfs:
            try:
                u = mda.Universe(str(prmtop), str(cand))
                dt = round(u.trajectory.dt, 3)  # ps
                n = len(u.trajectory)
                rec["frame_interval_ps"] = dt
                rec["frames"] = n
                rec["production_ns"] = round(dt * n / 1000.0, 2)
                rec["trajectory"] = cand.name
                opened = True
                del u
                break
            except Exception as e:  # noqa: BLE001
                rec["frames"] = f"ERR:{type(e).__name__}"
        if not opened and rec["frames"].startswith("ERR"):
            rec["frames"] += f" (tried {len(dcdfs)} files)"
    # boost log fields (skip comment lines)
    if logs:
        with open(logs[0]) as f:
            lines = [l.rstrip() for l in f if not l.lstrip().startswith(("#", "@"))]
        head = lines[0] if lines else ""
        ncol = len(head.split())
        rec["log_fields"] = f"{ncol} cols: {head[:100]}"
    rows.append(rec)

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"wrote {OUT} with {len(rows)} systems")
for r in rows:
    print(f"{r['system_id']:22s} {r['construct_id']} topo={r['topology_hash']} "
          f"frames={r['frames']:>6} dt={r['frame_interval_ps']}ps prod={r['production_ns']}ns "
          f"log={r['log_fields'][:40] if r['log_fields'] else 'NONE'}")
