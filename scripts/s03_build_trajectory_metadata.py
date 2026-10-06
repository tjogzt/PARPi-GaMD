#!/usr/bin/env python3
"""Build trajectory metadata table (reviewer Section 8 minimum verification package).

For every system directory under DATA_ROOT (top level) and the protonated
re-run roots (rerun_202609/runs_{nira,ruca}_prot), collect the topology hash,
trajectory frames, frame interval, boost log fields, restart history and
timestamps. Emits results/analysis/trajectory_metadata.csv (15 archival
systems; directory names of the protonated re-runs are canonicalised to
sys{1,2}_{nira,ruca}_prot).

Purpose:  Build the trajectory metadata table (frame counts, lengths, protocols) for all production runs.
Created:  2026-09-17 (header standardised 2026-10-05; coverage fix 2026-10-06)
Depends:  MDAnalysis, common.paths, hashlib
Run:      python3 scripts/s03_build_trajectory_metadata.py   (from the repository root)
"""
import hashlib
import os
import sys
import csv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths import data_root

DATA_ROOT = data_root()
OUT = Path(__file__).resolve().parents[1] / "results" / "analysis" / "trajectory_metadata.csv"

import MDAnalysis as mda  # noqa: E402


def md5(p, chunk=1 << 20):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()[:12]


def scan_dir(d, rows, system_id=None, seed="unknown"):
    """Collect one metadata record for a directory containing *.prmtop."""
    prmtops = list(d.glob("*.prmtop"))
    if not prmtops:
        return
    prmtop = prmtops[0]
    dcdfs = list(d.glob("*.dcd")) + list(d.glob("*.nc"))
    logs = list(d.glob("gamd.log"))
    if system_id is None:
        system_id = d.name
    compound = system_id.split("_", 1)[1] if "_" in system_id else "?"
    compound = compound.removesuffix("_prot")
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
        "random_seed": seed,
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
        if not opened and str(rec["frames"]).startswith("ERR"):
            rec["frames"] = f"{rec['frames']} (tried {len(dcdfs)} files)"
    # production segment: when the run records its production start step,
    # report production-only frames (start frame inferred at 2 fs MD steps,
    # inclusive slicing; prep ended at production-start-step). Display forms
    # for these re-run systems follow the archived convention: production_ns
    # spans the frame intervals ((n-1) x dt), and the frame interval is shown
    # as an integer picosecond value.
    is_prot = False
    pss = d / "production-start-step.txt"
    if pss.exists() and rec["frames"] and str(rec["frames"]).isdigit():
        is_prot = True
        dt = rec["frame_interval_ps"]
        start_steps = int(pss.read_text().strip())
        start_frame = round(start_steps * 0.002 / dt)
        n_prod = int(rec["frames"]) - start_frame + 1
        rec["frame_interval_ps"] = f"{dt:g}"
        rec["frames"] = n_prod
        rec["production_ns"] = round(dt * (n_prod - 1) / 1000.0, 2)
    # boost log fields (skip comment lines)
    if logs:
        with open(logs[0]) as f:
            lines = [l.rstrip() for l in f if not l.lstrip().startswith(("#", "@"))]
        head = lines[0] if lines else ""
        ncol = len(head.split())
        rec["log_fields"] = (f"{ncol} cols log col8 dV_D" if is_prot
                             else f"{ncol} cols: {head[:100]}")
    rows.append(rec)


rows = []
for d in sorted(DATA_ROOT.iterdir()):
    if d.is_dir():
        scan_dir(d, rows)

# Protonated re-run segments: canonical system ids differ from the on-disk
# directory layout (runs_{nira,ruca}_prot/{sys1,sys2}).
PROT_ROOTS = [
    (DATA_ROOT / "rerun_202609" / "runs_nira_prot" / "sys1", "sys1_nira_prot"),
    (DATA_ROOT / "rerun_202609" / "runs_nira_prot" / "sys2", "sys2_nira_prot"),
    (DATA_ROOT / "rerun_202609" / "runs_ruca_prot" / "sys1", "sys1_ruca_prot"),
    (DATA_ROOT / "rerun_202609" / "runs_ruca_prot" / "sys2", "sys2_ruca_prot"),
]
for d, sid in PROT_ROOTS:
    if d.is_dir():
        scan_dir(d, rows, system_id=sid, seed="42")

rows.sort(key=lambda r: r["system_id"])

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
