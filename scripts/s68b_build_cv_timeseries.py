#!/usr/bin/env python3
"""s68b_build_cv_timeseries.py — tidy CV1/CV2 series for the trajectory-QC figure.

Reads the deposited per-frame archives (50-ps cadence), maps frame index to
time (0.05 ns/frame) and truncates at 26 ns. Writes
data/analysis/s68_cv_timeseries.csv.
"""
import csv
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/analysis/s68_cv_timeseries.csv"

FILES = {
    "talazoparib": "orig_s2_talazoparib_frames.csv.gz",
    "niraparib": "sys2_niraparib_frames.csv.gz",
    "olaparib": "sys2_olaparib_frames.csv.gz",
    "rucaparib": "sys2_rucaparib_frames.csv.gz",
    "veliparib": "orig_s2_veliparib_frames.csv.gz",
    "niraparib_prot": "sys2_nira_prot_frames.csv.gz",
    "rucaparib_prot": "sys2_ruca_prot_frames.csv.gz",
    "veliparib_prot": "veli_prot_s2_frames.csv.gz",
}


def main():
    rows = []
    for system, fname in FILES.items():
        path = ROOT / "data/per_frame" / fname
        with gzip.open(path, "rt") as h:
            next(h)
            for line in h:
                p = line.replace(",", " ").split()
                if len(p) < 5:
                    continue
                frame = int(float(p[0]))
                t = round(frame * 0.05, 3)
                if t > 26.0:
                    continue
                rows.append((system, t, round(float(p[3]), 4), round(float(p[4]), 4)))
        print(f"{system}: {sum(1 for r in rows if r[0] == system)} rows")
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["system", "time_ns", "cv1", "cv2"])
        w.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
