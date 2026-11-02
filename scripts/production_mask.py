#!/usr/bin/env python3
"""production_mask.py — single source of truth for the GaMD production mask.

Only archives with an explicit production-start record are masked here.  The
veliparib protonated S2 archive (``veli_prot_s2_frames.csv.gz``) stores both
equilibration and production frames; production starts at MD step 6,000,000
(``production-start-step.txt`` in the run directory, ``extract_cv.log``:
"production_frames=521").  Analyses that slice ``cv[:280]`` without this mask
silently use 239 equilibration frames and only 41 production frames.

Archives that were already trimmed to production (frame indices 1..N, or
time-in-ns columns) are returned unchanged; do not apply the 6e6 threshold to
them, because their step columns mean different things.

Metadata caveat: the legacy column name ``step_ps`` is ambiguous.  In archives
written by s31_archive_frame_cvs.py it holds the 1-based frame index, whereas
in veli_prot_s2 it holds the MD step.  Never infer units from the column
name; use the per-archive PRODUCTION_STARTS entry or the time_ns column.
"""

PRODUCTION_STARTS = {
    "veli_prot_s2": 6_000_000,
}
MD_STEP_COLUMNS = ("step_ps", "step", "step_md")


def step_column(header):
    for name in MD_STEP_COLUMNS:
        if name in header:
            return name
    return None


def mask_production(name, header, arr):
    """Drop pre-production rows for archives in PRODUCTION_STARTS."""
    start = PRODUCTION_STARTS.get(name)
    if start is None:
        return arr, 0
    col = step_column(header)
    if col is None:
        raise KeyError(f"{name}: no MD-step column in header {header}")
    steps = arr[:, header.index(col)]
    keep = steps >= start
    dropped = int((~keep).sum())
    if dropped:
        print(f"[production_mask] {name}: dropped {dropped} pre-production rows "
              f"(< {start}); kept {int(keep.sum())}", flush=True)
    return arr[keep], dropped


def assert_veli_prot(n, expected=521):
    if n != expected:
        raise AssertionError(
            f"veli_prot production frames = {n}, expected {expected}; "
            "the 6,000,000-step mask is missing or the archive changed")
