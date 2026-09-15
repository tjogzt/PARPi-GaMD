"""Shared path configuration.

DATA_ROOT: environment variable pointing at the external trajectory/data
volume (the directory containing the sys1_*/sys2_*/md_analysis trees).
Scripts that read raw trajectories require it; scripts that only consume
results/ under the repository root do not.
"""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def data_root():
    """Return the external trajectory data directory (DATA_ROOT env var)."""
    v = os.environ.get("DATA_ROOT", "").strip()
    if not v:
        raise RuntimeError(
            "DATA_ROOT is not set. Export the trajectory data directory, e.g.:\n"
            "  export DATA_ROOT=/path/to/PARPi_data")
    return Path(v)


def analysis_dir():
    """results/analysis under the repository root."""
    return REPO_ROOT / "results" / "analysis"


def pyrew():
    """Path of the bundled PyReweighting-1D.py tool."""
    return REPO_ROOT / "tools" / "PyReweighting" / "PyReweighting-1D.py"
