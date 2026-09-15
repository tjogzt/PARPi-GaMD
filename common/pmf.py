"""PyReweighting PMF helpers.

Single implementation of the well-depth / reweighting code that was previously
duplicated across eight analysis scripts (s2_new_drugs_*.py, s1_convergence_blocks.py,
gen_new_pmf_files.py, gen_extension_convergence_csv.py, s1_length_sensitivity.py,
s2_weight_protocol_test.py).
"""
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from common.paths import pyrew


class ReweightResult:
    """Outcome of one PyReweighting run (C1-C3 well depths)."""

    def __init__(self, tmpdir, c1, c2, c3, stdout, stderr):
        self.tmpdir = Path(tmpdir)
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.stdout = stdout
        self.stderr = stderr


def wd_from_pmf(f):
    """Well depth (pmf.max() - pmf.min()) of a 2-column xvg file."""
    d = []
    for line in open(f):
        line = line.strip()
        if line and not line.startswith(("#", "@")):
            p = line.split()
            d.append([float(p[0]), float(p[1])])
    pmf = np.array(d)[:, 1]
    return float(pmf.max() - pmf.min())


def run_pyrew(cv, weights, tmp=None, wfmt="%.4f", cvfmt="%.4f"):
    """Run PyReweighting-1D.py (amdweight_CE, T = 300) and return C1-C3 well depths.

    `weights` may be 1-D or an (n, 3) array (w1, w2, dV columns). `tmp` is the
    scratch directory (created under the system temp dir when omitted).
    """
    tmp = Path(tmp) if tmp else Path(tempfile.mkdtemp())
    np.savetxt(tmp / "cv.dat", cv, fmt=cvfmt)
    np.savetxt(tmp / "weights.dat", weights, fmt=wfmt)
    r = subprocess.run(
        ["python3", str(pyrew()), "-input", "cv.dat", "-T", "300", "-disc", "0.1",
         "-Emax", "20", "-cutoff", "2", "-job", "amdweight_CE",
         "-weight", "weights.dat"],
        cwd=tmp, capture_output=True, text=True)
    if not (tmp / "pmf-c3-cv.dat.xvg").exists():
        raise RuntimeError(
            "PyReweighting failed: " + r.stdout[-300:] + " " + r.stderr[-300:])
    c1 = wd_from_pmf(tmp / "pmf-c1-cv.dat.xvg")
    c2 = wd_from_pmf(tmp / "pmf-c2-cv.dat.xvg")
    c3 = wd_from_pmf(tmp / "pmf-c3-cv.dat.xvg")
    return ReweightResult(tmp, c1, c2, c3, r.stdout, r.stderr)
