"""Kabsch optimal-rotation alignment.

Single implementation replacing the four independent copies that existed in
control_tala_rmsf.py, recompute_all_s2.py, 19-rmsf_all_systems.py and
recompute_rmsf_dccm.py (archived). The sign-corrected formulation below is
numerically identical to the earlier variants on proper-rotation data (the
reflection branch is unreachable for protein trajectories).
"""
import numpy as np


def kabsch_align(mobile, ref):
    """Align `mobile` onto `ref` by optimal rotation.

    Returns the centered aligned coordinates (mobile frame, mean-subtracted).
    """
    ref_c = ref - ref.mean(axis=0)
    mob_c = mobile - mobile.mean(axis=0)
    H = mob_c.T @ ref_c
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return mob_c @ (Vt.T @ np.diag([1.0, 1.0, d]) @ U.T)


def kabsch_rmsd(mobile, ref):
    """Return (RMSD, aligned coordinates translated back to the ref frame)."""
    ref_c = ref - ref.mean(axis=0)
    aligned = kabsch_align(mobile, ref)
    rmsd = float(np.sqrt(np.mean((aligned - ref_c) ** 2)))
    return rmsd, aligned + ref.mean(axis=0)
