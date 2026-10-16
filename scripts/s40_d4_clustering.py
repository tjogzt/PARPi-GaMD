#!/usr/bin/env python3
#
# Purpose:  D4-3: Unsupervised clustering of the 7 S2 systems in CV1xCV2 space,
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s40_d4_clustering.py   (from the repository root)

"""D4-3: Unsupervised clustering of the 7 S2 systems in CV1xCV2 space,
tested against the literature Type II/III labels (ARI).

For each system, the production-segment joint (CV1, CV2) samples are used as the
system's distribution; pairwise 2-Wasserstein distances -> spectral clustering
(k = 2..3) -> ARI vs Type II/III labels.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import sys
from pathlib import Path
from scipy.stats import gaussian_kde
from itertools import combinations

from common.paths import data_root

DATA = data_root()
A = Path("results/analysis")

systems = {
    "talazoparib": ("Type_II",  None,              None),
    "veliparib":    ("Type_III", None,              None),
    "AZD5305":      ("Unknown",  None,              None),
    "APO":          ("APO",      None,              None),
    "niraparib":    ("Type_III", "rerun_202609",    None),
    "olaparib":     ("Type_II",  "rerun_202609",    None),
    "rucaparib":    ("Type_III", "rerun_202609",    None),
}


def load_cv_pair(system):
    """Return (cv1, cv2) arrays over production frames."""
    if system in ("AZD5305", "veliparib"):
        d = A / "new_drugs_s2" / system
        cv1 = np.loadtxt(d / "analysis_CV1.dat")
        cv2 = np.loadtxt(d / "analysis_CV2.dat")
    else:
        cv1 = np.load(A / f"sys2_{system}_CV1_cv.npy")
        cv2 = np.load(A / f"sys2_{system}_CV2_cv.npy")
    n = min(len(cv1), len(cv2))
    return cv1[:n], cv2[:n]


def wasserstein2(a, b):
    """2-Wasserstein distance between two 2D point sets (mean/cov difference)."""
    ma, mb = a.mean(0), b.mean(0)
    Ca, Cb = np.cov(a.T), np.cov(b.T)
    dmean = np.linalg.norm(ma - mb)
    e, V = np.linalg.eigh(Ca)
    sq = V @ np.diag(np.sqrt(np.clip(e, 0, None))) @ V.T
    m = sq @ Cb @ sq
    e2, _ = np.linalg.eigh(m)
    dtrace = np.trace(Ca) + np.trace(Cb) - 2 * np.sum(np.sqrt(np.clip(e2, 0, None)))
    return dmean ** 2 + max(dtrace, 0.0)


def spectral_cluster(W, k):
    """Normalized spectral clustering of a similarity matrix."""
    from scipy.linalg import eigh
    d = W.sum(1)
    Dn = np.diag(1.0 / np.sqrt(np.clip(d, 1e-12, None)))
    L = Dn @ W @ Dn
    evals, evecs = eigh(L)
    X = evecs[:, -k:]
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    from scipy.cluster.vq import kmeans2
    _, labels = kmeans2(X, k, minit="++")
    return labels


def ari(labels_a, labels_b):
    from itertools import combinations as comb
    n = len(labels_a)
    pairs = list(comb(range(n), 2))
    a = sum(labels_a[i] == labels_a[j] for i, j in pairs)
    b = sum(labels_b[i] == labels_b[j] for i, j in pairs)
    both = sum((labels_a[i] == labels_a[j]) and (labels_b[i] == labels_b[j])
               for i, j in pairs)
    import math
    e = a * b / math.comb(n, 2)
    mx = (a + b) / 2
    return (both - e) / (mx - e) if mx != e else 0.0


names = list(systems)
pairs = {nm: load_cv_pair(nm) for nm in names}
# subsample to a common length for fairness
n_use = min(len(p[0]) for p in pairs.values())
rng = np.random.default_rng(7)
samples = {nm: rng.choice(len(pairs[nm][0]), n_use, replace=False) for nm in names}
pts = {nm: np.column_stack([pairs[nm][0][samples[nm]], pairs[nm][1][samples[nm]]])
       for nm in names}

W = np.zeros((len(names), len(names)))
for i, a in enumerate(names):
    for j, b in enumerate(names):
        W[i, j] = np.exp(-wasserstein2(pts[a], pts[b]) / 2.0)

true_lab = [{"Type_II": 0, "Type_III": 1}.get(systems[n][0], 2 if systems[n][0] == "Unknown" else 3)
            for n in names]
for k in (2, 3):
    lab = spectral_cluster(W, k)
    print(f"spectral k={k}: labels {dict(zip(names, lab.tolist()))}")
    print(f"  ARI vs TypeII/III = {ari(true_lab, lab):.3f}")

np.save(A / "d4_wasserstein_W.npy", W)
with open(A / "d4_clustering.txt", "w") as f:
    f.write("\n".join(f"{a} {b} {W[i,j]:.4f}" for i, a in enumerate(names)
                      for j, b in enumerate(names)) + "\n")
print("\nWROTE results/analysis/d4_wasserstein_W.npy, d4_clustering.txt")
