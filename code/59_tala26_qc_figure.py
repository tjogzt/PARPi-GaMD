#!/usr/bin/env python3
"""59_tala26_qc_figure.py — SI Figure S17: structural QC of the independent
26-ns talazoparib seed-2026 trajectory.
Panels: (A) CA RMSD; (B) Zn-SG minimum distances; (C) HD/ART per-residue
RMSF; (D) HD x ART DCCM block.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT.parent / "results" / "figures" / "Fig_SI_Tala26_Structure_QC"
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.labelsize": 8,
    "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.7, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
})

TALA = "#FF7F00"
BLUE, RED, GREEN = "#3D6BA8", "#C23531", "#4DAF4A"
q = pd.read_csv(ROOT / "data/analysis/s73_tala26_qc_frames.csv")
r = pd.read_csv(ROOT / "data/analysis/s73_tala26_rmsf_residue.csv")
dccm = np.loadtxt(ROOT / "data/analysis/s73_tala26_dccm_hd_art.csv", delimiter=",")

fig = plt.figure(figsize=(7.2, 6.0))
gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[1.15, 1],
                       height_ratios=[1, 1], wspace=0.34, hspace=0.36)

# A: RMSD
ax = fig.add_subplot(gs[0, 0])
ax.plot(q.time_ns, q.rmsd_ca, color=TALA, lw=1.0)
ax.axvline(14, color="0.35", ls="--", lw=0.7)
ax.set_xlabel("Production time (ns)")
ax.set_ylabel("CA RMSD (Å)")
ax.set_xlim(0, 26)
ax.set_ylim(0, 4.0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.text(0.02, 0.96, "A", transform=ax.transAxes, va="top", ha="left",
        fontweight="bold", fontsize=9)
ax.text(14.3, 3.72, "14 ns", fontsize=6.5, color="0.35")

# B: Zn-SG
ax = fig.add_subplot(gs[0, 1])
for col, lab, c in (("zn1_min", "Zn1", BLUE), ("zn2_min", "Zn2", RED),
                    ("zn3_min", "Zn3", GREEN)):
    ax.plot(q.time_ns, q[col], lw=0.8, color=c, label=lab)
ax.axhline(2.8, color="0.45", ls=":", lw=0.7)
ax.set_xlabel("Production time (ns)")
ax.set_ylabel("Min Zn–S$\\gamma$ (Å)")
ax.set_xlim(0, 26)
ax.set_ylim(1.75, 2.25)
ax.legend(frameon=False, ncol=3, loc="upper left", handlelength=1.0,
          columnspacing=0.8, borderpad=0.2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.text(0.02, 0.96, "B", transform=ax.transAxes, va="top", ha="left",
        fontweight="bold", fontsize=9)

# C: RMSF
ax = fig.add_subplot(gs[1, 0])
hd = r.dropna(subset=["rmsf_hd"])
art = r.dropna(subset=["rmsf_art"])
ax.plot(hd.resid, hd.rmsf_hd, color=BLUE, lw=0.9, label=f"HD (mean {hd.rmsf_hd.mean():.2f} Å)")
ax.plot(art.resid, art.rmsf_art, color=RED, lw=0.9, label=f"ART (mean {art.rmsf_art.mean():.2f} Å)")
ax.axvspan(662, 787, color=BLUE, alpha=0.06, lw=0)
ax.axvspan(788, 1014, color=RED, alpha=0.06, lw=0)
ax.set_xlabel("Residue (UniProt)")
ax.set_ylabel("RMSF (Å)")
ax.set_xlim(650, 1020)
ax.set_ylim(0, 9.0)
ax.legend(frameon=False, loc="upper left", handlelength=1.2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.text(0.02, 0.96, "C", transform=ax.transAxes, va="top", ha="left",
        fontweight="bold", fontsize=9)

# D: DCCM
ax = fig.add_subplot(gs[1, 1])
im = ax.imshow(dccm, cmap="RdBu_r", vmin=-1, vmax=1, origin="lower",
               aspect="auto", interpolation="nearest")
ax.set_xlabel("ART residue")
ax.set_ylabel("HD residue")
xt = np.linspace(0, dccm.shape[1] - 1, 5).astype(int)
yt = np.linspace(0, dccm.shape[0] - 1, 5).astype(int)
ax.set_xticks(xt); ax.set_xticklabels([str(788 + int(v)) for v in xt])
ax.set_yticks(yt); ax.set_yticklabels([str(662 + int(v)) for v in yt])
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label("Correlation", fontsize=7)
cb.ax.tick_params(labelsize=6, width=0.5, length=2)
ax.text(0.02, 0.96, "D", transform=ax.transAxes, va="top", ha="left",
        fontweight="bold", fontsize=9, color="black")

fig.subplots_adjust(left=0.09, right=0.97, bottom=0.09, top=0.97)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(str(OUT) + ".pdf", dpi=300)
fig.savefig(str(OUT) + ".png", dpi=300)
print("wrote", OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
