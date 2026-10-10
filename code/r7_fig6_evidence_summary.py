#!/usr/bin/env python3
"""r7_fig6_evidence_summary.py — Figure 6: diagnostic workflow + evidence-grade summary.

Schematic (no data): workflow strip + six-diagnostic verdict table + final verdict band.
All statements mirror the frozen manuscript language; numbers as in R6 master.
Output: results/figures/Fig6_Evidence_Summary.pdf + .png (172 mm wide).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

W_MM, H_MM = 172, 84
fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

C_BLUE, C_RED, C_SOFT = "#3D6BA8", "#C23531", "#F7F3EE"

def box(x, y, w, h, text, fc="white", ec=C_BLUE, fs=8.0, weight="normal", tc="black"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.35,rounding_size=1.2",
                                fc=fc, ec=ec, linewidth=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            fontweight=weight, color=tc, linespacing=1.25)

def arrow(x1, x2, y):
    ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>",
                                 mutation_scale=8, color="#555555", linewidth=0.9))

# ---- workflow strip ----
y0, h0 = 85.5, 11.5
box(0.5, y0, 17.0, h0, "Trajectories\n(12 systems, 4\nchemical-state panels)", fs=7.8)
box(20.5, y0, 19.5, h0, "Descriptors\nS1 span; S2 CV1/CV2\nspans; AAI", fs=7.8)
box(43.0, y0, 15.5, h0, "Association\nvs trapping rank\n(n = 5)", fs=7.8)
box(61.5, y0, 13.5, h0, "Six\ndiagnostics\n(D1\u2013D6)", fs=7.8)
box(78.0, y0, 21.5, h0, "Evidence grade\n\u201cnot established\u201d\nunder any state panel", fs=7.8, weight="bold")
arrow(17.8, 20.2, y0 + h0 / 2); arrow(40.3, 42.7, y0 + h0 / 2)
arrow(58.8, 61.2, y0 + h0 / 2); arrow(75.3, 77.7, y0 + h0 / 2)

# ---- verdict table ----
rows = [
    ("D1  Windowed estimator variance",
     "noise \u2265 signal (0.3\u20131.2\u00d7): ordering unresolvable at 26 ns"),
    ("D2  Estimator\u2013consistency matrix",
     "no cell significant (all p \u2265 0.27); sign unstable across estimators"),
    ("D3  Convergence\u2013timescale scaling",
     "no decay over 2\u201324 ns windows; no extrapolated sampling length"),
    ("D4  Ordinal sensitivity model",
     "no decision-level cell (P \u2264 0.78 across assumed bands)"),
    ("D5  Counterfactual sensitivity",
     "leave-one-out \u22120.63\u2026+0.63; single-point exclusion flips sign"),
    ("D6  Chemical-state audit",
     "protonation corrections move \u03c1 \u22120.67 \u2192 \u22120.05: state dominates"),
]
x0, wd = 0.5, 99.0
ytop, rh = 78.5, 9.2
ax.text(x0, 81.5, "Six diagnostics \u2192 outcomes in this study (n = 5, tie-aware exact permutation)",
        fontsize=8.4, fontweight="bold", va="bottom")
for i, (name, verdict) in enumerate(rows):
    y = ytop - (i + 1) * rh
    fc = C_SOFT if i % 2 == 0 else "white"
    ax.add_patch(FancyBboxPatch((x0, y + 0.6), wd, rh - 1.2, boxstyle="round,pad=0.2,rounding_size=0.8",
                                fc=fc, ec="#DDDDDD", linewidth=0.5))
    ax.text(x0 + 1.6, y + rh / 2, name, fontsize=8.2, fontweight="bold", va="center", color=C_BLUE)
    ax.text(x0 + 33.5, y + rh / 2, verdict, fontsize=8.2, va="center")

# ---- final verdict band ----
box(0.5, 1.5, 99.0, 9.0,
    "Verdict:  no robust association between the computed spans and trapping potency \u2014 not established under any\nchemically defensible ligand-state panel.  The transferable output is the audit workflow itself.",
    fc="#FBEDEB", ec=C_RED, fs=8.6, weight="bold", tc="#7A0F0B")

fig.savefig("results/figures/Fig6_Evidence_Summary.pdf")
fig.savefig("results/figures/Fig6_Evidence_Summary.png", dpi=220, facecolor="white")
print("WROTE Fig6_Evidence_Summary.pdf/png")
