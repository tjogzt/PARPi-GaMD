# -*- coding: utf-8 -*-
#
# Purpose:  JCIM TOC graphic, single proportional design scale S, current frozen wording.
# Created:  2026-10 (header standardised 2026-10-11)
# Run:      python3 code/56_toc_graphic.py   (from the repository root)

"""56_toc_graphic.py — JCIM TOC graphic, single proportional design scale S, frozen wording.

Rebuilt from the scratch generator toc_v4/v6 with corrections:
  * "Bayesian ordinal" -> "Ordinal sensitivity" (R6d terminology: the model is an
    ordinal sensitivity model; Bayesian framing withdrawn)
  * stale correlation "-0.82" -> current neutral-panel "-0.67" (R6 freeze)
  * left label kept as "six-diagnostic audit workflow" (R6 positioning:
    not pre-specified / not a registered plan)
  * 2026-10-13: type scale set to 8-11 pt at the base design grid (base type was
    3.4-3.9 pt and illegible; wording and numerals untouched, only line breaks).
  * 2026-10-14 (user directive): the design is drawn at ONE proportional scale S —
    canvas, geometry, strokes and type all scale together — so the graphic never
    looks sparse (canvas enlarged with type held) nor oversized (type enlarged with
    the set). S = 1.5 -> 123.9 x 66.7 mm canvas, type 12 / 14.25 / 16.5 pt; shown
    1:1 (12.4 cm) on the last page. Adjust S alone to rescale everything together.
Outputs: figures/TOC_Graphic_jcim.pdf + .png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

S = 1.0                      # official ACS TOC size (3.25 x 1.75 in); single proportional scale
BW, BH = 975, 525            # base design grid (px at 300 dpi = 3.25 x 1.75 in)
W0, H0 = round(BW * S), round(BH * S)
DPI = 300.0
fig = plt.figure(figsize=(W0 / DPI, H0 / DPI), dpi=DPI)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W0 / DPI); ax.set_ylim(0, H0 / DPI); ax.axis("off")
PX = lambda v: v * S / DPI           # base-grid x -> figure inches
PY = lambda v: (H0 - v * S) / DPI    # base-grid y (top-down) -> figure inches

blue = "#3D6BA8"; fill = "#F5F2EC"; gray = "#6F6F6F"; red = "#C23531"

# -- six boxes -- (base-grid geometry; type = 8 pt x S)
x0s = [8, 192]; rows = [(56, 194), (208, 346), (360, 498)]
labels = [["Windowed\nestimator\nvariance", "Estimator\nmatrix"],
          ["Timescale\nscaling", "Ordinal\nsensitivity"],
          ["Counter-\nfactuals", "Chemical-\nstate audit"]]
for r, (ya, yb) in enumerate(rows):
    for c in range(2):
        xa = x0s[c]
        ax.add_patch(FancyBboxPatch((PX(xa), PY(yb)), PX(176), PX(138),
                     boxstyle=f"round,pad=0,rounding_size={0.04*S}", fc=fill, ec=blue, lw=1.0*S))
        ax.text(PX(xa + 88), PY((ya + yb) / 2), labels[r][c], ha="center", va="center",
                fontsize=8.0*S, color=gray, family="Arial", linespacing=1.25)

# -- center column: process only, no sample size and no result numbers --
ax.text(PX(500), PY(90), "pre-publication\naudit", ha="center", va="center",
        fontsize=8.0*S, color=blue, family="Arial", linespacing=1.3)
ax.text(PX(500), PY(240), "state\nestimator\nsampling\nrank", ha="center", va="center",
        fontsize=8.0*S, color=red, family="Arial", linespacing=1.25)
ax.add_patch(FancyArrowPatch((PX(500), PY(300)), (PX(500), PY(352)), arrowstyle="-|>",
             mutation_scale=8*S, lw=2.0*S, color="#555555"))
ax.text(PX(500), PY(410), "audit before\ninference", ha="center", va="center",
        fontsize=8.0*S, color=blue, family="Arial", linespacing=1.3)

# -- right box: method deliverable, not a study verdict --
ax.add_patch(FancyBboxPatch((PX(630), PY(293)), PX(335), PX(203),
             boxstyle=f"round,pad=0,rounding_size={0.05*S}", fc=red, ec="none"))
ax.text(PX(797.5), PY(128), "auditable\nevidence chain", ha="center", va="center",
        fontsize=9.5*S, color="white", family="Arial", fontweight="bold", linespacing=1.2)
ax.text(PX(797.5), PY(222), "reusable workflow for\nenhanced-sampling\nassociation claims",
        ha="center", va="center", fontsize=8.0*S, color="white", family="Arial", linespacing=1.3)

# -- lower right: subject only --
ax.text(PX(797.5), PY(340), "PARP1 GaMD\ncase study", ha="center", va="center",
        fontsize=8.0*S, color=gray, family="Arial", linespacing=1.3)
ax.text(PX(797.5), PY(430), "open, reproducible\nworkflow", ha="center", va="center",
        fontsize=8.0*S, color=gray, family="Arial", linespacing=1.3)

# -- top-left label --
ax.text(PX(10), PY(29), "six-diagnostic audit workflow", ha="left", va="center",
        fontsize=8.0*S, color=blue, family="Arial")

fig.savefig("figures/TOC_Graphic_jcim.png", dpi=DPI)
fig.savefig("figures/TOC_Graphic_jcim.pdf")
print(f"WROTE figures/TOC_Graphic_jcim.pdf/.png  (S={S}, canvas {W0/DPI*25.4:.1f} x {H0/DPI*25.4:.1f} mm)")
