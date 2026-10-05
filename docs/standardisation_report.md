# Repository standardisation report (2026-10-05)

Scope: full standardisation of the distributable repository (analysis code,
data layer, figures, documentation) for the PARPi-GaMD submission. Performed on
the working tree; `docs/rename_log.tsv` carries the complete old→new mapping.

## 1. Naming and code organisation

- All code files renamed to `<number>_<snake_case_description>.<ext>`:
  - `code/` analysis chain → `NN_*` (22 renamed from `NN-*`, 10 promoted from
    unnumbered, 3 added generators, 3 merged pairs → single entry points).
  - `code/pipeline/` deployment chain → `pNN_*` (8 files; `environment_gamd_min.yml`
    kept as-is).
  - `scripts/` build chain → `sNN_*` (10 files; `rebuild_all.sh` kept as the
    documented entry point).
- Merges: `gen_extension_convergence_csv.py` + `gen_new_pmf_files.py` →
  `code/41_extension_panel_data_files.py`; `s1_convergence_blocks.py` +
  `s1_length_sensitivity.py` → `code/42_extension_s1_convergence_checks.py`;
  `s2_new_drugs_pmf.py` + `s2_new_drugs_pmf_dfw.py` →
  `code/43_extension_s2_well_depths.py`.
- Excluded from distribution (kept in the working tree): `s2_new_drugs_aai.py`
  (superseded values), `s2_weight_protocol_test.py` (one-off), and the former
  `code/archive/` tree (superseded scripts; no longer shipped).
- Every distributed code file carries a standard English header: Purpose,
  Author, Created, Inputs, Outputs, Depends, Run. Console output/messages are
  English throughout (deployment scripts translated from Chinese).

## 2. Hard-coded business values moved into the data layer

| Script | Former in-code constants | Now read from |
|---|---|---|
| `code/34_helix_rmsf_heatmap.R` | Table S9 RMSF matrix (7 × 5) | `data/01_curated/helix_rmsf_table_s9.csv` |
| `code/31_extension_convergence_figure.R` | S1 wells 54.0/48.7/52.2 ± 8.0/5.8/7.8 | `data/01_curated/extension_s1_wells.csv` (C1–C3 recomputed from the archived cv/weights with PyReweighting; pop. SD) |
| `code/21_spearman_correlation.R` | Two-state panel vectors (S1/S2 CV1/CV2, all three states) | `data/01_curated/two_state_panel_values.csv`, asserted on load against `data/analysis/s2_dbe_final.csv` + `s1_dbe_unified_wells.csv` |
| `code/39_s2_cv1_retention_figure.R` | Table S5 span bars; palette | `data/01_curated/s2cv1_spans_table_s5.csv`; `common/palette_si.csv` |
| `code/18_replicate_analysis.R` | Table-4 replicate anchor values | `data/01_curated/manuscript_anchor_values.csv` (assertions) |
| `code/13_mechanism_figure.R` | AAI-SD anchor values | `data/01_curated/manuscript_anchor_values.csv` (assertions) |
| `code/31_*` | Extension colours | `common/ligands.csv` |
| `code/16_2d_landscape.R` | Dead legacy estimator (298 K constant) | removed (documented in-place) |

Seeds (e.g. `set.seed(49)`) are retained by design.

## 3. Path model

- Inputs: `data/**` (curated `01_curated/`, analysis-input `analysis/`, replicate,
  per-frame and weight archives) + `common/` + `$DATA_ROOT` for raw trajectories.
- Outputs: figures → `figures/pdf|png/`; curated tables → in place in `data/`;
  regenerable intermediates → `results/` (git-ignored, auto-created).
- All `/tmp` hard-codes removed (runner scripts use `results/`).
- `common/paths.py::analysis_dir()` now points at `data/analysis`;
  `results_dir()` added.

## 4. Data layer additions

- `data/analysis/`: the files the standardised scripts consume but that lived
  only in the working tree — per-system CV series and PMF-c3 curves (`.xvg`),
  2-D grids (`.npz`), derived tables (`s2_dbe_final.csv`,
  `s1_dbe_unified_wells.csv`, `s1_table1_rows.csv`, `d4_consistency_matrix.csv`,
  `convergence_diag.csv`, `timescale_scaling.csv`, `s2_window_timescale.csv`,
  `p0_2_trunc_scan.csv`, `S1_S2_ratios.csv`, `pmf_features_summary.csv`,
  `s2_dbe_*.npy`), `rmsf_recomp/` profiles, `ref_ca.pdb`, plus `01_modeling/`
  (receptor PDB + docked SDF/mol2) for the deployment chain, and
  `extension_s1_convergence.csv` at `data/`.
- New curated tables (provenance): `two_state_panel_values.csv` (SI Tables
  S1/S5 C3 convention; asserted against the primary CSVs),
  `s2cv1_spans_table_s5.csv` (Table S5 spans), `extension_s1_wells.csv`
  (C1–C3 recomputed from `data/new_drugs_s1/` cv/weights — C3 reproduces the
  manuscript 54.0/48.7/52.2 exactly), `helix_rmsf_table_s9.csv` (SI Table S9),
  `manuscript_anchor_values.csv` (Table 1/4 anchors for in-script asserts).

## 5. Figures

- `figures/` now contains exactly the 27 in-manuscript figures
  (10 main + 17 SI), as `pdf/` (+ `png/` 150-dpi previews, including the
  Figure-1 structure crop asset used by `code/15_fig1_*`).
- Added with generators: `Fig_D4_ConsistencyMatrix`, `Fig_D4_Timescale`,
  `Fig_D4_WindowNoise` (generators `code/29`, `code/33`), and
  `Fig_SI_Sampling_Toolkit` (generator `code/40`).
- Removed from the distribution (~29 superseded/auxiliary PDFs plus duplicated
  PNGs; full list in the git diff): e.g. `01-07-*`, `Fig5_Pocket_Dynamics`,
  `Fig_Allosteric_Pathway`, `Fig_DCCM_Network`, `Fig_Graphical_Abstract`,
  `Fig_Ligand_Properties`, `Fig_QSAR_*`, `Fig_RMSF_*`, `Fig_S2_Retention_Paradox`,
  `Fig_Transition_Barriers`, `TOC_Graphic*`, `Fig_SI_PMF_Convergence` (auxiliary;
  its script `code/30` now writes to `results/analysis/`).
  All remain in the working tree under `PARPi_design/figures/`.

## 6. Verification log (standardised tree, 2026-10-05)

| Check | Result |
|---|---|
| Python `py_compile` (all 40 py files) | pass |
| Shell `bash -n` (all sh files) | pass |
| R `parse` (all R files) | pass |
| Residual `results/` refs in scripts | 0 (only intended `results/` scratch outputs) |
| CJK characters in code | 0 |
| `code/21` re-run vs shipped figure | pixel diff 0 (100 dpi) |
| `code/34` re-run vs shipped figure | pixel diff 0 |
| `code/39` re-run vs shipped figure | pixel diff 0 |
| `code/31` re-run vs shipped figure | pixel diff 0 (after population-SD correction) |
| `code/32` re-run vs curated CSV | max cell deviation 0.0000 (asserts pass) |
| `code/30` re-run | writes to `results/analysis/` only |
| Metadata checksum script `scripts/s10_verify_manifest.py` | coherent with the updated `data_manifest.md` (regenerable paths flagged) |


- **Figure-style revision (2026-10-05, later round)**: open-frame panel style
  (top/right borders removed, left/bottom axis lines retained) applied to the
  single-panel SI figures S3, S4, S5, S10, S11 via `theme_open` in
  `common/helpers.R`, and to S12 via `code/39`. Grid/facet figures
  (S1/S2/S6/S7/S8/S9/S13/S14 and the S15-S17 base-R panels) intentionally keep
  their existing frames. The SI was recompiled and the submission package
  updated accordingly.
  Main-manuscript round: open-frame applied to Figure 4 (Fig_Trapping_vs_Allostery,
  both panels) and Figure 3 (Fig_2D_Mechanism_Map); the D4 diagnostics (Figures
  5-7) were already in open/minimal styles; grid/composite main figures
  (Figures 1, 2, 8-10) intentionally unchanged.

## 7. Known limits

- Raw trajectories are not distributed; scripts that read them require
  `$DATA_ROOT` (documented per script header).
- `code/pipeline/*` and `scripts/s04–s09` run on the GPU host (conda env
  `gamd`); they are shipped for reproducibility of the simulation protocol, not
  for laptop execution.
- `README.md` and `data_manifest.md` were rewritten/updated to the new layout.

## 8. Post-standardisation re-run verification (2nd pass, 2026-10-05)

All data-only scripts were executed in the standardised tree. Findings and fixes:

- **09 (S1 cumulants)**: two panel inputs (`sys1_*_pmf_c1/c2.xvg`) were missing
  from the tree; the script silently plotted one order fewer. Inputs added →
  figure now reproduces pixel-identical.
- **16 (2D landscapes)**: the per-system CV archives (`sys2_*_cv.npy`, 14 files)
  were missing; added → both landscape figures reproduce pixel-identical.
- **28 (PCA)**: the earlier re-run of `41`/`44` had refreshed the pmf-c3 xvg
  snapshots with drifted values, shifting the PCA features. Snapshots restored
  to the manuscript-consistent versions → pixel-identical; the PCA CSV outputs
  now write to `data/pca_*.csv` (canonical location).
- **35**: output path corrected to `data/cumulant_wells_C3.csv`. The
  regenerated table now matches the current (S12-corrected) xvg inputs: the
  previously shipped rows for niraparib/olaparib/rucaparib (42.1/45.3/47.0)
  were pre-rebuild values; the corrected rows agree with Table S5 within
  rounding (ruca CV1: 59.465 from the restored archive vs Table S5 59.3).
- **29 (D4 consistency) and 38 (DCCM)**: current script re-runs (in both the
  working tree and the standardised tree) do not pixel-reproduce the shipped
  versions — the shipped files correspond to an earlier data revision. They
  were kept as shipped to preserve manuscript↔repository figure identity.
- **41/44 caveat**: their xvg regeneration follows the archived weight files;
  the snapshots shipped in `data/analysis/` are the manuscript-consistent
  versions. A note to this effect was added to both scripts' headers.
- Total: 7 figure chains pixel-verified at 0 diff (09/16/21/28/31/34/39), plus
  30/32/35 checked against their outputs; `scripts/s10_verify_manifest.py`
  passes 11/11 checksums with 0 mismatches.



## 9. Figure-review round 1 (2026-10-05, user-approved; four-team review)

Applied fixes from the full figure-set review (全篇图件审查; full report: workspace `docs/figures_review_round1.md`):
- **Base-R font scale** (S15 `code/31`, S16 `code/34`): device/par pointsize 8 -> 12 with explicit
  ``par(cex = 1)`` after ``layout()`` (base R resets/reduces the base cex on layout); all sub-8pt
  cex values raised to >= 8pt. S15 extras: panel-B note split to two lines (panel-letter collision),
  bottom margin raised 3.6 -> 4.2 lines and panel-C axis title raised (bottom-edge safety),
  senaparib drawn dashed (line-type redundancy for the red pair, matching the legend), title/label
  retune verified by a full 8-pt word scan.
- **S17** (`code/40`): ``par(cex = 1)`` after ``par(mfrow = c(1, 3))`` (base R reduces the base cex
  to 0.66 for >=3-panel layouts), legend cex 0.72 -> 0.75, significance-note moved and split
  (two lines, lower-left; no collision with the p=0.13 label or the page edge), data-label offsets
  retuned (dot/label clearance >= 2pt), panel-A axis title split to two lines (slot clipping at 8.8pt).
- **ggplot fixes**: S12 (`code/39`) and Fig8 (`code/37`) legend text/title set to 8pt; Fig2
  (`code/13`) 45-degree axis text 7 -> 8pt; Fig10 (`code/18`) subtitle 7 -> 8pt.
- **Axis-line unification**: ``theme_open`` (helpers.R) and S12's inline theme now ``black/0.5``
  (was grey20/0.4), matching Fig1/Fig2 panels and the documented axis standard.
- **Label harmonisation**: figure texts now use the en-dash ``Protein–DNA`` form matching the
  manuscript (scripts 10/16/21/39); Fig3 legend ``Unknown`` -> ``Unclassified`` (display-only via
  scale ``labels``) plus a marker-classes sentence added to the Fig3 caption; S17 caption and
  SI section 2.10 seed descriptions aligned with the data (talazoparib 3 x 200 ns; veliparib
  1 x 200 ns + 2 x 400 ns, REV-1/REV-2).
- **Housekeeping**: unreferenced ``Fig_S2_pmf_overlay_combined.pdf`` removed from this tree
  (see data_manifest.md note); all 27 shipped figures re-verified.

Verification: 8-pt scan passes on all 27 shipped figures (one documented false positive: a
rotated label fragment in Fig8 measured 7.8pt while its text span is 8.00pt); geometry audit
0 overlaps / 0 edge words on all figures; pixel-diffs of every touched figure restricted to the
intended regions; main manuscript (57 pp) and SI (25 pp) recompiled with the new figures and
feature-checked (caption sentences, in-figure notes, en-dash labels present in the compiled PDFs).

### 9b. Figure-review round 1 — follow-up batch (2026-10-05, continued)

Executed after the main round, same approval:

- **Numeric-precision unification (p-values -> 3 decimals).** All p-values in the manuscript
  and SI now use 3 decimal places, matching the sensitivity table and the computed sources
  (``data/analysis/d4_consistency_matrix.csv``, ``data/analysis/p0_2_trunc_scan.csv``):
  e.g. 0.77 -> 0.767, 0.63 -> 0.633, 0.50 -> 0.500, 0.27 -> 0.267, 0.13 -> 0.133, 0.65 -> 0.650,
  0.85 -> 0.850, 0.20 -> 0.200, 1.0 -> 1.000 (incl. the two ``(1.0)`` cells in the sensitivity
  table); the bound ``p >= 0.43`` was tightened to ``p >= 0.433``. Conventional inequality
  thresholds (``p < 0.05``, ``p >= 0.70``) and the structural probabilities
  (``1/120 ~= 0.008``, ``0.083``, ``0.356``) keep their existing forms.
- **Hyphen unification extended: ``HD-ART`` -> ``HD–ART`` in figure text.** All 10 shipped figures that displayed the
  hyphenated form (scripts 09/10/11/13/16) now use the en dash, matching the manuscript
  (main 39x, SI 23x ``HD--ART``); all remaining hyphen instances in scripts (comments, console
  strings) were also converted, and the one stray SI comment was fixed. Zero hyphen residues
  remain in the shipped figure PDFs.
- **S17 p-labels.** ``code/40`` now prints p-labels with ``%.3f`` (0.133 / 0.433 / 0.767 / 0.833 /
  1.000), consistent with the manuscript's p-value family. The rightmost label exceeded the
  panel region and had its last glyph clipped; it is now drawn with a small inward offset
  (dx = -1.8 data units for T > 25 ns) so the full string stays inside the panel.
- **Tree cleanup.** Re-running script 11 regenerates three non-shipped by-products
  (``Fig_Master_S1S2_HD_ART``, ``Fig_S1S2_hd_art_overlay``, ``Fig_S1S2_hd_art_panels``); removed
  from this tree after the run (same policy as ``Fig_S2_pmf_overlay_combined``). Stale header
  comment in ``code/40`` corrected (Figure S16 -> S17).
- **Re-verification.** 8-pt scan 27/28 (single documented false positive: rotated fragment in
  Fig8, span = 8.00 pt); geometry 0 overlaps / 0 edge words on all 27; pixel-diffs restricted to
  the intended regions; main (57 pp) and SI (25 pp) recompiled; packages refreshed; repository
  synced.

- **Fig10 legend label (found during post-sync verification).** The replicate-validation figure
  (``code/18``) displayed the legacy legend term ``Unknown``; it now shows ``Unclassified``
  (display-only via scale ``labels``, same treatment as Fig3). Re-rendered, redistributed;
  the main manuscript was recompiled (57 pp) and re-verified (main ``Unknown`` = 0).

- **Fig1 refinements (user-requested).** Panel C: the S1 box widened (x 0.545-1.455 -> 0.40-1.60)
  so the atom-count line sits fully inside the frame (measured margins 3.6 / 3.3 pt); Panel A: the
  ``CAT Pocket`` annotation and its leader line shifted right (780 / 870 -> 805 / 895 data units)
  to clear the adjacent ``Type I (EB-47)`` label (inter-label gap -2.2 -> +3.0 pt). Re-rendered,
  redistributed, manuscript recompiled (57 pp).
  Follow-up (same day, user-requested): the 3D structure in Panel B was sunk by scaling it
  0.92 about its bottom edge (centred), lifting the top of the protein clear of the note line
  (ink-top to note-bottom gap -2.9 -> +6.9 pt); labels were remapped with the same transform.
- **Fig2 Panel A bottom assembly (user-requested).** The 10-entry two-row legend was wider than
  the panel and its first key was clipped at the page edge; an asymmetric ``legend.margin``
  (left +18 pt) shifts the wrapped legend right (first key -5.5 -> +3.5 pt, right edge clears the
  neighbouring legend by 25 pt). The axis title and legend were pulled toward the axis
  (tick-to-title 26.3 -> 11.0 pt; title-to-legend 13.0 -> 11.0 pt) via explicit
  ``axis.title.x`` / ``axis.text.x`` / ``legend.margin`` margins. Re-rendered; 8-pt and geometry
  checks clean; manuscript recompiled (57 pp).
- **S17 panel spacing (user-requested).** The three-panel layout had oversized inter-panel margins
  (edge-to-edge gaps 79/76 pt vs ~30 pt in the ggplot figures). Margins unified to
  ``mar = c(4.2, 3.0, 2.5, 0.5)`` for all three panels: each plot region +~20 pt wider, gaps
  reduced to 58/54 pt (visual whitespace 22/18 pt), consistent with the journal style. Labels
  re-verified collision-free; all previous S17 fixes (3-decimal p labels, rightmost-label shift,
  two-line title) retained; SI recompiled (25 pp).
