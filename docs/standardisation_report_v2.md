# standardisation_report_v2.md — second-round reconciliation (2026-10-11)

> Scope: integrate every R6–R10 addition/revision into the standardised repository
> (`submission/PARPi-GaMD/`), keep every shipped artifact reproducible, and record
> the adjudication evidence. Supersedes/extends `standardisation_report.md` (round 1,
> 2026-10-05). Working archive: `review-stage/release_v2/` (not distributed).

## 1. Core principle — frozen-figure consistency

The shipped figures (manuscript + SI) are FROZEN and user-approved. Where a file
existed in two edited branches (the working analysis tree vs the standardised
snapshot), the revision that *generated the shipped figure* is authoritative; the
other branch donates format conventions only (7-field header, path conventions,
value-sourcing from curated CSVs when value-identical, English text).

Evidence used to identify the render stream per file (all recorded in the archive):

- **Embedded PDF CreationDate ≈ script mtime (second-level)** for the R6/R9
  re-renders (e.g. `Fig1` 10-10 01:52:20 vs `15_fig1_domain_arch.R` 01:52:19).
- **Figure text-layer forensics** (pdftotext): e.g. frozen Fig 1 contains the
  combined "CAT pocket:" label (working-tree revision); frozen 2D map contains two
  "(protonated)" entries; frozen S2/S1 panels contain en-dash labels and
  "Unclassified" (standardised-snapshot revision) — each verdict matched the
  CreationDate evidence.

## 2. Reconciled files (round-2 class A)

Eleven files were edited in both branches after the round-1 snapshot. Final content
= working-tree revision (render-consistent), carrying round-1 format conventions:

| file | decision points |
|---|---|
| 12_descriptive_analysis.R | R6 well-depth audit retained; CSVs written to `data/analysis/` |
| 13_mechanism_figure.R | R6-support values from `r6_panel_values.csv`; both map variants kept as shipped |
| 15_fig1_domain_arch.R | shipped geometry (combined CAT label) kept; crop read from `figures/png/` |
| 28_pca_features.R | R6 well-depth substitution kept; reads `data/analysis/Fig_Mechanism_Data.csv` |
| 29_d4_consistency_matrix.R | R6 rebuild (4-panel) kept; byproduct PDF/PNG → `results/figures/` |
| 31_extension_convergence_figure.R | audited extension values re-sourced from `extension_s1_wells.csv` (see §4) |
| 33_d4_diagnostics_figures.R | R9 "2.6–5.3" separations now DERIVED from `r6_panel_values.csv` (guard asserts) |
| 34_helix_rmsf_heatmap.R | matrix re-sourced from updated `helix_rmsf_table_s9.csv` (SI Table S9 values) |
| 38_dccm_figure.R | reads `data/analysis/rmsf_recomp/`; shipped axis breaks kept |
| 39_s2_cv1_retention_figure.R | spans from `s2cv1_spans_table_s5.csv`, palette from `common/palette_si.csv` (value-identical) |
| 40_sampling_toolkit_figure.R | SI figure number S17 (per si.tex); paths to `data/analysis/` |

Ten further "candidate diverge" files were adjudicated to the **standardised-snapshot
revision** (= the render stream for their figures): 09, 10, 11, 16, 18, 01, 30, 32,
05, 08. No change required (already in the tree).

Repo-stream style refinements that were never rendered (font pins, axis-limit
tweaks, legend moves, open-frame additions in 09/10/13/28/31/39/40) were deliberately
NOT adopted, to keep scripts faithful to shipped figures. Sub-8-pt element sizes
inherited from the frozen renders remain as shipped; a future figure revision round
should apply the `pointsize=12 + par(cex=1)` convention repo-wide.

## 3. Integrated files (round-2 class B/C)

- 20 files renamed to the standard scheme (see `rename_log_v2.tsv`); `33-d4…` merged
  into `33_d4…`; `analyze_nira_prot.py`/`analyze_ruca_prot.py` superseded by the
  merged `s17_protonated_dbe.py`.
- 31 workspace-only files integrated as `scripts/s11–s51`, `code/56_toc_graphic.py`
  and `code/pipeline/p09_extract_veli_prot_cv.py` (7-field headers added; map in
  `rename_log_v2.tsv`).
- `md/autodl_setup/` — curated deployment archive (21 English-ified files: system
  builders incl. tleap inputs, per-system deploy scripts for the protonated reruns,
  transfer helpers, environment specs). 22 superseded/one-off items excluded —
  equivalents live under `code/pipeline/` (p01–p09) or are archived off-repo; list
  in the working archive `excluded_autodl/`. SSH hosts/ports are not stored in the
  repository (aliases only).

## 4. Data-layer changes (values refreshed to the audited R6/R9 state)

| artifact | change |
|---|---|
| `data/01_curated/extension_s1_wells.csv` | values updated to the r6f-audited recompute: fluzoparib 3.01/8.59/22.77, pamiparib 4.96/7.99/16.26, senaparib 3.24/6.45/12.50 (C1/C2/C3; sample SDs 10.2/5.9/4.7 reproduce the shipped figure's error bars) |
| `data/01_curated/helix_rmsf_table_s9.csv` | updated to the rerun_202609 recomputation (niraparib 6.94/2.56/4.84/5.11/3.20; olaparib 5.36/…; rucaparib 4.49/…) = SI Table S9 values |
| `data/analysis/Fig_Mechanism_Data.csv` | moved from `data/` root; refreshed from the R6-era render output |
| `data/analysis/r6/` | new: r6c_estimatrix, r6_panel_values, r6c_window_timescale, r6c_sensitivity, r6c_trunc_scan, r6d_ordinal_sensitivity, r6c_s1_seeds_canonical, r6c_s1/s2(_cv2)_spanfix, r6f_extension_summary.txt, `series/` (1.2 MB canonical-window dumps for the estimator matrix) |
| `data/analysis/pmf_features_summary.csv`, `S1_S2_ratios.csv` | refreshed from the R6-era outputs |

Data files are shipped as archived inputs; their internal content is not rewritten
(the per-frame `*.csv.gz` archives retain their original column notes).

## 5. Paths / hardening

- No absolute paths in shipped code: `$DATA_ROOT`-based reads for trajectories
  (s37/s38/s39/s48/s24 de-absolutised; s38 audit dumps → `data/analysis/r6/series/`);
  `/tmp` uses replaced by `mktemp`/`tempfile` (s22/s50).
- All console output and comments English (project-wide scan = 0 CJK in shipped
  code/scripts/common/md; data archives excluded by design).
- Full syntax pass: R `parse`, Python `py_compile`, `bash -n` — all clean (see §6).

## 6. Verification log (executed 2026-10-11, isolated copies)

Environment: `/opt/anaconda3` python 3.12 (MDAnalysis/OpenMM/numpy), TinyTeX-R
cairo devices, `DATA_ROOT=/Volumes/tjogzt4T/PARPi_data`. Deliverables were never
touched: all runs executed in `review-stage/release_v2/verify_run{,_2}` copies.

**6.1 Syntax** — R `parse`, Python `py_compile`, `bash -n` across `code/`,
`code/pipeline/`, `scripts/`, `md/autodl_setup/`: 0 failures.

**6.2 Full chain — `bash scripts/rebuild_all.sh`** (5-step version): PASS.
1. s02 seed dataset — 8/8 InChIKey match.
2. s03 trajectory metadata — 15 archival systems.
3. s08 consolidated S2 DBE + AAI table — **golden check passed (33 rows identical
   to the frozen table)**.
4. Figures 13/21/12/16 — all render; figure content matches the frozen PDFs at
   the text layer (see 6.4 for the rasterization note).
5. s10 manifest verification — **23/23 checksums OK** (4 checksums re-synced to
   the current frozen files: association_input_table, Fig_PCA_Features [text-layer
   convention], s2_window_timescale, p0_2_trunc_scan).

**6.3 Spot regeneration (isolated copies vs frozen deliverables)**

| target | exit | text layer | raster note |
|---|---|---|---|
| Fig_SI_Extension_Convergence (31) | 0 | strings identical (order-level diffs only) | cairo 1.17.6→1.18.4 |
| Fig_SI_Helix_RMSF_Heatmap (34) | 0 | identical | 2.7% px (cairo) |
| Fig_DCCM_Allostery (38) | 0 | identical strings (emission order) | cairo rasterization |
| Fig_S2_CV1_Retention (39) | 0 | identical strings | cairo rasterization |
| Fig_SI_Sampling_Toolkit (40) | 0 | identical strings incl. p-values | cairo rasterization |
| Fig_SI_ToyBenchmark (52) | 0 | identical strings | cairo rasterization |
| Fig3_State_Matrix (53) / Fig4_Window_Production (54) | 0 | identical strings | cairo rasterization |
| r6c_estimatrix.csv (s38) | 0 | **cell-for-cell identical (65 rows, max |Δ| = 0)** | — |

**6.4 Renderer-version note.** The frozen PDFs were produced with cairo **1.17.6**;
the current environment ships cairo **1.18.4**. Re-renders are content-identical
(identical text/strings/page sizes/fonts) but differ at the glyph-rasterization
and text-emission-order level. Byte-identical PDF reproduction therefore requires
pinning the original renderer stack; content-level reproduction is verified.

**6.5 Fixes landed during batch F**
- 31: convergence-CSV path corrected to `data/extension_s1_convergence.csv`.
- 29/33: `dir.create("results/figures")` added (fresh-clone byproduct writes).
- s38: correct canonical-window dumps archived under `data/analysis/r6/series/`
  (r6c2/r6d_cv2/r6s1c/r6s1/veli_prot npz, 30 files + full audit set); spot run
  reproduces the frozen matrix exactly. `veli_prot/*.json` summaries archived.
- Inputs refreshed to the R9 frozen-render state: `p0_2_trunc_scan.csv`,
  `s2_window_timescale.csv`.
- `s08` rewritten (v3): reads the corrected-span data layer (r6c spanfix +
  panel + R9 veli summary); carries updated frozen-state asserts and a
  golden-file check; changelog updated.
- `rebuild_all.sh` rewritten to the post-R6 5-step chain; s04–s07/s09 marked
  superseded (provenance retained).

**6.6 Residual notes**
- Some frozen SI figures carry sub-8-pt element sizes (inherited, approved);
  a future figure-revision round should apply the pointsize=12 + cex-pin
  convention repo-wide.
- Data archives (per-frame `*.csv.gz`) retain their original non-English column
  notes; code/scripts/common/docs are fully English (0 CJK — mechanical scan).
