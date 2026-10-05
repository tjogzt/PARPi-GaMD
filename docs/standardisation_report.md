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

## 7. Known limits

- Raw trajectories are not distributed; scripts that read them require
  `$DATA_ROOT` (documented per script header).
- `code/pipeline/*` and `scripts/s04–s09` run on the GPU host (conda env
  `gamd`); they are shipped for reproducibility of the simulation protocol, not
  for laptop execution.
- `README.md` and `data_manifest.md` were rewritten/updated to the new layout.
