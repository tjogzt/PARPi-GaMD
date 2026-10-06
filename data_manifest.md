# data_manifest.md — manuscript-number provenance

> **Repository standardisation (2026-10-05).** All paths and script names below
> use the standardised layout: `code/NN_*`, `code/pipeline/pNN_*`,
> `scripts/sNN_*`; analysis inputs live under `data/` (see `docs/rename_log.tsv`
> for the full mapping). Regenerable outputs are written to `results/`
> (git-ignored) or in place in `data/`. New curated tables added by the
> standardisation (provenance in `docs/standardisation_report.md`):
> `data/01_curated/two_state_panel_values.csv`,
> `data/01_curated/s2cv1_spans_table_s5.csv`,
> `data/01_curated/extension_s1_wells.csv`,
> `data/01_curated/helix_rmsf_table_s9.csv`,
> `data/01_curated/manuscript_anchor_values.csv`.

Every number cited in the manuscript traces back to a generating script through
this table. Columns: artifact, generating script (and the exact line that writes
it), upstream data, downstream consumers, manuscript location, protocol tag,
generation date, status, and MD5 checksum of the current file.

| artifact | generating_script | script_line | upstream_data | downstream | manuscript_location | protocol_tag | generation_date | status | checksum |
|---|---|---|---|---|---|---|---|---|---|
| data/analysis/Fig_Mechanism_Data.csv | code/13_mechanism_figure.R | 238 | data/analysis/sys1_*_pmf_c3.xvg; data/analysis/pmf-c3-sys2_*_CV2_cv.dat.xvg | Fig_Mechanism_Master.pdf, Fig_2D_Mechanism_Map.pdf | Table 1 (tab:data): S1 well depths & AAI | C3 cumulant; S1 = DBE, S2 = DBE (textbook boost-energy weights) | 2026-09-15 | current | TBD (regenerable) |
| data/analysis/association_input_table.csv | code/21_spearman_correlation.R | 114 | data/01_curated/trapping_potency.csv + Fig_Mechanism_Data.csv | Fig_Trapping_vs_Allostery.pdf | Results "S1 well depth and AAI show a suggestive inverse relationship"; Fig. fig:trapping | exact permutation test, n = 5 | 2026-09-16 | current | 09e3473cb20b87db35f6cdb2f80807d7 |
| data/replicate_well_depths.csv | code/18_replicate_analysis.R | 252 | results/replicates/rep_{tala,veli,azd,eb47}_*/pmf.npy | Fig_Replicate_Validation.pdf | Table (tab:replicates): 23.8/7.1/6.8/23.1 (EB-47: N = 2, 7AAB co-crystal ligand) | histogram-reweighted CA CV | 2026-10-04 | current | 543ac93538a9a6cc1c93079efa587e2f |
| results/analysis/trajectory_metadata.csv (regenerable; archived copy at data/trajectory_metadata.csv) | scripts/s03_build_trajectory_metadata.py | 1 | DATA_ROOT/*/{*.prmtop,*.dcd,gamd.log} | Methods "Note on replicates"; reviewer Section-8 package | per-system topology hash, frame counts, intervals, log fields | direct file audit | 2026-09-17 | current | TBD |
| data/cumulant_wells_C3.csv | code/35_extract_cv1_wells.py | 48 | data/analysis/pmf-c3-sys2_*_CV{1,2}_cv.dat.xvg | AAI denominators; Table S5 | Methods AAI definition; Table S5 | C3 cumulant | 2026-09-16 | current | 3e96eb42915944c7ffb4a4f85f35e26f |
| data/analysis/rmsf_recomp/s2_rmsf_dccm_uniform.csv | code/36_recompute_s2_profiles.py | 114 | $DATA_ROOT sys2 DCDs (see script) | RMSF/DCCM figures (regen_*.R) | Results RMSF section; Table (tab:rmsf) | domain-aligned RMSF + whole-protein DCCM (v2: nira/ola/ruca unified to rerun_202609) | 2026-10-06 | current | a651bb12e0e800079e75eb1414ac1a17 |
| data/extension_s1_convergence.csv | code/41_extension_panel_data_files.py | 32 | data/analysis/new_drugs/sys1_*_cv.dat + weights | Fig_SI_Extension_Convergence.pdf | SI convergence (extension panel) | C3 cumulative blocks, 20 ns | 2026-09-16 | current | b31794a0217dc79727a43a45285bb866 |
| data/analysis/S1_PMF_c3_stats.csv | code/09_s1_pmf_figures.R | 179 | data/analysis/sys1_*_pmf_c3.xvg | S1 PMF figures | Table 1 S1 wells (C3 values) | C3 cumulant | 2026-09-16 | current | TBD (regenerable) |
| data/01_curated/trapping_potency.csv | curated input (Murai 2012/2014; Zandarashvili 2020) | — | literature (PMIDs in file) | code/21_spearman_correlation.R | Table 1 trapping column (x olaparib) | western-blot/chromatin trapping, x olaparib | 2026-09-16 | current | 3c1d6afe0269aa75e4d5a8197c368aba |
| data/01_curated/seed_ligands.csv | curated input (PubChem verified) | — | PubChem (CIDs in file) | scripts/s02_build_seed_dataset.py | ML seed source | PubChem canonical SMILES/InChIKey | 2026-09-16 | current | 4fe9fc1c797bab73d7aab97761650986 |
| data/01_curated/trapping_seed_dataset.csv | scripts/s02_build_seed_dataset.py | 111 | data/01_curated/seed_ligands.csv | ML pipeline; 26-qsar_analysis.R (retired) descriptors | ML seed dataset | RDKit descriptors (version-dependent; incl. n_heavy/n_rings) | 2026-09-16 | current | 9cf9e0aed38f59d4b80349a5a6a57964 |
| data/01_curated/prior_sys1_cv.csv | curated input (prior S1 study) | — | early S1 trajectories (superseded; see (archived; not distributed) README.md) | 06-sys1_pmf_figure.R (not distributed) Panel D | internal figure 06-sys1_pmf_comparison.pdf | HD-ART CV mean±SD (old runs) | 2026-09-16 | current | f08315117753aea5862e82bbf8b12400 |
| data/qsar_correlation.csv | 26-qsar_analysis.R (retired; not distributed) | 109 | Fig_Mechanism_Data.csv + pmf_features_summary.csv + s2_rmsf_dccm_uniform.csv + pca_struct_system_stats.csv + seed/trapping curated CSVs | Fig_QSAR.pdf (internal, not manuscript-cited) | — | Spearman structure-dynamics correlations | 2026-09-16 | current | 8a04ca39a0e71897c5c5b1198c916409 |
| AAI_with_uncertainty.csv (archived; not distributed) | superseded | — | old S2 CV2 values | — | — | superseded (see (archived; not distributed) README.md) | 2026-09-16 | superseded | see archive README |
| S2_RMSF_all_systems.csv (archived; not distributed) | superseded | — | pre-unification RMSF | — | — | superseded | 2026-09-16 | superseded | see archive README |
| DCCM_HD_ART_stats.csv (archived; not distributed) | superseded | — | pre-superposition DCCM | — | — | superseded | 2026-09-16 | superseded | see archive README |
| S1_S2_ratios.csv (archived; not distributed) | superseded | — | pre-rebuild S2 denominators | — | — | superseded | 2026-09-16 | superseded | see archive README |
| S2_CV1_retention.csv (archived; not distributed) | superseded | — | CV1/CV2 without CSV provenance | — | — | superseded (replaced by cumulant_wells_C3.csv) | 2026-09-16 | superseded | see archive README |
| the superseded Fig4B table (archived; not distributed) | superseded | — | included AZD5305 estimate row (n = 6) | — | — | superseded | 2026-09-16 | superseded | see archive README |
| data/pca_projections.csv | code/28_pca_features.R (feature PCA) | 60 | sys1_*_pmf_c3.xvg + pmf-c3-sys2_*_CV{1,2}_cv.dat.xvg (24-dim feature matrix) | Fig_PCA_Features.pdf | Results PCA paragraph (Fig. S10) | prcomp, standardized 24-dim features | 2026-09-16 | current | (see notes) |
| data/pca_eigenvalues.csv | code/28_pca_features.R (feature PCA) | 55 | (same feature matrix) | Fig_PCA_Features.pdf axis labels | Fig. S10 | prcomp | 2026-09-16 | current | (see notes) |
| data/pca_system_stats.csv | code/28_pca_features.R (feature PCA) | 61 | pca_projections.csv | Fig_PCA_Features.pdf labels | Fig. S10 | per-system PC coordinates | 2026-09-16 | current | (see notes) |
| data/analysis/pca_struct_*.csv | 25-pca_analysis.R (structural PCA; not distributed) | 112-218 | $DATA_ROOT md_analysis CA DCDs | internal (26 retired) | — | bio3d pca.xyz; not manuscript-cited | 2026-09-16 | current | (see notes) |
| figures/pdf/Fig_PCA_Features.pdf | code/28_pca_features.R | 40 | pca_projections.csv + pca_eigenvalues.csv | SI Fig. S10 | SI Fig. S10 | 7-point PCA, China palette, Type II/III labels, seed=49 | 2026-10-05 | current | 8517febde3dfd6578e928d09ee179bcf |
| data/pmf_convergence.csv | LEGACY — original script lost | — | replicate GaMD runs (pre-npy era) | Fig_SI_PMF_Convergence.pdf (legacy; not cited) | — (legacy) | cumulative local-barrier well depth; exact recipe not recoverable from surviving artifacts; terminal depths validated by pmf.npy/Table 4 chain | 2026-09-16 | legacy | 7de25addc4a4ec9d7ac2e9bf16951197 |
| data/cumulant_convergence.csv | code/32_cumulant_convergence.R (reconstructed producer) | 75 | data/analysis/sys1_*_pmf_c{1,2,3}.xvg | 13-mechanism_figure.R (AAI SD asserts), 18-replicate_analysis.R (Table 1 comparison) | Table 1 S1 ±SD column; Table S1 | C1-C3 cumulant spread; pmf_cor_C2C3 recomputed 2026-09-16 (Table S1 r column updated) | 2026-09-16 | current | (see notes) |
| data/analysis/Fig_SI_PMF_Convergence.pdf (auxiliary; not in SI) | code/30_convergence_figure.R | 52 | data/pmf_convergence.csv (legacy) | — (legacy; not cited) | — (legacy) | thin replicate lines + LOESS + SEM ribbon + 80 ns threshold | 2026-09-16 | current | TBD (regenerable) |

Notes:
- Checksums are MD5 of the file content at generation time; regenerate with `md5 -q <file>`.
- `script_line` refers to the line of the `write.csv`/CSV-writer statement in the current script.
- The extension-panel AAI values reported in the manuscript (fluzoparib 1.15, senaparib 1.06,
  pamiparib 0.94) are computed by the textbook-DBE S2 pipeline
  (scripts/s07_regenerate_s2_pmf_xvgs.py + scripts/s08_consolidate_s2_dbe.py, C3-based central value,
  SD from the C1-C3 spread) and are stored in data/analysis/s2_dbe_final.csv and
  Fig_Mechanism_Data.csv (wd_ratio column). The earlier the former s2_new_drugs_aai.py (retired, not distributed) values (1.92/1.74/1.86,
  force-scaled S2 scheme) are superseded; that script is marked LEGACY.
- AAI SD (Table 1 "AAI ±" column): canonical propagation aai_sd = wd_sd(C1-C3) / wd_S2,
  implemented as hard asserts in code/13_mechanism_figure.R (2026-09-16 recomputation;
  the previous manuscript column was partially based on pre-rebuild S2 denominators).
- n_states semantics: the "> 0.5 kcal depth" filter in extract_features
  (12_descriptive_analysis.R / 28_pca_features.R / 08_ml_pipeline.R) is a tautology
  (x <= x + 0.5) and never fires; n_states = count of ALL local minima. Values unchanged;
  comments updated 2026-09-16.
- xvg files have variable header lengths (typically 5 comment lines); scripts must use the
  comment-aware read_pmf in common/helpers.R — a fixed skip (e.g. skip=10) silently drops
  data rows (this corrupted the pre-2026-09-16 PCA feature matrix; fixed and re-derived).
- `code/10_s2_pmf_figures.R` additionally emits `figures/pdf/Fig_S2_pmf_overlay_combined.pdf`
  (a legacy single-file overlay variant, superseded by the split protein–DNA / HD–ART
  overlays). It is not referenced by any manuscript text and is intentionally excluded from
  this tree (the workspace keeps a historical copy; re-running the script re-creates it).
- Figure files were refreshed 2026-10-05 by the user-approved figure-review round
  (see docs/standardisation_report.md §9); checksums above reflect the current files.
- data/analysis/pca_*.csv are produced by code/28_pca_features.R (24-dim PMF feature PCA,
  manuscript-cited); pca_struct_*.csv by 25-pca_analysis.R (not distributed) (structural PCA, internal).