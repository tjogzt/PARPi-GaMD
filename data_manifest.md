# data_manifest.md — manuscript-number provenance

Note: in this repository the shipped artifacts live under `data/`
(the working tree stores them under `results/analysis` and `results/figures`).

Every number cited in the manuscript traces back to a generating script through
this table. Columns: artifact, generating script (and the exact line that writes
it), upstream data, downstream consumers, manuscript location, protocol tag,
generation date, status, and MD5 checksum of the current file.

| artifact | generating_script | script_line | upstream_data | downstream | manuscript_location | protocol_tag | generation_date | status | checksum |
|---|---|---|---|---|---|---|---|---|---|
| results/figures/Fig_Mechanism_Data.csv | code/13-mechanism_figure.R | 238 | results/analysis/sys1_*_pmf_c3.xvg; results/analysis/pmf-c3-sys2_*_CV2_cv.dat.xvg | Fig_Mechanism_Master.pdf, Fig_2D_Mechanism_Map.pdf | Table 1 (tab:data): S1 well depths & AAI | C3 cumulant; S1 = DBE, S2 = DFW | 2026-09-15 | current | 47f1f8bf8aff3f9eec3948c23f150b46 |
| results/figures/Fig4B_spearman_data.csv | code/21-spearman_correlation.R | 114 | data/01_curated/trapping_potency.csv + Fig_Mechanism_Data.csv | Fig_Trapping_vs_Allostery.pdf | Results "S1 well depth and AAI show a suggestive inverse relationship"; Fig. fig:trapping | exact permutation test, n = 5 | 2026-09-16 | current | 2dd78ca166fd9dbe54fc20cc88683692 |
| results/analysis/replicate_well_depths.csv | code/18-replicate_analysis.R | 205 | results/replicates/rep_{tala,veli,azd}_*/pmf.npy | Fig_Replicate_Validation.pdf | Table (tab:replicates): 23.7/7.1/6.8 | histogram-reweighted CA CV | 2026-09-16 | current | 8c0ff8840a40c0f4cb58cb22223347d2 |
| results/analysis/cumulant_wells_C3.csv | code/extract_cv1_wells.py | 48 | results/analysis/pmf-c3-sys2_*_CV{1,2}_cv.dat.xvg | AAI denominators; Table S5 | Methods AAI definition; Table S5 | C3 cumulant | 2026-09-16 | current | e8d3edefe45d9196cdbcd73c8ddfd3f8 |
| results/analysis/rmsf_recomp/s2_rmsf_dccm_uniform.csv | code/recompute_all_s2.py | 114 | $DATA_ROOT sys2 DCDs (see script) | RMSF/DCCM figures (regen_*.R) | Results RMSF section; Table (tab:rmsf) | domain-aligned RMSF + whole-protein DCCM | 2026-09-15 | current | b3df05e44d267af9348b8005d6e134ef |
| results/analysis/extension_s1_convergence.csv | code/gen_extension_convergence_csv.py | 32 | results/analysis/new_drugs/sys1_*_cv.dat + weights | Fig_SI_Extension_Convergence.pdf | SI convergence (extension panel) | C3 cumulative blocks, 20 ns | 2026-09-16 | current | b31794a0217dc79727a43a45285bb866 |
| results/figures/S1_PMF_c3_stats.csv | code/09-s1_pmf_figures.R | 179 | results/analysis/sys1_*_pmf_c3.xvg | S1 PMF figures | Table 1 S1 wells (C3 values) | C3 cumulant | 2026-09-16 | current | 462445e2442fbb219a3cbdd5753597ed |
| data/01_curated/trapping_potency.csv | curated input (Murai 2012/2014; Zandarashvili 2020) | — | literature (PMIDs in file) | code/21-spearman_correlation.R | Table 1 trapping column (x olaparib) | western-blot/chromatin trapping, x olaparib | 2026-09-16 | current | eabd57ea712d41d7499428991c8f129e |
| data/01_curated/seed_ligands.csv | curated input (PubChem verified) | — | PubChem (CIDs in file) | scripts/build_seed_dataset.py | ML seed source | PubChem canonical SMILES/InChIKey | 2026-09-16 | current | 32c83af459eb7ecf49edf62c74d5aadb |
| data/01_curated/trapping_seed_dataset.csv | scripts/build_seed_dataset.py | 111 | data/01_curated/seed_ligands.csv | ML pipeline | ML seed dataset | RDKit descriptors (version-dependent) | 2026-09-16 | current | 5e19f5c328bf1c6bdf6708e939eefbe4 |
| results/archive/AAI_with_uncertainty.csv | superseded | — | old S2 CV2 values | — | — | superseded (see results/archive/README.md) | 2026-09-16 | superseded | see archive README |
| results/archive/S2_RMSF_all_systems.csv | superseded | — | pre-unification RMSF | — | — | superseded | 2026-09-16 | superseded | see archive README |
| results/archive/DCCM_HD_ART_stats.csv | superseded | — | pre-superposition DCCM | — | — | superseded | 2026-09-16 | superseded | see archive README |
| results/archive/S1_S2_ratios.csv | superseded | — | pre-rebuild S2 denominators | — | — | superseded | 2026-09-16 | superseded | see archive README |
| results/archive/S2_CV1_retention.csv | superseded | — | CV1/CV2 without CSV provenance | — | — | superseded (replaced by cumulant_wells_C3.csv) | 2026-09-16 | superseded | see archive README |
| results/archive/Fig4B_spearman_data.csv | superseded | — | included AZD5305 estimate row (n = 6) | — | — | superseded | 2026-09-16 | superseded | see archive README |

Notes:
- Checksums are MD5 of the file content at generation time; regenerate with `md5 -q <file>`.
- `script_line` refers to the line of the `write.csv`/CSV-writer statement in the current script.
- The extension-panel AAI values (fluzoparib 1.92, pamiparib 1.74, senaparib 1.86) are
  computed by code/s2_new_drugs_aai.py (console output; C3-based central value, SD from the
  C1-C3 spread) and are stored in Fig_Mechanism_Data.csv (wd_ratio column).
