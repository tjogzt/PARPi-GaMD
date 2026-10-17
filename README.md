# PARPi-GaMD — dual-system GaMD analysis of PARP1 inhibitor allosteric encoding

Analysis pipeline for the dual-system Gaussian accelerated molecular dynamics
(GaMD) study of ten PARP1 inhibitors. System S1 (catalytic domain only) isolates
the allosteric response of the helical domain (HD); system S2 (full PARP1-DNA
complex) captures the trapping environment. The pipeline quantifies HD–ART sampled-range
spans and the S1/S2 span ratio, and implements a six-diagnostic robustness audit
of descriptor–phenotype rank associations (chemical-state assignment, estimator
dependence, sampling-window calibration with a planted-association positive
control, cross-docking validation, and small-panel association stability),
alongside residue-level RMSF, dynamic cross-correlations (DCCM), and PCA
landscapes.

## 1. Repository layout

| Directory | Role |
|---|---|
| `code/` | Analysis chain, numbered `NN_` (`01`–`56`): extraction, PMF reweighting, statistics, figures |
| `code/pipeline/` | GaMD deployment layer, numbered `pNN_`: system building (tleap), equilibration, config generation, CV extraction |
| `common/` | Shared configuration and helpers: `paths.py` (repository paths + `DATA_ROOT`), `pmf.py` (PyReweighting), `kabsch.py`, `helpers.R` (figure theme), `ligands.csv` (panel metadata, single source), `palette_si.csv` (SI figure palette) |
| `data/` | The archived data layer: `01_curated/` (literature/curated inputs), `analysis/` (analysis inputs: CV series, PMF curves, derived tables), plus replicate, per-frame and weight archives |
| `figures/` | The 27 in-manuscript figures: `pdf/` (vector) + `png/` (150-dpi previews); `TOC_Graphic_jcim.*` (submission TOC) |
| `md/` | `autodl_setup/` — archived deployment record (builders incl. tleap inputs, per-system deploy scripts, environment specs) |
| `scripts/` | Build / statistics chain, numbered `sNN_` (`s01`–`s53`) + `rebuild_all.sh`: seed dataset, DBE reweighting chain, metadata, validation, R6–R14 audit recomputations (incl. positive control, cross-docking) |
| `tools/` | Bundled tools: PyReweighting + the vendored GaMD engine (`gamd-openmm/`; provenance in `tools/gamd-openmm/VENDOR_INFO.md`) |
| `docs/` | `rename_log.tsv` + `rename_log_v2.tsv` (old→new mappings, rounds 1–2), `standardisation_report.md` + `standardisation_report_v2.md`, `workflow.md` (end-to-end reproduction), `comment_conventions.md`, `determinism_notes.md` |
| `results/` | Regenerable working outputs (not tracked; created on demand) |

Every manuscript number traces to its generating script via
[`data_manifest.md`](data_manifest.md) (artifact → script → manuscript location →
checksum).

## 2. Naming and layout conventions (repository standardisation, 2026-10-05)

- All code files are named `<number>_<snake_case_description>.<ext>`; the number
  encodes the pipeline order within its chain (`code/` = analysis, `code/pipeline/`
  = deployment, `scripts/` = build). The complete old→new mapping is in
  `docs/rename_log.tsv`.
- Every code file carries a standard English header: Purpose, Author, Created,
  Inputs, Outputs, Depends, Run.
- No business constants are hard-coded in analysis logic: all panel values,
  anchor assertions and style palettes are read from the `data/` layer or
  `common/*.csv` (see `docs/standardisation_report.md` for the value→file map).
- Scripts read inputs only from `data/`, `common/`, and `$DATA_ROOT`, and write
  either (a) curated outputs in place in `data/`, (b) figure PDFs/PNGs into
  `figures/`, or (c) regenerable intermediates into `results/` (git-ignored).
- Round 2 (2026-10-11) integrated the R6–R10 work: chains extended to the ranges
  above, audit tables added under `data/analysis/r6/`, and eleven dual-edited files
  reconciled by keeping the revision that generated the shipped (frozen) figures —
  policy and evidence in `docs/standardisation_report_v2.md`. Round 3 (2026-10-16) integrated the R12–R14 audit revisions: terminology unified throughout (sampled-range span; S1/S2 span ratio; audit-corrected), the positive-control benchmark (`scripts/s52_positive_control.py`) and the cross-docking RMSD validation (`scripts/s53_crossdock_rmsd.py`) added, the rerun-based 2D PMF grids for the three re-simulated S2 systems restored to the figure data layer, and every affected figure re-rendered on this tree with figure text layers verified against the manuscript renders. Round 4 (2026-10-17) applies the fourth external review: the positive-control benchmark uses the panel's exact trapping-rank tie structure (olaparib--rucaparib midrank 2.5) with a Δ = 0 null calibration on a five-level σ grid; the estimator-matrix script writes a reproducible copy (`data/analysis/r6/r6c_estimatrix_repro.csv`) that matches the published table identically; the occupied-support and effective-support definitions are stated operationally; and the Fig. 1 simulation-length annotation, the Fig. 6 diagnostic labels, and Fig. S13 were re-rendered accordingly. Round 5 (2026-10-17): the chemical-state panel table (Table~S8) now covers the complete eight-combination lattice (`data/analysis/r6/r9_state_lattice.csv`, `scripts/s54_state_combinations.py`; 32 cells, no nominally significant cell, primary S2 CV1 estimate negative in all eight combinations).

## 3. Citation and system requirements

If you use this code or data, please cite the corresponding manuscript
(submission details provided at publication) and the underlying software:
Amber (GaMD), MDAnalysis, PyReweighting, R/ggplot2.

- Python ≥ 3.10 with the packages in `requirements.txt`
- R ≥ 4.0 with the packages listed in `requirements.txt` (R section)
- AmberTools (cpptraj) for H-bond analysis (`CPPTRAJ` env var or PATH)
- PyReweighting (bundled under `tools/PyReweighting/`)
- GaMD engine: gamd-openmm (bundled under `tools/gamd-openmm/`; run via `PYTHONPATH` or `pip install`; MIT)
- ~1 TB storage for the raw trajectories (not distributed in this repository)

## 4. Installation

```bash
git clone <this repository>
pip install -r requirements.txt        # Python dependencies
# R: install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "tibble",
#                       "data.table", "bio3d", "igraph", "ggrepel", "reticulate"))
```

## 5. Datasets

- `data/replicates/` — histogram-reweighted replicate PMFs
  (`pmf.npy`/`rc.npy`/`dist_raw.npy`) for talazoparib, veliparib, AZD5305, and
  EB-47 (N = 2 replicates; co-crystal EB-47 ligand from PDB 7AAB on the 6VKK
  CAT-domain receptor)
- `data/per_frame/` — archived per-frame boost energies and collective-variable
  values (`*.csv.gz`; columns: frame, step, dV_D, dV_T, cv) for the base,
  re-run, and protonated systems (the "archived per-frame boost energies"
  referenced in the manuscript)
- `data/analysis_weights/` — per-frame PyReweighting weight inputs
  (`*_weights.dat`) for the S1 and S2 CV1/CV2 analyses
- `data/01_curated/` — trapping potencies, seed-ligand identities, derived
  descriptor datasets, and the curated panel tables (`two_state_panel_values.csv`,
  `s2cv1_spans_table_s5.csv`, `extension_s1_wells.csv`, `helix_rmsf_table_s9.csv`,
  `manuscript_anchor_values.csv`)
- `data/analysis/` — the analysis-input layer: per-system CV series
  (`analysis_CV*.dat`), PMF-c3 curves (`.xvg`, `.npz`), derived tables
  (`s2_dbe_final.csv`, `s1_dbe_unified_wells.csv`, `s1_table1_rows.csv`,
  `d4_consistency_matrix.csv`, …), and the recomputed RMSF/DCCM profiles
- `data/01_modeling/docking_system1/` — cross-docking validation inputs: docked
  poses (SDF + PDBQT), ligand templates, and the five co-crystal receptor
  structures (PDB 6VKK/5DS3/4R6E/4UND/2RD6) behind
  `data/analysis/crossdock_rmsd_results.csv` (Table S16)
- Raw trajectories are NOT distributed. Set the environment variable
  `DATA_ROOT` to the trajectory data directory before running extraction or
  RMSF/DCCM scripts:
  ```bash
  export DATA_ROOT=/path/to/PARPi_data
  ```

## 6. Pipeline quickstart (all commands from the repository root)

1. Build chain — `scripts/` (or `bash scripts/rebuild_all.sh`):
   `s01_prepare_ligands.py` → `s02_build_seed_dataset.py` →
   `s03_build_trajectory_metadata.py`; `s08_consolidate_s2_dbe.py` consolidates
   and golden-checks the S2 DBE + span-ratio table against the frozen state;
   `s10_verify_manifest.py` validates the data layer; `s52_positive_control.py`
   and `s53_crossdock_rmsd.py` regenerate the positive-control benchmark
   (`data/analysis/r6/c7g_positive_control.csv`) and the cross-docking RMSD
   table (`data/analysis/crossdock_rmsd_results.csv`, Table S16). The pre-R6 generators
   (`s04`–`s07`, `s09`) are superseded (kept for provenance).
2. Deployment / simulations — `code/pipeline/` (`p01`–`p08`) on the GPU
   workstation (conda env `gamd`); see the scripts' headers for the protocol
   parameters (300 K, lower-dual boost, seed 42, 50 ps frames).
3. Derived data files for the extension panel —
   `python3 code/41_extension_panel_data_files.py`.
4. Figures and statistics, e.g.:
   ```bash
   Rscript code/13_mechanism_figure.R       # core mechanism figure + data table
   Rscript code/21_spearman_correlation.R   # n = 5 exact permutation test
   Rscript code/18_replicate_analysis.R     # replicate table
   Rscript code/09_s1_pmf_figures.R         # S1 PMF panels + stats
   Rscript code/31_extension_convergence_figure.R   # SI Figure S15
   Rscript code/34_helix_rmsf_heatmap.R     # SI Figure S16
   ```

## 7. Reproducibility

Deterministic outputs (the ML pipeline seeds at 49) mean re-running a script
reproduces its outputs. During the 2026-10-05 standardisation the
manuscript-critical figure chains were verified by re-running them on the
standardised tree: `code/09`, `code/16`, `code/21`, `code/28`, `code/31`,
`code/34` and `code/39` reproduce their shipped figures pixel-identical, and
`code/32` reproduces `data/cumulant_convergence.csv` with maximum cell
deviation < 1e-4 (asserts in the script). `Fig_D4_ConsistencyMatrix` and
`Fig_DCCM_Allostery` correspond to an earlier data revision and are kept as
shipped for manuscript consistency. See `docs/standardisation_report.md` for
the full verification log and the regeneration caveats for `code/41`/`code/44`.
Round 2 (2026-10-11): every integrated/reconciled script passes syntax checks and
the submission tree replays its figure chains — the round-2 verification log is
appended to `docs/standardisation_report_v2.md`. Round 3 (2026-10-16): every R12–R14-revised figure chain was re-rendered on this tree and its text layer verified figure-by-figure against the manuscript renders; `scripts/s52_positive_control.py` reproduces `data/analysis/r6/c7g_positive_control.csv` byte-identically (seed 49), and `scripts/s53_crossdock_rmsd.py` reproduces `data/analysis/crossdock_rmsd_results.csv` (0.58/2.46/0.91/1.45/1.32 Å).

## 8. Key output map

| Output | Script | Description |
|---|---|---|
| `data/analysis/Fig_Mechanism_Data.csv` | `code/13_mechanism_figure.R` | S1/S2 sampled-range spans and the S1/S2 span ratio for all 10 systems |
| `data/analysis/r6/c7g_positive_control.csv` | `scripts/s52_positive_control.py` | planted-association positive control (125 cells: 5 σ × 25 Δ; tie-aware reference ranks, Δ = 0 null calibration; 1,500 worlds/cell) |
| `data/analysis/crossdock_rmsd_results.csv` | `scripts/s53_crossdock_rmsd.py` | cross-docking RMSD vs five co-crystal structures (Table S16) |
| `data/analysis/association_input_table.csv` | `code/21_spearman_correlation.R` | n = 5 correlation input table (no AZD5305 row) |
| `data/replicate_well_depths.csv` | `code/18_replicate_analysis.R` | per-replicate histogram-reweighted spans |
| `data/cumulant_wells_C3.csv` | `code/35_extract_cv1_wells.py` | S2 CV1/CV2 C3-cumulant spans |
| `data/analysis/rmsf_recomp/s2_rmsf_dccm_uniform.csv` | `code/36_recompute_s2_profiles.py` | domain-aligned RMSF + whole-protein DCCM summary |
| `figures/pdf/*.pdf` + `figures/png/*.png` | `code/*.R`, `code/50_fig1_structure.py` | the 27 in-manuscript figure files |

| `data/analysis/r6/r6c_estimatrix.csv` | `scripts/s38_estimator_matrix.py` | estimator × metric consistency matrix (canonical-window dumps in `data/analysis/r6/series/`) |
| `data/analysis/r6/r6_panel_values.csv` | R6 audit recompute | audited panel values used by codes 12/13/28/33 |

Superseded or one-off scripts are not distributed in this repository; the
archived working tree keeps them for provenance.
