# PARPi-GaMD — dual-system GaMD analysis of PARP1 inhibitor allosteric encoding

Analysis pipeline for the dual-system Gaussian accelerated molecular dynamics
(GaMD) study of ten PARP1 inhibitors. System S1 (catalytic domain only) isolates
the allosteric response of the helical domain (HD); system S2 (full PARP1-DNA
complex) captures the trapping environment. The pipeline quantifies HD-ART well
depths, the Allosteric Amplification Index (AAI), residue-level RMSF, dynamic
cross-correlations (DCCM), PCA landscapes, and replicate reproducibility.

## 1. Repository layout

| Directory | Role |
|---|---|
| `code/` | Analysis chain, numbered `NN_` (`01`–`48`): extraction, PMF reweighting, statistics, figures |
| `code/pipeline/` | GaMD deployment layer, numbered `pNN_`: system building (tleap), equilibration, config generation, CV extraction |
| `common/` | Shared configuration and helpers: `paths.py` (repository paths + `DATA_ROOT`), `pmf.py` (PyReweighting), `kabsch.py`, `helpers.R` (figure theme), `ligands.csv` (panel metadata, single source), `palette_si.csv` (SI figure palette) |
| `data/` | The archived data layer: `01_curated/` (literature/curated inputs), `analysis/` (analysis inputs: CV series, PMF curves, derived tables), plus replicate, per-frame and weight archives |
| `figures/` | The 27 in-manuscript figures: `pdf/` (vector) + `png/` (150-dpi previews) |
| `scripts/` | Build chain, numbered `sNN_` (`s01`–`s10`) + `rebuild_all.sh`: seed dataset, DBE reweighting chain, metadata, validation |
| `tools/` | Bundled PyReweighting tools |
| `docs/` | `rename_log.tsv` (full old→new name mapping), `standardisation_report.md` |
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

## 3. Citation and system requirements

If you use this code or data, please cite the corresponding manuscript
(submission details provided at publication) and the underlying software:
Amber (GaMD), MDAnalysis, PyReweighting, R/ggplot2.

- Python ≥ 3.10 with the packages in `requirements.txt`
- R ≥ 4.0 with the packages listed in `requirements.txt` (R section)
- AmberTools (cpptraj) for H-bond analysis (`CPPTRAJ` env var or PATH)
- PyReweighting (bundled under `tools/PyReweighting/`)
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
- Raw trajectories are NOT distributed. Set the environment variable
  `DATA_ROOT` to the trajectory data directory before running extraction or
  RMSF/DCCM scripts:
  ```bash
  export DATA_ROOT=/path/to/PARPi_data
  ```

## 6. Pipeline quickstart (all commands from the repository root)

1. Build chain — `scripts/` in order (or `bash scripts/rebuild_all.sh`):
   `s01_prepare_ligands.py` → `s02_build_seed_dataset.py` →
   `s03_build_trajectory_metadata.py`; the S2 DBE chain runs
   `s04`→`s05`→`s06`→`s07`→`s09`, consolidated by `s08_consolidate_s2_dbe.py`;
   `s10_verify_manifest.py` validates the data layer.
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

## 8. Key output map

| Output | Script | Description |
|---|---|---|
| `data/analysis/Fig_Mechanism_Data.csv` | `code/13_mechanism_figure.R` | S1/S2 well depths and AAI for all 10 systems |
| `data/analysis/association_input_table.csv` | `code/21_spearman_correlation.R` | n = 5 correlation input table (no AZD5305 row) |
| `data/replicate_well_depths.csv` | `code/18_replicate_analysis.R` | per-replicate histogram well depths |
| `data/cumulant_wells_C3.csv` | `code/35_extract_cv1_wells.py` | S2 CV1/CV2 C3-cumulant well depths |
| `data/analysis/rmsf_recomp/s2_rmsf_dccm_uniform.csv` | `code/36_recompute_s2_profiles.py` | domain-aligned RMSF + whole-protein DCCM summary |
| `figures/pdf/*.pdf` + `figures/png/*.png` | `code/*.R`, `code/50_fig1_structure.py` | the 27 in-manuscript figure files |

Superseded or one-off scripts are not distributed in this repository; the
archived working tree keeps them for provenance.
