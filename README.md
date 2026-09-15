# PARPi-GaMD — dual-system GaMD analysis of PARP1 inhibitor allosteric encoding

Analysis pipeline for the dual-system Gaussian accelerated molecular dynamics
(GaMD) study of ten PARP1 inhibitors. System S1 (catalytic domain only) isolates
the allosteric response of the helical domain (HD); system S2 (full PARP1-DNA
complex) captures the trapping environment. The pipeline quantifies HD-ART well
depths, the Allosteric Amplification Index (AAI), residue-level RMSF,
dynamic cross-correlations (DCCM), PCA landscapes, and replicate reproducibility.

## 1. Overview and output mapping

| Directory | Role |
|---|---|
| `code/` | Analysis chain: trajectory extraction, PMF reweighting, statistics, figures (numbered 01-30) |
| `code/archive/` | Superseded or non-reproducible scripts (replacement chain in its README) |
| `common/` | Shared helpers: `paths.py` (DATA_ROOT), `pmf.py` (PyReweighting), `kabsch.py` (alignment), `ligands.{csv,R,py}` (panel metadata, single source) |
| `data/` | Curated inputs (`01_curated/`) and replicate PMF data (`replicates/`) |
| `figures/` | Reproducible figure PDFs and PNG previews |
| `scripts/` (working tree only) | Build chain: PDB retrieval, ligand preparation, docking, seed-dataset construction |

Every manuscript number traces to its generating script via
[`data_manifest.md`](data_manifest.md) (artifact → script line → manuscript
location → checksum).

## 2. Citation

If you use this code or data, please cite the corresponding manuscript
(submission details provided at publication) and the underlying software:
Amber (GaMD), MDAnalysis, PyReweighting, R/ggplot2.

## 3. System requirements

- Python ≥ 3.10 with the packages in `requirements.txt`
- R ≥ 4.0 with the packages listed in `requirements.txt` (R section)
- AmberTools (cpptraj) for H-bond analysis (`CPPTRAJ` env var or PATH)
- PyReweighting (bundled under `tools/PyReweighting/` in the working tree)
- ~1 TB storage for the raw trajectories (not distributed in this repository)

## 4. Installation

```bash
git clone <this repository>
pip install -r requirements.txt        # Python dependencies
# R: install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "tibble",
#                       "bio3d", "igraph", "ggrepel"))
```

## 5. Datasets

- `data/replicates/` — histogram-reweighted replicate PMFs
  (`pmf.npy`/`rc.npy`/`dist_raw.npy`) for talazoparib, veliparib, AZD5305 and
  EB-47 (receptor 6VKK)
- `data/01_curated/` — trapping potencies (`trapping_potency.csv`), seed-ligand
  identities (`seed_ligands.csv`, PubChem-verified), derived descriptor dataset
  (`trapping_seed_dataset.csv`)
- Raw trajectories are NOT distributed. Set the environment variable
  `DATA_ROOT` to the trajectory data directory before running extraction or
  RMSF/DCCM scripts:
  ```bash
  export DATA_ROOT=/path/to/PARPi_data
  ```

## 6. Pipeline quickstart

1. `scripts/` (working tree): PDB retrieval → receptor/ligand preparation →
   docking → `build_seed_dataset.py` (curated seed dataset).
2. GaMD simulation and trajectory extraction on the compute cluster
   (`code/pipeline/` deploy scripts).
3. PMF reweighting: `code/gen_new_pmf_files.py`,
   `code/extract_cv1_wells.py` (C3 CV1/CV2 wells → `cumulant_wells_C3.csv`).
4. Figures and statistics (run from the repository root):
   ```bash
   Rscript code/13-mechanism_figure.R    # Fig_Mechanism_Data.csv (Table 1 chain)
   Rscript code/21-spearman_correlation.R  # n=5 exact permutation test
   Rscript code/18-replicate_analysis.R  # replicate Table
   Rscript code/09-s1_pmf_figures.R      # S1 PMF panels + stats
   ```

## 7. Regeneration scripts

Analysis scripts are numbered `01`-`30` in `code/`; the numbering is stable and
old versions are archived rather than renumbered. Scripts read inputs only from
`results/`, `data/`, `common/`, and `$DATA_ROOT`; they write to `results/`
(working tree). Deterministic outputs (no unseeded randomness; the ML pipeline
seeds at 49) mean re-running a script reproduces its outputs byte-for-byte —
verified for the manuscript-critical chains (see `data_manifest.md`).

## 8. Output description

| Output | Script | Description |
|---|---|---|
| `Fig_Mechanism_Data.csv` | `code/13-mechanism_figure.R` | S1/S2 well depths and AAI for all 10 systems |
| `Fig4B_spearman_data.csv` | `code/21-spearman_correlation.R` | n = 5 correlation data (no AZD5305 row) |
| `replicate_well_depths.csv` | `code/18-replicate_analysis.R` | per-replicate histogram well depths |
| `cumulant_wells_C3.csv` | `code/extract_cv1_wells.py` | S2 CV1/CV2 C3-cumulant well depths |
| `rmsf_recomp/s2_rmsf_dccm_uniform.csv` | `code/recompute_all_s2.py` | domain-aligned RMSF + whole-protein DCCM summary |
| figure PDFs | `code/0X-*.R`, `code/15-fig1_structure.py` | manuscript/SI figures |

Superseded intermediates were moved to `results/archive/` (working tree) with a
replacement-chain README; they are not part of this repository.
