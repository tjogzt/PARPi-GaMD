# workflow.md — end-to-end reproduction workflow

All commands run from the repository root. See `README.md` for layout and
`data_manifest.md` for the artifact→script provenance map.

## 0. Prerequisites

- Python ≥ 3.10 (`pip install -r requirements.txt`), R ≥ 4.0 with ggplot2/patchwork/
  data.table/etc., AmberTools (cpptraj) on PATH.
- Raw trajectories: set `export DATA_ROOT=/path/to/PARPi_data` (needed only by
  trajectory-level extraction/RMSF/DCCM scripts; all figure scripts below run
  purely from the `data/` layer).
- Bundled tools: `tools/PyReweighting/` (reweighting), `tools/gamd-openmm/` (engine).

## 1. Data-layer validation

```bash
python3 scripts/s10_verify_manifest.py     # checksums + row counts of the data layer
```

## 2. Build / recompute chain (scripts/)

```bash
bash scripts/rebuild_all.sh                # s02/s03 + s08 golden consolidation + figures + s10 manifest
# R6–R9 audit recomputations (canonical-window tables):
python3 scripts/s44_window_timescale.py    # → results/analysis/... (see script header)
python3 scripts/s45_timescale_scaling.py
python3 scripts/s38_estimator_matrix.py    # → r6c_estimatrix.csv (uses data/analysis/r6/series/)
python3 scripts/s37_bayesian_sensitivity_extended.py
python3 scripts/s39_convergence_recompute.py   # → r6f extension summary + convergence csv
python3 scripts/s48_toy_benchmark.py       # synthetic C7 benchmark
```

## 3. Figure regeneration (code/)

```bash
# main chain (shipped figures reproduce the frozen PDFs in figures/pdf/)
Rscript code/09_s1_pmf_figures.R
Rscript code/10_s2_pmf_figures.R
Rscript code/11_s1_s2_combined_pmf.R
Rscript code/12_descriptive_analysis.R
Rscript code/13_mechanism_figure.R
Rscript code/15_fig1_domain_arch.R
Rscript code/16_2d_landscape.R
Rscript code/18_replicate_analysis.R
Rscript code/21_spearman_correlation.R
Rscript code/29_d4_consistency_matrix.R
Rscript code/31_extension_convergence_figure.R
Rscript code/33_d4_diagnostics_figures.R
Rscript code/34_helix_rmsf_heatmap.R
Rscript code/38_dccm_figure.R
Rscript code/39_s2_cv1_retention_figure.R
Rscript code/40_sampling_toolkit_figure.R
# composite figures (consume the panels above)
Rscript code/52_fig_c7_benchmark.R
Rscript code/53_fig3_state_matrix.R
Rscript code/54_fig4_window_production.R
python3 code/55_fig6_evidence_summary.py
python3 code/56_toc_graphic.py
```

## 4. Statistics / audit outputs

```bash
Rscript code/28_pca_features.R             # PCA feature matrix (Fig S10)
python3 scripts/s17_protonated_dbe.py --lig nira   # protonated-state DBE (and --lig ruca)
python3 scripts/s20_eb47_primary_dbe.py
python3 scripts/s21_new_replicate_dbe.py --tag <TAG>
```
All scripts carry hard asserts (`stopifnot`) pinning the manuscript values they
feed; a passing run is the verification.

## 5. Deployment (GPU instances)

`code/pipeline/` (p01–p09) for the standardised build/deploy chain; the archived
record of the executions is in `md/autodl_setup/` (builders, per-system deploy
scripts, environment specs). Protocol: 300 K, lower-dual boost, seed 42, 50-ps
frames; protonated-state reruns share the identical protocol (see deploy headers).

## 6. Conventions recap

- Naming: `code/NN_`, `code/pipeline/pNN_`, `scripts/sNN_` (order encodes the chain).
- Header: 7 fields (Purpose / Created / Inputs / Outputs / Depends / Run).
- Inputs only from `data/`, `common/`, `$DATA_ROOT`; outputs: curated → `data/`
  in place, figures → `figures/`, regenerable intermediates → `results/`
  (git-ignored).
- No business constants in code — values live in `data/`/`common/*.csv` and are
  asserted, not typed inline. Seeds fixed (`set.seed(49)` / explicit seed args).
- The 2026-10-11 round-2 reconciliation policy and verification log:
  `docs/standardisation_report_v2.md`.
