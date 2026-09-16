# code/archive — superseded scripts (replacement chain)

Scripts moved here on 2026-09-16 by the pre-submission code audit. Nothing in
this directory is part of the current pipeline. Two kinds of files live here:
(1) the old `scripts/`-side copies that lost the dual-directory deduplication,
and (2) dead/superseded analysis scripts.

| Archived file | Why archived | Current authoritative path |
|---|---|---|
| `scripts_side_09-s1_pmf_figures.R` | Pre-extension version (7 systems); `code/09-s1_pmf_figures.R` adds the extension panel | `code/09-s1_pmf_figures.R` |
| `scripts_side_13-mechanism_figure.R` | Older AAI annotations ("Neutral (Type II)" / "Pro-release (Type III)"); `code/` version corrects the AAI semantics and adds the extension panel | `code/13-mechanism_figure.R` |
| `scripts_side_18-replicate_analysis.R` | Pre-audit version: unclosed `ifelse` parens (syntax error), wrong replicate counts, EB-47 no-data rows | `code/18-replicate_analysis.R` (fixed; reproduces Table 4) |
| `scripts_side_21-spearman_correlation.R` | Embedded trapping table incl. AZD5305 trapping = 100 estimate; n = 6 test; hard-coded conclusion text | `code/21-spearman_correlation.R` (reads curated CSVs; n = 5 exact permutation test) |
| `recompute_rmsf_dccm.py` | One-off S2 RMSF/DCCM rebuild for AZD5305/veliparib (superseded protocol) | `code/recompute_all_s2.py` (unified domain-aligned RMSF + whole-protein DCCM) |
| `s2_rerun_pmf.py` | One-off S2 DFW well-depth rebuild for AZD5305/veliparib | Outputs live in `results/analysis/new_drugs_s2/`; values consumed via `results/figures/Fig_Mechanism_Data.csv` |
| `extract_pmf_v3.py` | Hard-coded `/root/autodl-tmp` cloud-machine paths (not reproducible locally) | `code/05-sys1_cv_extract.py`, `code/gen_new_pmf_files.py` |
| `02-pmf_comparison.R` | Retired (2026-09-16 four-team review): output not cited in manuscript; superseded by 09/10 | — |
| `03-rmsd_rmsf.py` / `04-rmsd_rmsf_figure.R` | Retired: pre-unification RMSF; superseded by `recompute_all_s2.py` + `regen_rmsf_figure.R` | `code/recompute_all_s2.py` |
| `06-sys1_pmf_figure.R` | Retired: internal S1 CV figure (06-sys1_pmf_comparison.pdf not in manuscript); 09/10 are the canonical PMF figures | — |
| `07-master_figure.R` | Retired: superseded by 11-s1_s2_combined_pmf.R (Fig_Master_S1S2_HD_ART) | `code/11-s1_s2_combined_pmf.R` |
| `14-rmsf_diff_analysis.R` | Retired: superseded by domain-aligned recompute_all_s2.py + regen_rmsf_figure.R | `code/recompute_all_s2.py` |
| `17-toc_graphic.R` | Retired: TOC draft graphic, not in manuscript | — |
| `19-rmsf_all_systems.py` / `19-rmsf_type_comparison.R` | Retired: pre-superposition RMSF; superseded by recompute_all_s2.py | `code/recompute_all_s2.py` |
| `20-dccm_network.R` | Retired: DCCM network figure not in manuscript; whole-protein DCCM = recompute_all_s2.py + regen_dccm_figure.R | `code/regen_dccm_figure.R` |
| `22-cat_pocket_mapping.R` | Retired: pocket mapping not in manuscript | — |
| `23-allosteric_pathway.R` | Retired: pathway analysis not in manuscript | — |
| `24-transition_barriers.R` | Retired: transition barrier analysis not in manuscript (S1 barriers now from 12's pmf_features_summary.csv) | `code/12-descriptive_analysis.R` |
| `26-qsar_analysis.R` | Retired: internal QSAR correlation (Fig_QSAR not in manuscript). Refactored once to externalize hardcoded data, then retired; qsar_correlation.csv remains as an internal artifact | — |
| `27-hbond_protein_dna.R` | Retired: H-bond analysis not in manuscript | — |

Small python utility scripts that remain in `code/` are upstream pipeline steps
(trajectory prep, extraction) and are kept in place for reproducibility; see
data_manifest.md for the canonical producer of each cited artifact.
