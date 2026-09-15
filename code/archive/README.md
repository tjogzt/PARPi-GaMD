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
